package varstore

import (
	"strconv"
	"strings"

	"github.com/compgenlab/cghts/vcf"
)

// This file holds the one authoritative reading of a VCF genotype. Both the
// converter and the VCF-backed Store go through it, so a query against a VCF
// and the same query against a Parquet store converted from it cannot drift.

// SplitGT recodes a raw GT string for a single ALT allele, given that allele's
// 1-based index in the record's ALT list. It returns the recoded biallelic
// genotype and whether the sample carries this particular allele.
//
// Alleles other than the focal one and the reference become "." rather than
// "0". A 1/2 sample genuinely carries both alt alleles, so it must appear as a
// carrier in each split row; calling the other allele reference would invent a
// reference observation that the data does not support.
// gtAlleles splits a raw GT into its alleles, normalizing the phased separator
// away. Seven places did this inline, in two files, character for character --
// including the two in vcfstore.go that this file's header says it exists to
// keep from drifting.
func gtAlleles(gt string) []string {
	return strings.Split(strings.ReplaceAll(gt, "|", "/"), "/")
}

func SplitGT(gt string, altIdx int) (recoded string, carrier bool) {
	if gt == "" {
		return ".", false
	}
	sep := "/"
	if strings.Contains(gt, "|") {
		sep = "|"
	}
	parts := gtAlleles(gt)
	out := make([]string, 0, len(parts))
	for _, a := range parts {
		switch {
		case a == "." || a == "":
			out = append(out, ".")
		case a == "0":
			out = append(out, "0")
		default:
			n, err := strconv.Atoi(a)
			if err != nil {
				out = append(out, ".")
				continue
			}
			if n == altIdx {
				out = append(out, "1")
				carrier = true
			} else {
				out = append(out, ".")
			}
		}
	}
	return strings.Join(out, sep), carrier
}

// IsAltCarrier reports whether a raw GT carries any non-reference allele.
func IsAltCarrier(gt string) bool {
	for _, a := range gtAlleles(gt) {
		if a != "" && a != "." && a != "0" {
			return true
		}
	}
	return false
}

// HasCall reports whether a raw GT carries at least one called allele.
//
// Depth alone does not make a site callable: a "./." genotype at DP 30 means
// the caller saw reads and still declined to call, which is not the positive
// observation a callable region is meant to assert.
func HasCall(gt string) bool {
	for _, a := range gtAlleles(gt) {
		if a != "" && a != "." {
			return true
		}
	}
	return false
}

// IsHomRef reports whether a raw GT is an explicit all-reference call. This is
// the positive observation that separates a non-carrier from a sample that was
// simply never assayed, so "./." and an empty GT deliberately do not qualify.
func IsHomRef(gt string) bool {
	seen := false
	for _, a := range gtAlleles(gt) {
		if a == "" {
			continue
		}
		if a != "0" {
			return false
		}
		seen = true
	}
	return seen
}

// SplitAD returns AD[0] and AD[altIdx] from a raw AD field, or Missing.
//
// Depth is taken per allele rather than summed across alt alleles: at a
// multiallelic site the depth supporting allele 1 says nothing about allele 2.
func SplitAD(ad string, altIdx int) (adRef, adAlt int32) {
	adRef, adAlt = Missing, Missing
	if ad == "" || ad == "." {
		return
	}
	f := strings.Split(ad, ",")
	if len(f) > 0 {
		adRef = atoi32(f[0])
	}
	if altIdx < len(f) {
		adAlt = atoi32(f[altIdx])
	}
	return
}

func atoi32(s string) int32 {
	n, err := strconv.Atoi(strings.TrimSpace(s))
	if err != nil {
		return Missing
	}
	return int32(n)
}

// AddAlleleCounts tallies one raw GT into allele counts, returning that
// genotype's contribution to AN and incrementing acc[k] once per copy of ALT
// allele k+1.
//
// These are allele counts, not sample counts: a 1/1 genotype contributes 2 to
// AC and 2 to AN where it would count as only one carrier. Missing alleles
// contribute to neither. acc is written in place so a whole cohort can be
// tallied without allocating per sample.
func AddAlleleCounts(gt string, acc []int32) (an int32) {
	for _, a := range gtAlleles(gt) {
		if a == "" || a == "." {
			continue
		}
		an++
		n, err := strconv.Atoi(a)
		if err != nil || n <= 0 { // 0 is the reference allele
			continue
		}
		if n-1 < len(acc) {
			acc[n-1]++
		}
	}
	return an
}

// SampleFields is the subset of a per-sample FORMAT column this package reads.
type SampleFields struct {
	GT string
	DP int32
	GQ int32
	AD string
}

// ReadSample pulls the GT/DP/GQ/AD fields for one sample of a record. Absent
// fields come back as Missing (or "" for GT/AD).
func ReadSample(rec *vcf.VcfRecord, i int) (SampleFields, error) {
	s := SampleFields{DP: Missing, GQ: Missing}
	attrs, err := rec.Sample(i)
	if err != nil {
		return s, err
	}
	if v, ok := attrs.Get("GT"); ok {
		s.GT = v.String()
	}
	if v, ok := attrs.Get("DP"); ok {
		s.DP = atoi32(v.String())
	}
	if v, ok := attrs.Get("GQ"); ok {
		s.GQ = atoi32(v.String())
	}
	if v, ok := attrs.Get("AD"); ok {
		s.AD = v.String()
	}
	return s, nil
}

// RecordLoci returns one Locus per ALT allele of a record, in ALT order.
//
// The chromosome keeps the record's own spelling, deliberately: a store records
// the naming its source used rather than rewriting it, and comparisons are made
// by canonical identity at query time instead (see CanonKey). This comment used
// to claim the name was normalized here, which is not what the code has ever
// done -- a trap for the converter that calls it.
func RecordLoci(rec *vcf.VcfRecord) []Locus {
	alts := rec.Alt()
	out := make([]Locus, 0, len(alts))
	for _, a := range alts {
		out = append(out, Locus{
			Chrom: rec.Chrom,
			Pos:   int32(rec.Pos),
			Ref:   rec.Ref,
			Alt:   a,
		})
	}
	return out
}

// CallFor builds the Call for one sample at one ALT allele of a record, or
// reports false when the sample does not carry that allele.
func CallFor(rec *vcf.VcfRecord, sample string, sf SampleFields, altIdx int, alt string) (Call, bool) {
	gt, carrier := SplitGT(sf.GT, altIdx)
	if !carrier {
		return Call{}, false
	}
	adRef, adAlt := SplitAD(sf.AD, altIdx)
	return Call{
		SampleID: sample,
		Chrom:    rec.Chrom,
		Pos:      int32(rec.Pos),
		Ref:      rec.Ref,
		Alt:      alt,
		// The record's reference extent, so a VCF-backed row carries the same
		// span a converted store's does. Leaving it zero meant every deletion,
		// symbolic ALT and long REF read back as one base wide from a VCF and
		// its true width from the store built out of that same VCF -- and the
		// two backends are supposed to answer identically.
		RefEnd: int32(rec.RefSpanEnd()),
		GT:     gt,
		DP:     sf.DP,
		ADRef:  adRef,
		ADAlt:  adAlt,
		GQ:     sf.GQ,
	}, true
}

// IsHomAlt reports whether a RECODED genotype is homozygous for the focal
// allele: every called allele is the focal one.
//
// Recoded, which matters -- a store's GT has already been rewritten per ALT, so
// the focal allele is 1, reference is 0, and any other alternate is ".". A 1/2
// sample therefore reads as "1/." on the focal allele's row.
func IsHomAlt(gt string) bool {
	alleles := gtAlleles(gt)
	if len(alleles) == 0 {
		return false
	}
	for _, a := range alleles {
		if a != "1" {
			return false
		}
	}
	return true
}

// IsHet reports whether a recoded genotype is a CLEAN heterozygote: exactly the
// focal allele and the reference.
//
// Deliberately excludes 1/. -- a sample carrying the focal allele and some other
// alternate. Its AD holds reads for the reference and the focal allele, and the
// third allele's reads are somewhere else or nowhere, so the alt fraction is not
// a quantity an allele-balance threshold can be applied to. Gating it on a
// number that does not describe it would discard real compound carriers.
func IsHet(gt string) bool {
	alleles := gtAlleles(gt)
	if len(alleles) != 2 {
		return false
	}
	var ref, focal int
	for _, a := range alleles {
		switch a {
		case "0":
			ref++
		case "1":
			focal++
		default:
			return false
		}
	}
	return ref == 1 && focal == 1
}
