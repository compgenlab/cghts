package gtstore

import (
	"fmt"
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
// forEachAllele calls fn with each allele of a raw GT, in order, splitting on
// either separator ("/" unphased, "|" phased) and stopping early when fn
// returns false. It yields what splitting on both separators would -- an empty
// GT is one empty allele, "0/" is "0" then "" -- as substrings of gt.
//
// It replaced gtAlleles, which built that list with ReplaceAll and Split. Every
// genotype passes through here three or four times during a conversion, so at
// tens of thousands of samples those two allocations per call were a large
// share of what the garbage collector spent its time on.
func forEachAllele(gt string, fn func(a string) bool) {
	start := 0
	for i := 0; i < len(gt); i++ {
		if gt[i] == '/' || gt[i] == '|' {
			if !fn(gt[start:i]) {
				return
			}
			start = i + 1
		}
	}
	fn(gt[start:])
}

// gtAlleles returns forEachAllele's alleles as a slice, for the callers that
// need to count them first; the per-genotype paths use forEachAllele directly.
func gtAlleles(gt string) []string {
	var out []string
	forEachAllele(gt, func(a string) bool { out = append(out, a); return true })
	return out
}

// alleleIndex parses one allele as a non-negative integer without allocating;
// ok is false for ".", "" and anything strconv.Atoi would reject.
func alleleIndex(a string) (n int, ok bool) {
	if a == "" {
		return 0, false
	}
	if len(a) > 9 { // could overflow the loop below; rare enough to defer
		v, err := strconv.Atoi(a)
		return v, err == nil
	}
	for i := 0; i < len(a); i++ {
		c := a[i]
		if c < '0' || c > '9' {
			v, err := strconv.Atoi(a) // signs and the like: defer to Atoi
			return v, err == nil
		}
		n = n*10 + int(c-'0')
	}
	return n, true
}

// carriesAllele reports whether a raw GT carries ALT allele altIdx (1-based),
// what SplitGT's second result says, without building the recoded genotype.
func carriesAllele(gt string, altIdx int) bool {
	found := false
	forEachAllele(gt, func(a string) bool {
		if a == "." || a == "" || a == "0" {
			return true
		}
		if n, ok := alleleIndex(a); ok && n == altIdx {
			found = true
			return false
		}
		return true
	})
	return found
}

func SplitGT(gt string, altIdx int) (recoded string, carrier bool) {
	if gt == "" {
		return ".", false
	}
	sep := byte('/')
	if strings.IndexByte(gt, '|') >= 0 {
		sep = '|'
	}
	var b strings.Builder
	b.Grow(len(gt))
	first := true
	forEachAllele(gt, func(a string) bool {
		if !first {
			b.WriteByte(sep)
		}
		first = false
		switch {
		case a == "." || a == "":
			b.WriteByte('.')
		case a == "0":
			b.WriteByte('0')
		default:
			n, ok := alleleIndex(a)
			if ok && n == altIdx {
				b.WriteByte('1')
				carrier = true
			} else {
				b.WriteByte('.')
			}
		}
		return true
	})
	return b.String(), carrier
}

// IsAltCarrier reports whether a raw GT carries any non-reference allele.
func IsAltCarrier(gt string) bool {
	found := false
	forEachAllele(gt, func(a string) bool {
		found = a != "" && a != "." && a != "0"
		return !found
	})
	return found
}

// HasCall reports whether a raw GT carries at least one called allele.
//
// Depth alone does not make a site callable: a "./." genotype at DP 30 means
// the caller saw reads and still declined to call, which is not the positive
// observation a callable region is meant to assert.
func HasCall(gt string) bool {
	found := false
	forEachAllele(gt, func(a string) bool {
		found = a != "" && a != "."
		return !found
	})
	return found
}

// IsHomRef reports whether a raw GT is an explicit all-reference call. This is
// the positive observation that separates a non-carrier from a sample that was
// simply never assayed, so "./." and an empty GT deliberately do not qualify.
func IsHomRef(gt string) bool {
	seen, other := false, false
	forEachAllele(gt, func(a string) bool {
		if a == "" {
			return true
		}
		if a != "0" {
			other = true
			return false
		}
		seen = true
		return true
	})
	return seen && !other
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
	forEachAllele(gt, func(a string) bool {
		if a == "" || a == "." {
			return true
		}
		an++
		n, ok := alleleIndex(a)
		if !ok || n <= 0 { // 0 is the reference allele
			return true
		}
		if n-1 < len(acc) {
			acc[n-1]++
		}
		return true
	})
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
//
// It resolves the FORMAT positions on every call; a loop over a record's
// samples should build one SampleReader and call Read instead.
func ReadSample(rec *vcf.VcfRecord, i int) (SampleFields, error) {
	if n := rec.NumSamples(); i < 0 || i >= n {
		return SampleFields{DP: Missing, GQ: Missing},
			fmt.Errorf("vcf: sample index %d out of range (%d samples)", i, n)
	}
	return NewSampleReader(rec).Read(i), nil
}

// SampleReader reads GT/DP/GQ/AD from one record's samples, resolving where
// those keys sit in the FORMAT column once rather than per sample.
//
// It reads each field straight out of the record line (vcf.SampleField), never
// building the per-sample map that vcf.Sample does. That map, built for every
// one of tens of thousands of samples and holding every FORMAT key -- PL, GP,
// F1R2 and the rest, of which this package reads four -- was most of the time a
// wide conversion spent, much of it in the garbage collector.
type SampleReader struct {
	rec            *vcf.VcfRecord
	gt, dp, gq, ad int
}

// NewSampleReader prepares to read rec's samples.
func NewSampleReader(rec *vcf.VcfRecord) SampleReader {
	return SampleReader{
		rec: rec,
		gt:  rec.FormatIndex("GT"),
		dp:  rec.FormatIndex("DP"),
		gq:  rec.FormatIndex("GQ"),
		ad:  rec.FormatIndex("AD"),
	}
}

// Read returns sample i's fields, as ReadSample would. A key the record's
// FORMAT does not declare is Missing (DP/GQ) or "" (GT/AD); one it declares but
// the sample dropped from the end of its column is ".", which DP/GQ read as
// Missing.
func (r SampleReader) Read(i int) SampleFields {
	s := SampleFields{DP: Missing, GQ: Missing}
	if r.gt >= 0 {
		s.GT = r.rec.SampleField(i, r.gt)
	}
	if r.dp >= 0 {
		s.DP = atoi32(r.rec.SampleField(i, r.dp))
	}
	if r.gq >= 0 {
		s.GQ = atoi32(r.rec.SampleField(i, r.gq))
	}
	if r.ad >= 0 {
		s.AD = r.rec.SampleField(i, r.ad)
	}
	return s
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
	// Most genotypes in a wide cohort carry nothing; deciding that without
	// building the recoded genotype keeps them from allocating.
	if !carriesAllele(sf.GT, altIdx) {
		return Call{}, false
	}
	gt, _ := SplitGT(sf.GT, altIdx)
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
