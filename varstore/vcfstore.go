package varstore

import (
	"fmt"
	"io"
	"math"
	"strings"

	"github.com/compgenlab/cghts/htsio"
	"github.com/compgenlab/cghts/vcf"
)

// VcfStore is a Store backed by a VCF file.
//
// A joint-called VCF needs no sidecars to classify: it carries an explicit
// genotype for every sample at every record, so a 0/0 is a positive reference
// observation and a ./. is a positive statement of absence. That is precisely
// what the Parquet store has to reconstruct from its sites and regions files.
type VcfStore struct {
	path    string
	samples []string
}

// OpenVcf opens a VCF store. Region-scoped queries additionally require a
// tabix index next to the file.
func OpenVcf(path string) (*VcfStore, error) {
	r, err := vcf.NewVcfFile(path)
	if err != nil {
		return nil, err
	}
	defer r.Close()
	h, err := r.Header()
	if err != nil {
		return nil, err
	}
	return &VcfStore{path: path, samples: h.Samples()}, nil
}

// Samples returns the header sample list.
func (s *VcfStore) Samples() ([]string, error) { return s.samples, nil }

// scan walks records, optionally restricted to a span via the tabix index.
//
// strictRef distinguishes an explicit user assertion from a question. A
// --region names a contig the caller believes is there, so an unresolvable name
// is an error. A variant lookup merely asks whether something is present, and a
// contig the file does not have is simply an absence -- answered with no
// records, and reported by the caller as not-assayed.
func (s *VcfStore) scan(span *Span, strictRef bool, fn func(*vcf.VcfRecord) (bool, error)) error {
	if span != nil {
		ir, err := vcf.NewIndexedVcfReader(s.path)
		if err != nil {
			return fmt.Errorf("--region requires a tabix-indexed VCF: %w", err)
		}
		defer ir.Close()
		end := int(span.End)
		if end < 0 {
			end = math.MaxInt32
		}
		ref, err := resolveRef(ir, span.Chrom)
		if err != nil {
			if strictRef {
				return err
			}
			return nil // contig absent from this file: no records, not a failure
		}
		seq, err := ir.Query(ref, int(span.Start), end)
		if err != nil {
			return err
		}
		for rec, qerr := range seq {
			if qerr != nil {
				return qerr
			}
			cont, err := fn(rec)
			if err != nil {
				return err
			}
			if !cont {
				return nil
			}
		}
		return nil
	}

	r, err := vcf.NewVcfFile(s.path)
	if err != nil {
		return err
	}
	defer r.Close()
	if _, err := r.Header(); err != nil {
		return err
	}
	for {
		rec, err := r.NextRecord()
		if err == io.EOF {
			return nil
		}
		if err != nil {
			return err
		}
		cont, err := fn(rec)
		if err != nil {
			return err
		}
		if !cont {
			return nil
		}
	}
}

// spanFor returns a one-position span for a locus, letting an indexed VCF seek
// instead of scanning. Falls back to a full scan when unindexed. The span is
// 0-based half-open, so a 1-based locus position becomes [Pos-1, Pos).
func (s *VcfStore) spanFor(l Locus) *Span {
	if !fileExists(s.path+".tbi") && !fileExists(s.path+".csi") {
		return nil
	}
	return &Span{Chrom: l.Chrom, Start: l.Pos - 1, End: l.Pos}
}

// Carriers returns the gated ALT calls at a locus.
func (s *VcfStore) Carriers(l Locus, g Gate) ([]Call, error) {
	var out []Call
	err := s.scan(s.spanFor(l), false, func(rec *vcf.VcfRecord) (bool, error) {
		if !SameChrom(rec.Chrom, l.Chrom) || int32(rec.Pos) != l.Pos || rec.Ref != l.Ref {
			return true, nil
		}
		altIdx := altIndex(rec, l.Alt)
		if altIdx < 0 {
			return true, nil
		}
		for i, name := range s.samples {
			if i >= rec.NumSamples() {
				break
			}
			sf, err := ReadSample(rec, i)
			if err != nil {
				return false, err
			}
			c, ok := CallFor(rec, name, sf, altIdx, l.Alt)
			if ok && g.Admits(c) {
				out = append(out, c)
			}
		}
		return true, nil
	})
	return out, err
}

// Variants returns the gated ALT calls for one sample.
func (s *VcfStore) Variants(sample string, span *Span, g Gate) ([]Call, error) {
	idx := -1
	for i, n := range s.samples {
		if n == sample {
			idx = i
			break
		}
	}
	if idx < 0 {
		return nil, fmt.Errorf("sample not found: %s", sample)
	}
	var out []Call
	err := s.scan(span, true, func(rec *vcf.VcfRecord) (bool, error) {
		if idx >= rec.NumSamples() {
			return true, nil
		}
		sf, err := ReadSample(rec, idx)
		if err != nil {
			return false, err
		}
		if !IsAltCarrier(sf.GT) {
			return true, nil
		}
		for j, alt := range rec.Alt() {
			c, ok := CallFor(rec, sample, sf, j+1, alt)
			if ok && g.Admits(c) {
				out = append(out, c)
			}
		}
		return true, nil
	})
	return out, err
}

// Classify resolves every sample at a locus directly from the genotypes.
func (s *VcfStore) Classify(l Locus, g Gate) ([]SampleState, error) {
	states := make(map[string]SampleState, len(s.samples))
	found := false

	err := s.scan(s.spanFor(l), false, func(rec *vcf.VcfRecord) (bool, error) {
		if !SameChrom(rec.Chrom, l.Chrom) || int32(rec.Pos) != l.Pos || rec.Ref != l.Ref {
			return true, nil
		}
		altIdx := altIndex(rec, l.Alt)
		if altIdx < 0 {
			return true, nil
		}
		found = true
		for i, name := range s.samples {
			if i >= rec.NumSamples() {
				break
			}
			sf, err := ReadSample(rec, i)
			if err != nil {
				return false, err
			}
			st := SampleState{SampleID: name}
			if c, ok := CallFor(rec, name, sf, altIdx, l.Alt); ok {
				cc := c
				st.Call = &cc
				if g.Admits(c) {
					st.State = StateCarrier
				} else {
					st.State = StateUncertain
				}
			} else if IsHomRef(sf.GT) || IsAltCarrier(sf.GT) {
				// An explicit genotype that is not this allele: the sample was
				// assayed here and does not carry it -- but only if the call
				// clears the gate. A 0/0 at DP 3 under --min-dp 10 is not a
				// reference observation we are willing to make, and treating it
				// as one is what would let a poorly covered sample quietly
				// enlarge the non-carrier denominator.
				if g.Admits(Call{DP: sf.DP, GQ: sf.GQ}) {
					st.State = StateNonCarrier
				} else {
					st.State = StateNotAssayed
				}
			} else {
				st.State = StateNotAssayed
			}
			states[name] = st
		}
		return false, nil
	})
	if err != nil {
		return nil, err
	}

	out := make([]SampleState, 0, len(s.samples))
	for _, name := range s.samples {
		if st, ok := states[name]; ok {
			out = append(out, st)
			continue
		}
		// The record was absent entirely: nothing was interrogated here.
		_ = found
		out = append(out, SampleState{SampleID: name, State: StateNotAssayed})
	}
	return out, nil
}

// SiteKnown reports whether the VCF actually contains a record for this exact
// REF/ALT at this position. A plain VCF asserts nothing about anything else --
// an unreported base was not observed to be reference, it was simply not
// reported -- so this is the boundary of what the store can answer.
func (s *VcfStore) SiteKnown(l Locus) (bool, error) {
	found := false
	err := s.scan(s.spanFor(l), false, func(rec *vcf.VcfRecord) (bool, error) {
		if SameChrom(rec.Chrom, l.Chrom) && int32(rec.Pos) == l.Pos &&
			rec.Ref == l.Ref && altIndex(rec, l.Alt) > 0 {
			found = true
			return false, nil
		}
		return true, nil
	})
	return found, err
}

// Close is a no-op; VcfStore opens the file per query.
func (s *VcfStore) Close() error { return nil }

// altIndex returns the 1-based index of alt in the record's ALT list, or -1.
func altIndex(rec *vcf.VcfRecord, alt string) int {
	for i, a := range rec.Alt() {
		if a == alt {
			return i + 1
		}
	}
	return -1
}

// resolveRef maps a chromosome name onto whatever the tabix index actually
// calls that contig.
//
// A query must not have to know the file's naming convention. ParseLocus keeps
// the user's spelling and comparisons are canonical, but tabix looks names up
// verbatim, so "chr22" against an Ensembl-named index -- or "22" against a
// UCSC-named one -- fails at the index rather than returning nothing. cghts's
// ContigConverter resolves by canonical identity across UCSC, Ensembl and NCBI
// RefSeq spellings, built from the index's own reference list.
func resolveRef(ir *vcf.IndexedVcfReader, name string) (string, error) {
	if ir.HasRef(name) {
		return name, nil
	}
	refs := ir.RefNames()
	if got, ok := vcf.NewContigConverter(refs).Resolve(name); ok {
		return got, nil
	}
	shown := refs
	if len(shown) > 8 {
		shown = shown[:8]
	}
	more := ""
	if len(refs) > len(shown) {
		more = fmt.Sprintf(" (and %d more)", len(refs)-len(shown))
	}
	return "", fmt.Errorf("unknown reference %q; the index has %s%s",
		name, strings.Join(shown, ", "), more)
}

// ParseSpan parses a 1-based inclusive "chrom:start-end" or bare "chrom".
func ParseSpan(region string) (*Span, error) {
	if region == "" {
		return nil, nil
	}
	ref, start, end, err := htsio.ParseRegion(region)
	if err != nil {
		return nil, err
	}
	if end < 0 {
		end = math.MaxInt32
	}
	return &Span{Chrom: ref, Start: int32(start), End: int32(end)}, nil
}
