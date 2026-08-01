package vcf

import "github.com/compgenlab/cghts/internal/vcfspan"

// RefSpan returns the reference bases a record covers, in 0-based half-open
// coordinates.
//
// This is not the same question as AltPositions answers. That resolves an SV's
// *partner breakpoint* -- where the other end of a translocation or the mate of a
// breakend lies -- which may be on another chromosome and may precede POS. This
// reports how much of the reference this one record describes, which is what an
// interval query has to compare against.
//
// The span is the widest of len(REF), INFO/END, INFO/SVLEN for symbolic alternates
// measured on the reference, and FORMAT/LEN for reference blocks; see
// internal/vcfspan for the precedence and why each is scoped the way it is. The
// same code backs the tabix index writer, so an index and a query agree on how far
// a record reaches.
func (r *VcfRecord) RefSpan() (start, end int) {
	start = r.ZeroBasedStart()
	return start, vcfspan.End(recordFields{r}, start)
}

// RefSpanEnd is RefSpan's end alone, for callers that already know the start.
func (r *VcfRecord) RefSpanEnd() int {
	_, end := r.RefSpan()
	return end
}

// IsRefBlock reports whether this record describes only reference -- a gVCF block
// -- rather than a variant.
//
// The test is that *no* alternate is a real allele: every one is <*> or <NON_REF>,
// or ALT is a bare ".". The distinction matters because GATK writes the block
// allele alongside a genuine alternate at variant sites, as in ALT "G,<NON_REF>".
// Such a record does describe a variant, so treating the presence of <NON_REF> as
// sufficient would make a caller listing variants drop real calls.
//
// Worth having separately from RefSpan because the two questions diverge: a caller
// listing variants skips these, a caller asking what was interrogated wants exactly
// these, and both need to recognise them the same way.
func (r *VcfRecord) IsRefBlock() bool {
	if alt := r.AltOrig(); alt == "." || alt == "" {
		return true
	}
	alts := r.Alt()
	if len(alts) == 0 {
		return true
	}
	for _, a := range alts {
		if !IsRefBlockAlt(a) {
			return false
		}
	}
	return true
}

// IsRefBlockAlt reports whether one alternate is a gVCF reference-block allele,
// <*> or <NON_REF>, rather than a real alternate. Exported because a caller walking
// a record's alternates needs to skip these individually -- a "G,<NON_REF>" record
// has one of each.
func IsRefBlockAlt(alt string) bool { return vcfspan.IsRefBlockAlt(alt) }

// BlockDepth returns the depth floor a record vouches for at sample i, and whether
// anything is known.
//
// For a gVCF reference block this is `FORMAT/MIN_DP`: the minimum depth anywhere in
// the block, and therefore the only depth claim that holds across the whole span.
// `FORMAT/DP` is the fallback, but it is a **weaker** statement -- the depth at POS
// alone, which says nothing about the rest of the block -- so a caller gating a
// span-wide assertion on it is overclaiming. Reported separately from DP for that
// reason rather than folded into one accessor.
func (r *VcfRecord) BlockDepth(i int) (int32, bool) {
	if v, ok := r.sampleInt(i, "MIN_DP"); ok {
		return v, true
	}
	return r.sampleInt(i, "DP")
}

// BlockRGQ returns `FORMAT/RGQ` for sample i: the confidence that the sample really
// is reference here, which is what GQ is for a genotype call. GATK writes it on
// reference blocks in place of GQ, so a caller gating reference calls on quality has
// to look here rather than at GQ.
func (r *VcfRecord) BlockRGQ(i int) (int32, bool) {
	if v, ok := r.sampleInt(i, "RGQ"); ok {
		return v, true
	}
	return r.sampleInt(i, "GQ")
}

// sampleInt reads one integer FORMAT value, treating absent, empty and "." alike as
// unknown. Absence is not zero: a zero depth is a real claim and unknown depth is
// not, and conflating them makes a gate silently admit everything.
func (r *VcfRecord) sampleInt(i int, key string) (int32, bool) {
	s, err := r.Sample(i)
	if err != nil {
		return 0, false
	}
	v, ok := s.Get(key)
	if !ok || v.IsMissing() || v.IsEmpty() {
		return 0, false
	}
	n, err := v.Int()
	if err != nil {
		return 0, false
	}
	return int32(n), true
}

// recordFields adapts a parsed record to the shared span rules. It reads through
// the record's own lazy accessors rather than re-splitting the raw line, so a
// record built in memory works as well as one read from a file.
type recordFields struct{ r *VcfRecord }

func (f recordFields) Ref() string { return f.r.Ref }

func (f recordFields) Alts() []string { return f.r.Alt() }

func (f recordFields) Info(key string) (string, bool) {
	v, ok := f.r.Info().Get(key)
	if !ok {
		return "", false
	}
	return v.String(), true
}

func (f recordFields) SampleValues(key string) []string {
	n := f.r.NumSamples()
	if n == 0 {
		return nil
	}
	out := make([]string, 0, n)
	for i := 0; i < n; i++ {
		s, err := f.r.Sample(i)
		if err != nil {
			continue
		}
		if v, ok := s.Get(key); ok {
			if str := v.String(); str != "" && str != "." {
				out = append(out, str)
			}
		}
	}
	return out
}
