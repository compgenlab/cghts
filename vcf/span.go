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

// IsRefBlock reports whether this record is a gVCF reference block rather than a
// variant: a <*> or <NON_REF> alternate, or a bare "." where no alternate was
// reported at all.
//
// Worth having separately from RefSpan because the two questions diverge. A caller
// listing variants wants to skip these; a caller asking what was interrogated wants
// exactly these. Both need to recognise them the same way.
func (r *VcfRecord) IsRefBlock() bool {
	alt := r.AltOrig()
	if alt == "." || alt == "" {
		return true
	}
	for _, a := range r.Alt() {
		if vcfspan.IsRefBlockAlt(a) {
			return true
		}
	}
	return false
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
