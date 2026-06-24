package annotate

import "github.com/compgenlab/hts/vcf"

// base provides alt-coordinate resolution for annotators that query external
// files by genomic locus. The query chromosome/position can be taken from the
// record directly or, when configured, from INFO fields (--alt-chrom/--alt-pos/
// --end-pos), which is useful for structural variants. It ports ngsutilsj
// AbstractBasicAnnotator.getChrom/getPos/getEndPos and implements [CoordAware].
type base struct {
	closeNoop
	altChrom string
	altPos   string
	endPos   string
}

// SetAltChrom configures the INFO field used as the query chromosome.
func (b *base) SetAltChrom(key string) { b.altChrom = key }

// SetAltPos configures the INFO field used as the query position.
func (b *base) SetAltPos(key string) { b.altPos = key }

// SetEndPos configures the INFO field used as the query end position.
func (b *base) SetEndPos(key string) { b.endPos = key }

// Chrom resolves the query chromosome. ok is false when an alt-chrom INFO field
// is configured but absent on the record (the caller should skip annotation).
func (b *base) Chrom(rec *vcf.VcfRecord) (string, bool) {
	if b.altChrom == "" {
		return rec.Chrom, true
	}
	v, ok := rec.InfoValue(b.altChrom)
	if !ok {
		return "", false
	}
	return v.String(), true
}

// Pos resolves the 1-based query position. With no alt-pos field, it is the
// record position for an SNV and position+1 for a deletion (the variant is the
// next base). With an alt-pos field, it is that INFO value.
func (b *base) Pos(rec *vcf.VcfRecord) (int, bool) {
	if b.altPos == "" {
		if len(rec.Ref) == 1 {
			return rec.Pos, true
		}
		return rec.Pos + 1, true
	}
	v, ok := rec.InfoValue(b.altPos)
	if !ok {
		return 0, false
	}
	n, err := v.Int()
	if err != nil {
		return 0, false
	}
	return n, true
}

// EndPos resolves the 1-based query end position. With no end-pos field, it is
// the query position for an SNV and pos-1+len(ref) for a deletion. With an
// end-pos field, it is that INFO value.
func (b *base) EndPos(rec *vcf.VcfRecord) (int, bool) {
	if b.endPos == "" {
		if len(rec.Ref) == 1 {
			return b.Pos(rec)
		}
		p, ok := b.Pos(rec)
		if !ok {
			return 0, false
		}
		return p - 1 + len(rec.Ref), true
	}
	v, ok := rec.InfoValue(b.endPos)
	if !ok {
		return 0, false
	}
	n, err := v.Int()
	if err != nil {
		return 0, false
	}
	return n, true
}
