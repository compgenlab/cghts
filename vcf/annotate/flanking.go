package annotate

import (
	"fmt"
	"strings"

	"github.com/compgenlab/hts/seqio"
	"github.com/compgenlab/hts/support/sequtils"
	"github.com/compgenlab/hts/vcf"
)

// FlankingOptions configures a [FlankingBases] annotator.
type FlankingOptions struct {
	Filename string // reference FASTA (indexed, optionally bgzip-compressed)
	Size     int    // flanking bases on each side (default 1)
}

// FlankingBases adds the reference context around an SNV (CG_FLANKING) and the
// normalized substitution (CG_FLANKING_SUB, e.g. A[C>A]A). The substitution is
// reported on the pyrimidine strand: when the REF is A or G the context and
// alleles are reverse-complemented. Indels and variants too close to a sequence
// end are skipped. It ports ngsutilsj FlankingBases (--flanking).
type FlankingBases struct {
	ref  seqio.ReferenceReader
	size int
}

// NewFlankingBases opens the reference FASTA and returns the annotator.
func NewFlankingBases(opts FlankingOptions) (*FlankingBases, error) {
	size := opts.Size
	if size <= 0 {
		size = 1
	}
	ref, err := seqio.OpenReference(opts.Filename)
	if err != nil {
		return nil, fmt.Errorf("annotate: open %s: %w", opts.Filename, err)
	}
	return &FlankingBases{ref: ref, size: size}, nil
}

// SetupHeader declares CG_FLANKING and CG_FLANKING_SUB.
func (a *FlankingBases) SetupHeader(h *vcf.VcfHeader) error {
	h.AddInfo(infoDef("CG_FLANKING", "1", "String", fmt.Sprintf("+/- %d bp flanking the variant (no indels)", a.size)))
	h.AddInfo(infoDef("CG_FLANKING_SUB", "A", "String", "Substitution caused by variant "))
	return nil
}

// Annotate fetches the flanking reference context for an SNV.
func (a *FlankingBases) Annotate(rec *vcf.VcfRecord) error {
	if len(rec.Ref) > 1 { // skip indels
		return nil
	}
	length, ok := a.ref.SequenceLength(rec.Chrom)
	if !ok {
		return nil
	}
	start := rec.Pos - 1 - a.size
	end := rec.Pos + a.size
	if start < 0 || end > length {
		// no flanking available past the first/last base
		return nil
	}
	seq, err := a.ref.GetSequenceRange(rec.Chrom, start, end)
	if err != nil {
		return err
	}
	refSeq := string(seq)
	rec.AddInfo("CG_FLANKING", refSeq)

	work := refSeq
	revcomp := false
	switch strings.ToUpper(rec.Ref) {
	case "A", "G":
		work = sequtils.ReverseComplement(refSeq)
		revcomp = true
	}
	pre := work[:a.size]
	varBase := work[a.size : a.size+1]
	post := work[a.size+1:]

	var outs []string
	for _, alt := range rec.Alt() {
		if len(alt) > 1 {
			continue
		}
		if revcomp {
			alt = sequtils.ReverseComplement(alt)
		}
		outs = append(outs, pre+"["+varBase+">"+alt+"]"+post)
	}
	if len(outs) > 0 {
		rec.AddInfo("CG_FLANKING_SUB", strings.Join(outs, ","))
	}
	return nil
}

// Close releases the reference reader.
func (a *FlankingBases) Close() error { return a.ref.Close() }
