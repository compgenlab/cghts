package annotate

import (
	"fmt"
	"strconv"
	"strings"

	"github.com/compgenlab/hts/htsio/bbi"
	"github.com/compgenlab/hts/vcf"
)

// BigWigOptions configures a [BigWigAnnotator]. A bigWig holds one value per base
// (no ref/alt), so the annotation is a single numeric INFO/FORMAT value at the
// variant position. The input file is a local, self-indexed .bw (no tabix).
type BigWigOptions struct {
	Name     string // INFO/FORMAT key to add
	Filename string // bigWig (.bw) file
	Sample   string // "" = INFO; otherwise a FORMAT field for this sample
	NoHeader bool   // do not add a ##INFO/##FORMAT def

	// Aggregate selects how multiple overlapping values are combined (a point
	// query normally yields one): "first" (default), "max", or "mean".
	Aggregate string

	// AutoConvert matches contig names across UCSC/Ensembl/NCBI naming instead of
	// requiring an exact-string match.
	AutoConvert bool
}

// BigWigAnnotator adds a numeric INFO or FORMAT annotation from a bigWig file —
// the base-resolution value at the variant position. It mirrors [TabixAnnotator]
// but reads a self-indexed BBI file.
type BigWigAnnotator struct {
	base
	opts      BigWigOptions
	reader    *bbi.Reader
	sampleIdx int
	conv      *vcf.ContigConverter
}

// NewBigWigAnnotator opens the bigWig file and returns the annotator.
func NewBigWigAnnotator(opts BigWigOptions) (*BigWigAnnotator, error) {
	r, err := bbi.Open(opts.Filename)
	if err != nil {
		return nil, fmt.Errorf("annotate: open %s: %w", opts.Filename, err)
	}
	if r.Kind() != bbi.BigWig {
		r.Close()
		return nil, fmt.Errorf("annotate: %s is not a bigWig file", opts.Filename)
	}
	a := &BigWigAnnotator{opts: opts, reader: r, sampleIdx: -1}
	if opts.AutoConvert {
		a.EnableContigMatching()
	}
	return a, nil
}

// EnableContigMatching turns on cross-scheme contig-name matching, built from the
// bigWig's contig names. Implements [ContigMatcher].
func (a *BigWigAnnotator) EnableContigMatching() {
	a.conv = vcf.NewContigConverter(a.reader.RefNames())
}

// SetupHeader resolves the sample (for FORMAT) and adds the ##INFO/##FORMAT def.
func (a *BigWigAnnotator) SetupHeader(h *vcf.VcfHeader) error {
	if a.opts.Sample != "" {
		a.sampleIdx = h.SampleIndex(a.opts.Sample)
		if a.sampleIdx < 0 {
			return fmt.Errorf("annotate: missing sample: %s", a.opts.Sample)
		}
	}
	if a.opts.NoHeader {
		return nil
	}
	desc := "bigWig value at position"
	if a.opts.Sample != "" {
		h.AddFormat(formatDefSrc(a.opts.Name, "1", "Float", desc, a.opts.Filename))
	} else {
		h.AddInfo(infoDefSrc(a.opts.Name, "1", "Float", desc, a.opts.Filename))
	}
	return nil
}

// Annotate queries the bigWig for the variant position and adds the value.
func (a *BigWigAnnotator) Annotate(rec *vcf.VcfRecord) error {
	chrom, ok := a.Chrom(rec)
	if !ok {
		return nil
	}
	pos, ok := a.Pos(rec)
	if !ok {
		return nil
	}
	endpos, ok := a.EndPos(rec)
	if !ok {
		return nil
	}
	if a.conv != nil {
		if chrom, ok = a.conv.Resolve(chrom); !ok {
			return nil
		}
	} else if !a.reader.HasRef(chrom) {
		return nil
	}
	seq, err := a.reader.Query(chrom, pos-1, endpos) // 0-based half-open
	if err != nil {
		return err
	}
	var vals []float64
	for r, err := range seq {
		if err != nil {
			return err
		}
		vals = append(vals, r.Value)
	}
	if len(vals) == 0 {
		return nil
	}
	out := formatFloat(aggregateFloat(vals, a.opts.Aggregate))
	if a.opts.Sample != "" {
		return rec.AddFormat(a.sampleIdx, a.opts.Name, out)
	}
	rec.AddInfo(a.opts.Name, out)
	return nil
}

// Close releases the bigWig reader.
func (a *BigWigAnnotator) Close() error { return a.reader.Close() }

// aggregateFloat combines overlapping bigWig values per the Aggregate mode.
func aggregateFloat(vals []float64, mode string) float64 {
	switch mode {
	case "max":
		m := vals[0]
		for _, v := range vals[1:] {
			if v > m {
				m = v
			}
		}
		return m
	case "mean":
		var s float64
		for _, v := range vals {
			s += v
		}
		return s / float64(len(vals))
	default: // "first"
		return vals[0]
	}
}

// formatFloat renders a bigWig value. bigWig stores float32, so it is formatted
// at float32 precision (bitSize 32) — this yields the clean decimal the file
// author wrote (e.g. "0.42") rather than its float64-widened artifact.
func formatFloat(f float64) string {
	return strings.TrimSuffix(strconv.FormatFloat(f, 'g', -1, 32), ".0")
}
