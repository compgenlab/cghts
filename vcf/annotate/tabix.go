package annotate

import (
	"fmt"
	"strconv"
	"strings"

	"github.com/compgenlab/hts/htsio/tabix"
	"github.com/compgenlab/hts/vcf"
)

// TabixOptions configures a [TabixAnnotator]. Columns are 1-based; Col=0 means
// the annotation is a presence flag (no value). The input file must be
// BGZF-compressed and tabix-indexed (.tbi/.csi). A BED annotation is just
// Col=4 (the BED name column).
type TabixOptions struct {
	Name     string // INFO/FORMAT key to add
	Filename string // tabix-indexed file
	Sample   string // "" = INFO; otherwise a FORMAT field for this sample
	Col      int    // 1-based value column; 0 = presence flag
	AltCol   int    // 1-based ALT-match column; 0 = none
	RefCol   int    // 1-based REF-match column; 0 = none
	IsNumber bool   // declare the value Float (required by Max)
	Collapse bool   // join unique values with ","
	First    bool   // keep only the first value
	Max      bool   // keep the numeric maximum (".0"-trimmed)
	Extend   int    // widen the query by N bases on each side
	NoHeader bool   // do not add a ##INFO/##FORMAT def
}

// TabixAnnotator adds INFO or FORMAT annotations from a tabix-indexed file. It
// generalizes BED annotation (use the name column) to any column, with optional
// alt/ref-allele exact matching and value aggregation. It ports ngsutilsj
// TabixAnnotation (and BEDAnnotation, as a Col=4 preset).
type TabixAnnotator struct {
	base
	opts      TabixOptions
	reader    *tabix.Reader
	col       int // 0-based value column; -1 = flag
	altCol    int // 0-based; -1 = none
	refCol    int // 0-based; -1 = none
	sampleIdx int
}

// NewTabixAnnotator opens the tabix-indexed file and returns the annotator.
func NewTabixAnnotator(opts TabixOptions) (*TabixAnnotator, error) {
	r, err := tabix.NewReader(opts.Filename)
	if err != nil {
		return nil, fmt.Errorf("annotate: open %s: %w", opts.Filename, err)
	}
	return &TabixAnnotator{
		opts:      opts,
		reader:    r,
		col:       opts.Col - 1,
		altCol:    opts.AltCol - 1,
		refCol:    opts.RefCol - 1,
		sampleIdx: -1,
	}, nil
}

// SetupHeader resolves the sample (for FORMAT) and adds the ##INFO/##FORMAT def.
func (a *TabixAnnotator) SetupHeader(h *vcf.VcfHeader) error {
	if a.opts.Sample != "" {
		a.sampleIdx = h.SampleIndex(a.opts.Sample)
		if a.sampleIdx < 0 {
			return fmt.Errorf("annotate: missing sample: %s", a.opts.Sample)
		}
	}
	if a.opts.NoHeader {
		return nil
	}
	if a.col < 0 {
		h.AddInfo(infoDefSrc(a.opts.Name, "0", "Flag", "Present in Tabix file", a.opts.Filename))
		return nil
	}
	typ := "String"
	if a.opts.IsNumber {
		typ = "Float"
	}
	desc := fmt.Sprintf("Column %d from file", a.col+1)
	if a.opts.Sample != "" {
		h.AddFormat(formatDefSrc(a.opts.Name, ".", typ, desc, a.opts.Filename))
	} else {
		h.AddInfo(infoDefSrc(a.opts.Name, ".", typ, desc, a.opts.Filename))
	}
	return nil
}

// Annotate queries the tabix file for overlapping rows and adds the annotation.
func (a *TabixAnnotator) Annotate(rec *vcf.VcfRecord) error {
	chrom, ok := a.Chrom(rec)
	if !ok {
		return nil
	}
	var pos, endpos int
	if a.refCol >= 0 {
		// Matching against a ref/alt column: use the variant position directly
		// (equivalent to a VCF comparison), no SNV/deletion adjustment.
		pos, endpos = rec.Pos, rec.Pos
	} else {
		var ok1, ok2 bool
		if pos, ok1 = a.Pos(rec); !ok1 {
			return nil
		}
		if endpos, ok2 = a.EndPos(rec); !ok2 {
			return nil
		}
	}

	seq, err := a.reader.Query(chrom, pos-1-a.opts.Extend, endpos+a.opts.Extend)
	if err != nil {
		return err
	}

	var vals []string
	found := false
	for tr, err := range seq {
		if err != nil {
			return err
		}
		fields := strings.Split(tr.Line, "\t")
		if a.altCol >= 0 {
			altOk := false
			if a.altCol < len(fields) {
				for _, alt := range rec.Alt() {
					if alt == fields[a.altCol] {
						altOk = true
						break
					}
				}
			}
			if !altOk {
				continue
			}
		}
		if a.refCol >= 0 {
			if !(a.refCol < len(fields) && rec.Ref == fields[a.refCol]) {
				continue
			}
		}
		found = true
		if a.col >= 0 && a.col < len(fields) && fields[a.col] != "" {
			vals = append(vals, fields[a.col])
		}
	}

	if !found {
		return nil
	}
	if a.col < 0 {
		rec.AddInfoFlag(a.opts.Name)
		return nil
	}
	if len(vals) == 0 {
		return nil
	}
	out, err := a.aggregate(vals)
	if err != nil {
		return err
	}
	if a.opts.Sample != "" {
		return rec.AddFormat(a.sampleIdx, a.opts.Name, out)
	}
	rec.AddInfo(a.opts.Name, out)
	return nil
}

func (a *TabixAnnotator) aggregate(vals []string) (string, error) {
	switch {
	case a.opts.Collapse:
		return strings.Join(uniqueStrings(vals), ","), nil
	case a.opts.First:
		return vals[0], nil
	case a.opts.Max:
		m, err := strconv.ParseFloat(vals[0], 64)
		if err != nil {
			return "", fmt.Errorf("annotate: non-numeric value %q for max: %w", vals[0], err)
		}
		for _, v := range vals[1:] {
			d, err := strconv.ParseFloat(v, 64)
			if err != nil {
				return "", fmt.Errorf("annotate: non-numeric value %q for max: %w", v, err)
			}
			if d > m {
				m = d
			}
		}
		s := strconv.FormatFloat(m, 'f', -1, 64)
		return strings.TrimSuffix(s, ".0"), nil
	default:
		return strings.Join(vals, ","), nil
	}
}

// Close releases the tabix reader.
func (a *TabixAnnotator) Close() error { return a.reader.Close() }

// uniqueStrings returns the values with duplicates removed, preserving order.
func uniqueStrings(vals []string) []string {
	seen := make(map[string]bool, len(vals))
	out := make([]string, 0, len(vals))
	for _, v := range vals {
		if !seen[v] {
			seen[v] = true
			out = append(out, v)
		}
	}
	return out
}
