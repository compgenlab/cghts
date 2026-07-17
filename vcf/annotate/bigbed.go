package annotate

import (
	"fmt"
	"strings"

	"github.com/compgenlab/cghts/htsio/bbi"
	"github.com/compgenlab/cghts/vcf"
)

// BigBedOptions configures a [BigBedAnnotator]. Columns are 1-based over the BED
// row (1=chrom, 2=start, 3=end, 4=name, …); Col=0 means a presence flag. The
// input file is a local, self-indexed .bb (no tabix).
type BigBedOptions struct {
	Name     string // INFO/FORMAT key to add
	Filename string // bigBed (.bb) file
	Sample   string // "" = INFO; otherwise a FORMAT field for this sample
	Col      int    // 1-based value column; 0 = presence flag
	IsNumber bool   // declare the value Float
	Collapse bool   // join unique values with ","
	First    bool   // keep only the first value
	Extend   int    // widen the query by N bases on each side
	NoHeader bool   // do not add a ##INFO/##FORMAT def

	// AutoConvert matches contig names across UCSC/Ensembl/NCBI naming.
	AutoConvert bool
}

// BigBedAnnotator adds an INFO or FORMAT annotation from a bigBed file — a column
// from the BED interval(s) overlapping the variant. It mirrors [TabixAnnotator]
// but reads a self-indexed BBI file.
type BigBedAnnotator struct {
	base
	opts      BigBedOptions
	reader    *bbi.Reader
	col       int // 0-based value column; -1 = flag
	sampleIdx int
	conv      *vcf.ContigConverter
}

// NewBigBedAnnotator opens the bigBed file and returns the annotator.
func NewBigBedAnnotator(opts BigBedOptions) (*BigBedAnnotator, error) {
	r, err := bbi.Open(opts.Filename)
	if err != nil {
		return nil, fmt.Errorf("annotate: open %s: %w", opts.Filename, err)
	}
	if r.Kind() != bbi.BigBed {
		r.Close()
		return nil, fmt.Errorf("annotate: %s is not a bigBed file", opts.Filename)
	}
	a := &BigBedAnnotator{opts: opts, reader: r, col: opts.Col - 1, sampleIdx: -1}
	if opts.AutoConvert {
		a.EnableContigMatching()
	}
	return a, nil
}

// EnableContigMatching turns on cross-scheme contig-name matching. Implements
// [ContigMatcher].
func (a *BigBedAnnotator) EnableContigMatching() {
	a.conv = vcf.NewContigConverter(a.reader.RefNames())
}

// SetupHeader resolves the sample (for FORMAT) and adds the ##INFO/##FORMAT def.
func (a *BigBedAnnotator) SetupHeader(h *vcf.VcfHeader) error {
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
		h.AddInfo(infoDefSrc(a.opts.Name, "0", "Flag", "Present in bigBed file", a.opts.Filename))
		return nil
	}
	typ := "String"
	if a.opts.IsNumber {
		typ = "Float"
	}
	desc := fmt.Sprintf("Column %d from bigBed", a.col+1)
	if a.opts.Sample != "" {
		h.AddFormat(formatDefSrc(a.opts.Name, ".", typ, desc, a.opts.Filename))
	} else {
		h.AddInfo(infoDefSrc(a.opts.Name, ".", typ, desc, a.opts.Filename))
	}
	return nil
}

// Annotate queries the bigBed for overlapping intervals and adds the annotation.
func (a *BigBedAnnotator) Annotate(rec *vcf.VcfRecord) error {
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
	seq, err := a.reader.Query(chrom, pos-1-a.opts.Extend, endpos+a.opts.Extend)
	if err != nil {
		return err
	}
	var vals []string
	found := false
	for r, err := range seq {
		if err != nil {
			return err
		}
		found = true
		if a.col >= 0 {
			fields := strings.Split(r.Line, "\t")
			if a.col < len(fields) && fields[a.col] != "" {
				vals = append(vals, fields[a.col])
			}
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
	if a.opts.Sample != "" {
		return rec.AddFormat(a.sampleIdx, a.opts.Name, a.aggregate(vals))
	}
	rec.AddInfo(a.opts.Name, a.aggregate(vals))
	return nil
}

// Close releases the bigBed reader.
func (a *BigBedAnnotator) Close() error { return a.reader.Close() }

func (a *BigBedAnnotator) aggregate(vals []string) string {
	switch {
	case a.opts.Collapse:
		return strings.Join(uniqueStrings(vals), ",")
	case a.opts.First:
		return vals[0]
	default:
		return strings.Join(vals, ",")
	}
}
