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
	ColName  string // value column by header name (overrides Col when set)
	AltName  string // ALT-match column by header name (overrides AltCol)
	RefName  string // REF-match column by header name (overrides RefCol)
	IsNumber bool   // declare the value Float (required by Max)
	Collapse bool   // join unique values with ","
	First    bool   // keep only the first value
	Max      bool   // keep the numeric maximum (".0"-trimmed)
	Extend   int    // widen the query by N bases on each side
	NoHeader bool   // do not add a ##INFO/##FORMAT def

	// AutoConvert matches contig names across UCSC/Ensembl/NCBI naming (human
	// primary contigs 1-22,X,Y,MT) instead of requiring an exact-string match.
	AutoConvert bool
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
	conv      *vcf.ContigConverter // non-nil when contig-name matching is enabled
}

// NewTabixAnnotator opens the tabix-indexed file and returns the annotator. Any
// column specified by name (ColName/AltName/RefName) is resolved against the
// file's header.
func NewTabixAnnotator(opts TabixOptions) (*TabixAnnotator, error) {
	r, err := tabix.NewReader(opts.Filename)
	if err != nil {
		return nil, fmt.Errorf("annotate: open %s: %w", opts.Filename, err)
	}
	for _, res := range []struct {
		name string
		col  *int
	}{
		{opts.ColName, &opts.Col},
		{opts.AltName, &opts.AltCol},
		{opts.RefName, &opts.RefCol},
	} {
		if res.name == "" {
			continue
		}
		n, err := r.ColumnByName(res.name)
		if err != nil {
			r.Close()
			return nil, fmt.Errorf("annotate: %w", err)
		}
		*res.col = n
	}
	a := &TabixAnnotator{
		opts:      opts,
		reader:    r,
		col:       opts.Col - 1,
		altCol:    opts.AltCol - 1,
		refCol:    opts.RefCol - 1,
		sampleIdx: -1,
	}
	if opts.AutoConvert {
		a.EnableContigMatching()
	}
	return a, nil
}

// EnableContigMatching turns on cross-scheme contig-name matching (UCSC/Ensembl/
// NCBI) by building a converter from the source file's contig names. It
// implements [ContigMatcher] so the CLI's --auto-convert flag can enable it.
func (a *TabixAnnotator) EnableContigMatching() {
	a.conv = vcf.NewContigConverter(a.reader.RefNames())
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

	if a.conv != nil {
		if chrom, ok = a.conv.Resolve(chrom); !ok {
			return nil // no contig in this file shares the query's identity
		}
	} else if !a.reader.HasRef(chrom) {
		return nil // contig not in this file → nothing to annotate
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
	return aggregateVals(vals, a.opts.Collapse, a.opts.First, a.opts.Max)
}

// aggregateVals reduces the matched column values per the collapse/first/max mode
// (default: join with ","). Shared by TabixAnnotator and TabixAnnotationGroup.
func aggregateVals(vals []string, collapse, first, max bool) (string, error) {
	switch {
	case collapse:
		return strings.Join(uniqueStrings(vals), ","), nil
	case first:
		return vals[0], nil
	case max:
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

// TabixFieldOptions is one field rule within a [TabixAnnotationGroup]: the per-field
// knobs of [TabixOptions]. The region-determining columns (AltCol/RefCol) and Extend
// are shared by the group (they are source-level in practice), so they are not here.
type TabixFieldOptions struct {
	Name     string // INFO/FORMAT key to add
	Sample   string // "" = INFO; otherwise a FORMAT field for this sample
	Col      int    // 1-based value column; 0 = presence flag
	ColName  string // value column by header name (overrides Col)
	IsNumber bool   // declare the value Float
	Collapse bool   // join unique values with ","
	First    bool   // keep only the first value
	Max      bool   // keep the numeric maximum
	NoHeader bool   // do not add a def
}

// TabixGroupOptions configures a [TabixAnnotationGroup]: one tabix-indexed source and
// N field rules sharing a single reader, one query, and the same alt/ref match columns
// and Extend per input record.
type TabixGroupOptions struct {
	Filename    string
	AltCol      int // 1-based ALT-match column; 0 = none (shared by all fields)
	RefCol      int // 1-based REF-match column; 0 = none (shared by all fields)
	AltName     string
	RefName     string
	Extend      int  // widen the query by N bases on each side (shared)
	AutoConvert bool // cross-scheme contig-name matching
	Fields      []TabixFieldOptions
}

type tabixRule struct {
	name      string
	sample    string
	sampleIdx int
	col       int // 0-based value column; -1 = flag
	isNumber  bool
	collapse  bool
	first     bool
	max       bool
	noHeader  bool
}

// TabixAnnotationGroup adds several INFO/FORMAT annotations from one tabix source
// using a single reader and one query per input record — the multi-field analogue of
// [TabixAnnotator]. The alt/ref match and query region are computed once and shared;
// each rule contributes its own value column and aggregation. All rules must use the
// same match columns / Extend (true for a cganno source), so the region is uniform.
type TabixAnnotationGroup struct {
	base
	filename string
	reader   *tabix.Reader
	altCol   int // 0-based; -1 = none
	refCol   int // 0-based; -1 = none
	extend   int
	conv     *vcf.ContigConverter
	rules    []tabixRule
}

// NewTabixAnnotationGroup opens the source once and returns a grouped annotator. The
// shared match columns and each field's value column may be given by header name.
func NewTabixAnnotationGroup(opts TabixGroupOptions) (*TabixAnnotationGroup, error) {
	r, err := tabix.NewReader(opts.Filename)
	if err != nil {
		return nil, fmt.Errorf("annotate: open %s: %w", opts.Filename, err)
	}
	for _, res := range []struct {
		name string
		col  *int
	}{
		{opts.AltName, &opts.AltCol},
		{opts.RefName, &opts.RefCol},
	} {
		if res.name == "" {
			continue
		}
		n, err := r.ColumnByName(res.name)
		if err != nil {
			r.Close()
			return nil, fmt.Errorf("annotate: %w", err)
		}
		*res.col = n
	}
	g := &TabixAnnotationGroup{
		filename: opts.Filename,
		reader:   r,
		altCol:   opts.AltCol - 1,
		refCol:   opts.RefCol - 1,
		extend:   opts.Extend,
	}
	for _, f := range opts.Fields {
		col := f.Col
		if f.ColName != "" {
			n, err := r.ColumnByName(f.ColName)
			if err != nil {
				r.Close()
				return nil, fmt.Errorf("annotate: %w", err)
			}
			col = n
		}
		g.rules = append(g.rules, tabixRule{
			name: f.Name, sample: f.Sample, sampleIdx: -1,
			col: col - 1, isNumber: f.IsNumber,
			collapse: f.Collapse, first: f.First, max: f.Max, noHeader: f.NoHeader,
		})
	}
	if opts.AutoConvert {
		g.EnableContigMatching()
	}
	return g, nil
}

// EnableContigMatching turns on cross-scheme contig-name matching for the group.
func (g *TabixAnnotationGroup) EnableContigMatching() {
	g.conv = vcf.NewContigConverter(g.reader.RefNames())
}

// SetupHeader resolves each rule's sample (for FORMAT) and adds its def.
func (g *TabixAnnotationGroup) SetupHeader(h *vcf.VcfHeader) error {
	for i := range g.rules {
		r := &g.rules[i]
		if r.sample != "" {
			r.sampleIdx = h.SampleIndex(r.sample)
			if r.sampleIdx < 0 {
				return fmt.Errorf("annotate: missing sample: %s", r.sample)
			}
		}
		if r.noHeader {
			continue
		}
		if r.col < 0 {
			h.AddInfo(infoDefSrc(r.name, "0", "Flag", "Present in Tabix file", g.filename))
			continue
		}
		typ := "String"
		if r.isNumber {
			typ = "Float"
		}
		desc := fmt.Sprintf("Column %d from file", r.col+1)
		if r.sample != "" {
			h.AddFormat(formatDefSrc(r.name, ".", typ, desc, g.filename))
		} else {
			h.AddInfo(infoDefSrc(r.name, ".", typ, desc, g.filename))
		}
	}
	return nil
}

// Annotate queries the source once and applies every rule. The query region and the
// alt/ref match are computed once per input record (they are shared); each matching
// row contributes to every rule's values. Matches [TabixAnnotator.Annotate].
func (g *TabixAnnotationGroup) Annotate(rec *vcf.VcfRecord) error {
	chrom, ok := g.Chrom(rec)
	if !ok {
		return nil
	}
	var pos, endpos int
	if g.refCol >= 0 {
		pos, endpos = rec.Pos, rec.Pos
	} else {
		var ok1, ok2 bool
		if pos, ok1 = g.Pos(rec); !ok1 {
			return nil
		}
		if endpos, ok2 = g.EndPos(rec); !ok2 {
			return nil
		}
	}
	if g.conv != nil {
		if chrom, ok = g.conv.Resolve(chrom); !ok {
			return nil
		}
	} else if !g.reader.HasRef(chrom) {
		return nil
	}
	seq, err := g.reader.Query(chrom, pos-1-g.extend, endpos+g.extend)
	if err != nil {
		return err
	}

	vals := make([][]string, len(g.rules))
	found := false
	for tr, err := range seq {
		if err != nil {
			return err
		}
		fields := strings.Split(tr.Line, "\t")
		if g.altCol >= 0 {
			altOk := false
			if g.altCol < len(fields) {
				for _, alt := range rec.Alt() {
					if alt == fields[g.altCol] {
						altOk = true
						break
					}
				}
			}
			if !altOk {
				continue
			}
		}
		if g.refCol >= 0 {
			if !(g.refCol < len(fields) && rec.Ref == fields[g.refCol]) {
				continue
			}
		}
		found = true
		for i := range g.rules {
			r := &g.rules[i]
			if r.col >= 0 && r.col < len(fields) && fields[r.col] != "" {
				vals[i] = append(vals[i], fields[r.col])
			}
		}
	}

	if !found {
		return nil
	}
	for i := range g.rules {
		r := &g.rules[i]
		if r.col < 0 { // flag
			rec.AddInfoFlag(r.name)
			continue
		}
		if len(vals[i]) == 0 {
			continue
		}
		out, err := aggregateVals(vals[i], r.collapse, r.first, r.max)
		if err != nil {
			return err
		}
		if r.sample != "" {
			if err := rec.AddFormat(r.sampleIdx, r.name, out); err != nil {
				return err
			}
			continue
		}
		rec.AddInfo(r.name, out)
	}
	return nil
}

// Close releases the shared tabix reader.
func (g *TabixAnnotationGroup) Close() error { return g.reader.Close() }

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
