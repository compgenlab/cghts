package annotate

import (
	"fmt"
	"sort"
	"strings"

	"github.com/compgenlab/cghts/vcf"
)

// VcfOptions configures a [VcfAnnotation]. The source file must be a
// BGZF-compressed, tabix/CSI-indexed VCF.
type VcfOptions struct {
	Name     string // INFO key to add (or "@ID" to copy the ID column)
	Filename string // tabix-indexed source VCF
	Field    string // source INFO field to copy; "" = presence flag; "@ID" = copy the source ID as the value
	Exact    bool   // require REF and an ALT allele to match (forced for @ID)
	Passing  bool   // only consider source records that pass filters
	Unique   bool   // de-duplicate (sorted) when multiple values match
	NoHeader bool   // do not add a ##INFO def

	// AutoConvert matches contig names across UCSC/Ensembl/NCBI naming (human
	// primary contigs 1-22,X,Y,MT) instead of requiring an exact-string match.
	AutoConvert bool
}

// VcfAnnotation annotates a record from a tabix-indexed source VCF: it adds a
// flag, copies an INFO field, or copies the ID, matching by position (and
// optionally REF+ALT). It ports ngsutilsj VCFAnnotation (--vcf/--vcf-flag/
// --vcf-id).
type VcfAnnotation struct {
	base
	opts   VcfOptions
	reader *vcf.IndexedVcfReader
	isID   bool                 // Name == "@ID": copy the source ID into the record ID
	conv   *vcf.ContigConverter // non-nil when contig-name matching is enabled
}

// NewVcfAnnotation opens the source VCF and returns the annotator.
func NewVcfAnnotation(opts VcfOptions) (*VcfAnnotation, error) {
	r, err := vcf.NewIndexedVcfReader(opts.Filename)
	if err != nil {
		return nil, fmt.Errorf("annotate: open %s: %w", opts.Filename, err)
	}
	isID := opts.Name == "@ID"
	if isID {
		opts.Exact = true // ID copy is exact-match only
	}
	a := &VcfAnnotation{opts: opts, reader: r, isID: isID}
	if opts.AutoConvert {
		a.EnableContigMatching()
	}
	return a, nil
}

// EnableContigMatching turns on cross-scheme contig-name matching (UCSC/Ensembl/
// NCBI) by building a converter from the source VCF's contig names. It
// implements [ContigMatcher] so the CLI's --auto-convert flag can enable it.
func (a *VcfAnnotation) EnableContigMatching() {
	a.conv = vcf.NewContigConverter(a.reader.RefNames())
}

// SetupHeader adds the ##INFO def (none for @ID or NoHeader).
func (a *VcfAnnotation) SetupHeader(h *vcf.VcfHeader) error {
	if a.opts.NoHeader || a.isID {
		return nil
	}
	suffix := vcfModifierSuffix(a.opts.Passing, a.opts.Exact, a.opts.Unique)
	if a.opts.Field == "" { // flag
		h.AddInfo(infoDefSrc(a.opts.Name, "0", "Flag", "Present in VCF file"+suffix, a.opts.Filename))
	} else {
		h.AddInfo(infoDefSrc(a.opts.Name, "1", "String", a.opts.Field+" from VCF file"+suffix, a.opts.Filename))
	}
	return nil
}

// Annotate queries the source VCF at the variant position and applies the
// annotation. ID-copy and flag are simple: on the first matching source record
// they set the value and return.
func (a *VcfAnnotation) Annotate(rec *vcf.VcfRecord) error {
	chrom, ok := a.Chrom(rec)
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
	// Exact-position query (no deletion shift): a VCF:VCF comparison matches the
	// variant position exactly.
	seq, err := a.reader.Query(chrom, rec.Pos-1, rec.Pos)
	if err != nil {
		return err
	}

	var vals []string
	for src, err := range seq {
		if err != nil {
			return err
		}
		if src.Pos != rec.Pos {
			continue
		}
		if a.opts.Passing && src.IsFiltered() {
			continue
		}
		if a.opts.Exact && !altRefMatch(src, rec) {
			continue
		}

		if a.isID { // --vcf-id: copy the source ID and finish
			rec.SetID(src.ID())
			return nil
		}
		if a.opts.Field == "" { // --vcf-flag: mark present and finish
			rec.AddInfoFlag(a.opts.Name)
			return nil
		}
		if a.opts.Field == "@ID" {
			if id := src.ID(); id != "" {
				vals = append(vals, id)
			}
		} else if v, ok := src.InfoValue(a.opts.Field); ok {
			if s := v.String(); s != "" {
				vals = append(vals, s)
			}
		}
	}

	if len(vals) > 0 {
		if a.opts.Unique {
			vals = sortedUnique(vals)
		}
		rec.AddInfo(a.opts.Name, strings.Join(vals, ","))
	}
	return nil
}

// Close releases the source reader.
func (a *VcfAnnotation) Close() error { return a.reader.Close() }

// VcfFieldOptions is one field rule within a [VcfAnnotationGroup]: the same per-field
// knobs as [VcfOptions] minus the shared source file and contig-matching.
type VcfFieldOptions struct {
	Name     string // INFO key to add (or "@ID" to copy the ID column)
	Field    string // source INFO field to copy; "" = presence flag; "@ID" = copy the source ID as the value
	Exact    bool   // require REF and an ALT allele to match (forced for @ID)
	Passing  bool   // only consider source records that pass filters
	Unique   bool   // de-duplicate (sorted) when multiple values match
	NoHeader bool   // do not add a ##INFO def
}

// VcfGroupOptions configures a [VcfAnnotationGroup]: one tabix-indexed source VCF and
// N field rules that share a single reader and one query per record.
type VcfGroupOptions struct {
	Filename    string
	AutoConvert bool // cross-scheme contig-name matching (see [VcfOptions.AutoConvert])
	Fields      []VcfFieldOptions
}

// vcfRule is a resolved [VcfFieldOptions] (isID pre-computed).
type vcfRule struct {
	name     string
	field    string
	exact    bool
	passing  bool
	unique   bool
	noHeader bool
	isID     bool
}

// VcfAnnotationGroup annotates a record with several fields from one source VCF using
// a single reader and a single positional query per input record — the multi-field
// analogue of [VcfAnnotation]. It exists so a source referenced by N annotations opens
// one reader and scans/parses the region once instead of N times.
type VcfAnnotationGroup struct {
	base
	filename string
	reader   *vcf.IndexedVcfReader
	conv     *vcf.ContigConverter // non-nil when contig-name matching is enabled
	rules    []vcfRule
}

// NewVcfAnnotationGroup opens the source VCF once and returns a grouped annotator for
// all the given field rules.
func NewVcfAnnotationGroup(opts VcfGroupOptions) (*VcfAnnotationGroup, error) {
	r, err := vcf.NewIndexedVcfReader(opts.Filename)
	if err != nil {
		return nil, fmt.Errorf("annotate: open %s: %w", opts.Filename, err)
	}
	g := &VcfAnnotationGroup{filename: opts.Filename, reader: r}
	for _, f := range opts.Fields {
		isID := f.Name == "@ID"
		if isID {
			f.Exact = true // ID copy is exact-match only
		}
		g.rules = append(g.rules, vcfRule{
			name: f.Name, field: f.Field, exact: f.Exact, passing: f.Passing,
			unique: f.Unique, noHeader: f.NoHeader, isID: isID,
		})
	}
	if opts.AutoConvert {
		g.EnableContigMatching()
	}
	return g, nil
}

// EnableContigMatching turns on cross-scheme contig-name matching for the group,
// mirroring [VcfAnnotation.EnableContigMatching].
func (g *VcfAnnotationGroup) EnableContigMatching() {
	g.conv = vcf.NewContigConverter(g.reader.RefNames())
}

// SetupHeader adds one ##INFO def per rule (none for @ID or NoHeader).
func (g *VcfAnnotationGroup) SetupHeader(h *vcf.VcfHeader) error {
	for i := range g.rules {
		r := &g.rules[i]
		if r.noHeader || r.isID {
			continue
		}
		suffix := vcfModifierSuffix(r.passing, r.exact, r.unique)
		if r.field == "" {
			h.AddInfo(infoDefSrc(r.name, "0", "Flag", "Present in VCF file"+suffix, g.filename))
		} else {
			h.AddInfo(infoDefSrc(r.name, "1", "String", r.field+" from VCF file"+suffix, g.filename))
		}
	}
	return nil
}

// Annotate queries the source once and applies every rule, computing the per-source
// predicates (filtered, REF/ALT match) once per source record. It matches
// [VcfAnnotation.Annotate] field-for-field; the only difference is that @ID/flag rules
// set-once-and-continue (so co-grouped value fields still accumulate) rather than
// returning on the first match.
func (g *VcfAnnotationGroup) Annotate(rec *vcf.VcfRecord) error {
	chrom, ok := g.Chrom(rec)
	if !ok {
		return nil
	}
	if g.conv != nil {
		if chrom, ok = g.conv.Resolve(chrom); !ok {
			return nil
		}
	} else if !g.reader.HasRef(chrom) {
		return nil
	}
	seq, err := g.reader.Query(chrom, rec.Pos-1, rec.Pos)
	if err != nil {
		return err
	}

	vals := make([][]string, len(g.rules))
	done := make([]bool, len(g.rules)) // @ID / flag rules: set once

	for src, err := range seq {
		if err != nil {
			return err
		}
		if src.Pos != rec.Pos {
			continue
		}
		filtered := src.IsFiltered()
		arMatch := altRefMatch(src, rec)
		for i := range g.rules {
			r := &g.rules[i]
			if done[i] {
				continue
			}
			if r.passing && filtered {
				continue
			}
			if r.exact && !arMatch {
				continue
			}
			if r.isID {
				rec.SetID(src.ID())
				done[i] = true
				continue
			}
			if r.field == "" { // flag
				rec.AddInfoFlag(r.name)
				done[i] = true
				continue
			}
			if r.field == "@ID" {
				if id := src.ID(); id != "" {
					vals[i] = append(vals[i], id)
				}
			} else if v, ok := src.InfoValue(r.field); ok {
				if s := v.String(); s != "" {
					vals[i] = append(vals[i], s)
				}
			}
		}
	}

	for i := range g.rules {
		r := &g.rules[i]
		if len(vals[i]) == 0 {
			continue
		}
		v := vals[i]
		if r.unique {
			v = sortedUnique(v)
		}
		rec.AddInfo(r.name, strings.Join(v, ","))
	}
	return nil
}

// Close releases the shared source reader.
func (g *VcfAnnotationGroup) Close() error { return g.reader.Close() }

// altRefMatch reports whether the source and target records share REF and at
// least one ALT allele.
func altRefMatch(src, rec *vcf.VcfRecord) bool {
	if src.Ref != rec.Ref {
		return false
	}
	for _, a1 := range src.Alt() {
		for _, a2 := range rec.Alt() {
			if a1 == a2 {
				return true
			}
		}
	}
	return false
}

func vcfModifierSuffix(passing, exact, unique bool) string {
	var parts []string
	if passing {
		parts = append(parts, "passing")
	}
	if exact {
		parts = append(parts, "exact")
	}
	if unique {
		parts = append(parts, "unique")
	}
	if len(parts) == 0 {
		return ""
	}
	return " (" + strings.Join(parts, ",") + ")"
}

func sortedUnique(vals []string) []string {
	set := make(map[string]bool, len(vals))
	for _, v := range vals {
		set[v] = true
	}
	out := make([]string, 0, len(set))
	for v := range set {
		out = append(out, v)
	}
	sort.Strings(out)
	return out
}
