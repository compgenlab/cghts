package annotate

import (
	"fmt"
	"sort"
	"strings"

	"github.com/compgenlab/hts/vcf"
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
	suffix := vcfModifierSuffix(a.opts)
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

func vcfModifierSuffix(o VcfOptions) string {
	var parts []string
	if o.Passing {
		parts = append(parts, "passing")
	}
	if o.Exact {
		parts = append(parts, "exact")
	}
	if o.Unique {
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
