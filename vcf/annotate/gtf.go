package annotate

import (
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/compgenlab/cghts/bed"
	"github.com/compgenlab/cghts/gtf"
	"github.com/compgenlab/cghts/vcf"
)

// GeneModel is the query surface a [GtfAnnotator] needs from a GTF.
//
// Both of the gtf package's sources satisfy it — the whole-file
// [gtf.AnnotationSource] and the tabix-backed [gtf.IndexedAnnotationSource] —
// which is what lets the annotator hold either. That mattered: loading a GTF
// fully into memory is fine for one small file and is not fine for a human
// GENCODE, and a caller annotating a whole genome against two of them ran out
// of memory. The choice belongs here rather than in each caller, and a caller
// that worked around it by writing its own annotator ended up with a second
// implementation of the region logic.
type GeneModel interface {
	FindGenes(ref string, start, end int) []*gtf.Gene
	FindGenicRegionForPos(ref string, pos int, strand bed.Strand, geneID string) gtf.GenicRegion
	RefNames() []string
	HasRef(ref string) bool
}

// biotypeProvider is implemented by a gene model that can say up front whether
// its GTF carries biotypes.
//
// Only the whole-file model can: it has read everything. The indexed one has
// read nothing until it is queried, and scanning the file to find out would cost
// exactly what indexing it saved.
type biotypeProvider interface {
	Provides(key string) bool
}

// GtfOptions configures a [GtfAnnotator].
//
// By default the GTF is queried through its tabix index when one sits beside it
// and read fully into memory when none does — see NewGtfAnnotator.
type GtfOptions struct {
	Prefix string // INFO key prefix; defaults to "GTF_"

	// Names overrides the INFO key each logical field is written under, keyed by
	// the GtfGeneSymbol/GtfGeneID/… constants. A named field ignores Prefix —
	// the caller has given the whole key, not a suffix to decorate.
	Names FieldNames

	Filename string // GTF file (optionally .gz)

	// Source, when set, is used instead of opening Filename — so the gene model
	// can be one the caller already built, indexed or not, or one served from
	// somewhere that is not a local file. Filename is then only a label,
	// recorded as the annotation's provenance in the header.
	//
	// Ownership transfers: Close releases it if it holds anything.
	Source GeneModel

	// Fields selects which logical fields to emit, from the GtfGeneSymbol/
	// GtfGeneID/… constants. Empty emits all of them, which is what a caller
	// asking for "the GTF annotations" means.
	//
	// Selecting matters because these seven are one annotator's output but not
	// one annotation: a caller whose configuration names two of them should get
	// two INFO fields, not seven with five it never asked for cluttering every
	// record of a whole-genome VCF.
	Fields []string

	// InMemory forces the whole-file gene model even when an index is present.
	// The indexed reader is not safe for concurrent use, so a caller sharing one
	// annotator across goroutines wants this; everything else wants the default.
	InMemory     bool
	RequiredTags []string // keep only features carrying every tag (the --gtf-tag filter)

	// AutoConvert matches contig names across UCSC/Ensembl/NCBI naming (human
	// primary contigs 1-22,X,Y,MT) instead of requiring an exact-string match.
	AutoConvert bool
}

// GtfAnnotator overlays gene annotations from a GTF onto VCF records: for every
// gene overlapping a variant it writes the gene name(s), ID(s), strand(s),
// biotype(s), and a genic-region classification (coding_exon / UTR / intron /
// nc_exon / …), plus the coding and non-coding gene names split out. It ports
// ngsutilsj's vcf-annotate --gtf (the GTFGene annotator).
//
// INFO fields added (default prefix GTF_):
//
//	GTF_GENE      gene name(s)
//	GTF_GENEID    gene ID(s)
//	GTF_STRAND    gene strand(s)
//	GTF_BIOTYPE   gene biotype(s)        (only when the GTF supplies biotypes)
//	GTF_REGION    genic region code(s)
//	GTF_CODING    name(s) of overlapping coding genes      (only when present)
//	GTF_NONCODING name(s) of overlapping non-coding genes  (only when present)
//
// Multiple overlapping genes are comma-joined in parallel across the fields.
type GtfAnnotator struct {
	base
	opts   GtfOptions
	prefix string
	src    GeneModel
	conv   *vcf.ContigConverter // non-nil when contig-name matching is enabled
}

// NewGtfAnnotator opens the GTF's gene model and returns the annotator.
//
// The tabix index is preferred when one sits beside the file, because it bounds
// memory to the genes actually queried; without one the whole GTF is read in,
// which for a human annotation set is over a gigabyte. Index it with a GFF-preset
// tabix index to get the cheap path — bgzip the GTF and `tabix -p gff` it.
func NewGtfAnnotator(opts GtfOptions) (*GtfAnnotator, error) {
	src := opts.Source
	if src == nil {
		var err error
		if src, err = openGeneModel(opts); err != nil {
			return nil, fmt.Errorf("annotate: %w", err)
		}
	}
	prefix := opts.Prefix
	if prefix == "" {
		prefix = "GTF_"
	}
	a := &GtfAnnotator{opts: opts, prefix: prefix, src: src}
	if opts.AutoConvert {
		a.EnableContigMatching()
	}
	return a, nil
}

// EnableContigMatching turns on cross-scheme contig-name matching (UCSC/Ensembl/
// NCBI) using the GTF's contig names. Implements [ContigMatcher] for
// --auto-convert.
func (a *GtfAnnotator) EnableContigMatching() {
	a.conv = vcf.NewContigConverter(a.src.RefNames())
}

// SetupHeader adds the ##INFO defs. CG_BIOTYPE is declared only when the GTF
// actually supplies biotypes; CG_CODING/CG_NONCODING are always declared (they
// are written per-record when an overlapping gene of that kind exists).
func (a *GtfAnnotator) SetupHeader(h *vcf.VcfHeader) error {
	for _, f := range []struct{ logical, desc string }{
		{GtfGeneSymbol, "Gene name"},
		{GtfGeneID, "Gene ID"},
		{GtfStrand, "Gene strand"},
		{GtfBiotype, "Gene biotype"},
		{GtfRegion, "Genic region"},
		{GtfCoding, "Coding gene name"},
		{GtfNoncoding, "Non-coding gene name"},
	} {
		if !a.emits(f.logical) {
			continue
		}
		if f.logical == GtfBiotype && !a.providesBiotype() {
			continue
		}
		h.AddInfo(infoDefSrc(a.key(f.logical), ".", "String", f.desc, a.opts.Filename))
	}
	return nil
}

// emits reports whether a logical field is one this annotator writes.
func (a *GtfAnnotator) emits(logical string) bool {
	if len(a.opts.Fields) == 0 {
		return true
	}
	for _, f := range a.opts.Fields {
		if strings.EqualFold(f, logical) {
			return true
		}
	}
	return false
}

// providesBiotype reports whether the BIOTYPE field should be declared.
//
// True when the model cannot say, which is the safe direction: a declared field
// that is never written is a header line nobody reads, while a field written
// without one is a record that a strict parser is entitled to reject.
func (a *GtfAnnotator) providesBiotype() bool {
	if p, ok := a.src.(biotypeProvider); ok {
		return p.Provides("biotype")
	}
	return true
}

// SetFieldNames chooses the INFO key each logical field is written under.
func (a *GtfAnnotator) SetFieldNames(n FieldNames) { a.opts.Names = n }

// key is the INFO id one logical field is written under: the caller's name for
// it, or the prefix and the historical suffix.
func (a *GtfAnnotator) key(logical string) string {
	return a.opts.Names.nameOr(logical, a.prefix+logical)
}

// Annotate finds the genes overlapping the variant position and writes the
// gene/region INFO fields. Variants are unstranded, so regions are always sense
// codes (matching GTFGene.annotate).
func (a *GtfAnnotator) Annotate(rec *vcf.VcfRecord) error {
	chrom, ok := a.Chrom(rec)
	if !ok {
		return nil
	}
	pos, ok := a.Pos(rec)
	if !ok {
		return nil
	}
	pos0 := pos - 1 // 1-based → 0-based

	if a.conv != nil {
		if chrom, ok = a.conv.Resolve(chrom); !ok {
			return nil
		}
	} else if !a.src.HasRef(chrom) {
		return nil
	}

	genes := a.src.FindGenes(chrom, pos0, pos0+1)
	if len(genes) == 0 {
		return nil
	}

	var names, ids, strands, biotypes, regions, coding, noncoding []string
	hasBiotype := false
	for _, g := range genes {
		names = append(names, g.GeneName)
		ids = append(ids, g.GeneID)
		strands = append(strands, string(g.Strand))
		if g.BioType != "" {
			biotypes = append(biotypes, g.BioType)
			hasBiotype = true
		} else {
			biotypes = append(biotypes, ".")
		}
		regions = append(regions, a.src.FindGenicRegionForPos(chrom, pos0, bed.StrandNone, g.GeneID).Code)
		if g.IsCoding() {
			coding = append(coding, g.GeneName)
		} else {
			noncoding = append(noncoding, g.GeneName)
		}
	}

	a.write(rec, GtfGeneSymbol, names)
	a.write(rec, GtfGeneID, ids)
	a.write(rec, GtfStrand, strands)
	if hasBiotype {
		a.write(rec, GtfBiotype, biotypes)
	}
	a.write(rec, GtfRegion, regions)
	a.write(rec, GtfCoding, coding)
	a.write(rec, GtfNoncoding, noncoding)
	return nil
}

// Close is a no-op: the gene model lives in memory.
func (a *GtfAnnotator) Close() error {
	// The indexed model holds a tabix reader; the in-memory one holds nothing.
	if c, ok := a.src.(io.Closer); ok {
		return c.Close()
	}
	return nil
}

// openGeneModel chooses how to read the GTF: through its index when it has one,
// wholly into memory when it does not.
//
// Falling back rather than failing, because an unindexed GTF is a working GTF —
// it just costs memory proportional to the file instead of to the query. A small
// annotation set is entirely reasonable to load, and refusing one because it has
// no .tbi would break every caller that has been passing plain files.
func openGeneModel(opts GtfOptions) (GeneModel, error) {
	if !opts.InMemory && hasTabixIndex(opts.Filename) {
		return gtf.NewIndexedAnnotationSource(opts.Filename, opts.RequiredTags)
	}
	return gtf.NewAnnotationSource(opts.Filename, opts.RequiredTags)
}

// hasTabixIndex reports whether a tabix index sits beside the file.
//
// Checked by looking rather than by trying and falling back: opening a bgzipped
// GTF without an index succeeds and then fails on the first query, which would
// surface as "no genes found" rather than as a missing index.
func hasTabixIndex(path string) bool {
	for _, ext := range []string{".tbi", ".csi"} {
		if st, err := os.Stat(path + ext); err == nil && !st.IsDir() {
			return true
		}
	}
	return false
}

// write adds one field, if it is one this annotator emits and there is anything
// to say. Overlapping genes are comma-joined in parallel across the fields, so
// the nth value of each belongs to the nth gene.
func (a *GtfAnnotator) write(rec *vcf.VcfRecord, logical string, vals []string) {
	if len(vals) == 0 || !a.emits(logical) {
		return
	}
	rec.AddInfo(a.key(logical), strings.Join(vals, ","))
}
