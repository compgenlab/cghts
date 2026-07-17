package annotate

import (
	"fmt"
	"strings"

	"github.com/compgenlab/cghts/bed"
	"github.com/compgenlab/cghts/gtf"
	"github.com/compgenlab/cghts/vcf"
)

// GtfOptions configures a [GtfAnnotator]. The GTF file is read fully into memory
// (a gene model with an interval index); it may be gzip/bgzip-compressed.
type GtfOptions struct {
	Prefix       string   // INFO key prefix; defaults to "GTF_"
	Filename     string   // GTF file (optionally .gz)
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
	src    *gtf.AnnotationSource
	conv   *vcf.ContigConverter // non-nil when contig-name matching is enabled
}

// NewGtfAnnotator loads the GTF into memory and returns the annotator.
func NewGtfAnnotator(opts GtfOptions) (*GtfAnnotator, error) {
	src, err := gtf.NewAnnotationSource(opts.Filename, opts.RequiredTags)
	if err != nil {
		return nil, fmt.Errorf("annotate: %w", err)
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
	p := a.prefix
	h.AddInfo(infoDefSrc(p+"GENE", ".", "String", "Gene name", a.opts.Filename))
	h.AddInfo(infoDefSrc(p+"GENEID", ".", "String", "Gene ID", a.opts.Filename))
	h.AddInfo(infoDefSrc(p+"STRAND", ".", "String", "Gene strand", a.opts.Filename))
	if a.src.Provides("biotype") {
		h.AddInfo(infoDefSrc(p+"BIOTYPE", ".", "String", "Gene biotype", a.opts.Filename))
	}
	h.AddInfo(infoDefSrc(p+"REGION", ".", "String", "Genic region", a.opts.Filename))
	h.AddInfo(infoDefSrc(p+"CODING", ".", "String", "Coding gene name", a.opts.Filename))
	h.AddInfo(infoDefSrc(p+"NONCODING", ".", "String", "Non-coding gene name", a.opts.Filename))
	return nil
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

	p := a.prefix
	rec.AddInfo(p+"GENE", strings.Join(names, ","))
	rec.AddInfo(p+"GENEID", strings.Join(ids, ","))
	rec.AddInfo(p+"STRAND", strings.Join(strands, ","))
	if hasBiotype {
		rec.AddInfo(p+"BIOTYPE", strings.Join(biotypes, ","))
	}
	rec.AddInfo(p+"REGION", strings.Join(regions, ","))
	if len(coding) > 0 {
		rec.AddInfo(p+"CODING", strings.Join(coding, ","))
	}
	if len(noncoding) > 0 {
		rec.AddInfo(p+"NONCODING", strings.Join(noncoding, ","))
	}
	return nil
}

// Close is a no-op: the gene model lives in memory.
func (a *GtfAnnotator) Close() error { return nil }
