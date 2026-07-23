// Package gtf parses a GTF gene-annotation file into an in-memory gene model
// (genes → transcripts → exons/CDS/codons) with an interval index for position
// lookup, and classifies genomic positions into genic regions (coding exon,
// UTR, intron, junction, …). It is a Go port of ngsutilsj's GTFAnnotationSource
// / GenicRegion, reproducing the same biotype derivation and region-code
// precedence so results match.
//
// The model is reusable beyond VCF annotation: the position/region classifiers
// and gene iteration are the basis for read-level region counting (e.g. RNA-seq
// gene coverage in a future sam-stats / BAM read counter). Coordinates are
// 0-based half-open throughout; the GTF file's 1-based inclusive coordinates are
// converted on parse (start-1, end unchanged).
package gtf

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/compgenlab/cghts/bed"
)

// binSize is the fixed width of the interval index bins (ports
// AbstractAnnotationSource.BINSIZE).
const binSize = 100_000

type refBin struct {
	ref string
	bin int
}

// AnnotationSource is an in-memory GTF gene model with a bin interval index. It
// ports GTFAnnotationSource together with its AbstractAnnotationSource bin index.
// Build one with NewAnnotationSource.
type AnnotationSource struct {
	genes      []*Gene
	bins       map[refBin][]*Gene
	refNames   []string
	hasBioType bool
	hasStatus  bool
}

// NewAnnotationSource reads a GTF file (optionally gzip/bgzip-compressed, by a
// ".gz" suffix) into an in-memory gene model. When requiredTags is non-empty,
// only feature rows carrying every listed `tag` attribute are kept (AND
// semantics) — the ngsutilsj --gtf-tag filter.
func NewAnnotationSource(path string, requiredTags []string) (*AnnotationSource, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("gtf: open %s: %w", path, err)
	}
	defer f.Close()

	var r io.Reader = f
	if strings.HasSuffix(path, ".gz") {
		gz, gzErr := gzip.NewReader(bufio.NewReader(f))
		if gzErr != nil {
			return nil, fmt.Errorf("gtf: gzip %s: %w", path, gzErr)
		}
		defer gz.Close()
		r = gz
	}

	src, err := parse(r, requiredTags)
	if err != nil {
		return nil, fmt.Errorf("gtf: read %s: %w", path, err)
	}
	return src, nil
}

func parse(r io.Reader, requiredTags []string) (*AnnotationSource, error) {
	src := &AnnotationSource{bins: make(map[refBin][]*Gene)}

	// Genes are keyed by (ref, gene_id) so a gene_id reused on another
	// chromosome (some miRNAs) stays a distinct gene — the same outcome as the
	// Java per-chromosome cache flush, but without assuming chromosome-grouped
	// input.
	cache := make(map[string]*Gene)

	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
	for sc.Scan() {
		line := sc.Text()
		if len(line) == 0 || line[0] == '#' {
			continue
		}
		cols := strings.Split(line, "\t")
		if len(cols) < 9 {
			continue
		}
		chrom := cols[0]
		recordType := cols[2]
		start, err := strconv.Atoi(cols[3])
		if err != nil {
			return nil, fmt.Errorf("bad start %q: %w", cols[3], err)
		}
		start-- // 1-based inclusive → 0-based half-open
		end, err := strconv.Atoi(cols[4])
		if err != nil {
			return nil, fmt.Errorf("bad end %q: %w", cols[4], err)
		}
		strand := bed.ParseStrand(cols[6])

		a := parseAttributes(cols[8])
		if a.geneName == "" && a.gene != "" {
			a.geneName = a.gene // RefSeq GTFs use "gene" rather than "gene_name"
		}
		if len(requiredTags) > 0 && !hasAllTags(a.tags, requiredTags) {
			continue
		}
		if !src.hasBioType && a.bioType != "" {
			src.hasBioType = true
		}
		if !src.hasStatus && a.status != "" {
			src.hasStatus = true
		}

		key := chrom + "\t" + a.geneID
		gene, ok := cache[key]
		if !ok {
			// Seed the gene from this row so genes that have no separate
			// transcript/exon lines (some RefSeq genes) still get a span.
			gene = newGene(a.geneID, a.geneName, chrom, start, end, strand, a.bioType, a.status)
			cache[key] = gene
			src.genes = append(src.genes, gene)
		} else {
			// Backfill gene-level attributes that only appear on some rows.
			// RefSeq carries gene_biotype on the "gene" feature line only, while
			// transcript/exon rows omit it; if such a row seeded the gene first,
			// the biotype (and name/status) would otherwise be lost regardless of
			// row order. GENCODE repeats these on every row, so it is unaffected.
			if gene.BioType == "" && a.bioType != "" {
				gene.BioType = a.bioType
			}
			if gene.Status == "" && a.status != "" {
				gene.Status = a.status
			}
			if gene.GeneName == "" && a.geneName != "" {
				gene.GeneName = a.geneName
			}
		}

		switch recordType {
		case "exon":
			gene.addExon(a.transcriptID, start, end, a.attrs)
		case "CDS":
			gene.addCDS(a.transcriptID, start, end, a.attrs)
		case "stop_codon":
			gene.addStopCodon(a.transcriptID, start, end, a.attrs)
		case "start_codon":
			gene.addStartCodon(a.transcriptID, start, end, a.attrs)
		}
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}

	src.buildIndex()
	return src, nil
}

// attributes holds the recognized fields parsed out of a GTF attribute column.
type attributes struct {
	geneID       string
	geneName     string
	gene         string // RefSeq "gene" key
	transcriptID string
	bioType      string
	status       string
	attrs        [][2]string // leftover attributes, file order
	tags         []string    // values of "tag" attributes
}

// parseAttributes parses the GTF 9th column. Only gene_type/gene_biotype feed
// the biotype (no transcript_biotype/source-column fallback), matching the Java
// implementation.
func parseAttributes(col string) attributes {
	var a attributes
	for _, attr := range splitQuoted(col, ';') {
		attr = strings.TrimSpace(attr)
		if attr == "" {
			continue
		}
		k, v, _ := strings.Cut(attr, " ")
		v = strings.Trim(strings.TrimSpace(v), "\"")
		switch k {
		case "gene_id":
			a.geneID = v
		case "gene_name":
			a.geneName = v
		case "gene":
			a.gene = v
		case "transcript_id":
			a.transcriptID = v
		case "gene_type", "gene_biotype":
			a.bioType = v
		case "gene_status":
			a.status = v
		default:
			a.attrs = append(a.attrs, [2]string{k, v})
			if k == "tag" {
				a.tags = append(a.tags, v)
			}
		}
	}
	return a
}

// splitQuoted splits s on sep, ignoring separators inside double quotes. Ports
// StringUtils.quotedSplit.
func splitQuoted(s string, sep byte) []string {
	var out []string
	var b strings.Builder
	inQuote := false
	for i := 0; i < len(s); i++ {
		c := s[i]
		switch {
		case c == '"':
			inQuote = !inQuote
			b.WriteByte(c)
		case c == sep && !inQuote:
			out = append(out, b.String())
			b.Reset()
		default:
			b.WriteByte(c)
		}
	}
	if b.Len() > 0 {
		out = append(out, b.String())
	}
	return out
}

func hasAllTags(have, required []string) bool {
	for _, r := range required {
		found := false
		for _, h := range have {
			if h == r {
				found = true
				break
			}
		}
		if !found {
			return false
		}
	}
	return true
}

// buildIndex registers every gene in each interval bin it overlaps and sorts
// each bin's slice by start so lookups can early-break. Ports
// AbstractAnnotationSource.addAnnotation (done once, after parsing).
func (s *AnnotationSource) buildIndex() {
	refSeen := make(map[string]bool)
	for _, g := range s.genes {
		if !refSeen[g.Ref] {
			refSeen[g.Ref] = true
			s.refNames = append(s.refNames, g.Ref)
		}
		lo, hi := binsOf(g.Start, g.End)
		for b := lo; b <= hi; b++ {
			rb := refBin{g.Ref, b}
			s.bins[rb] = append(s.bins[rb], g)
		}
	}
	for rb := range s.bins {
		genes := s.bins[rb]
		sort.Slice(genes, func(i, j int) bool {
			if genes[i].Start != genes[j].Start {
				return genes[i].Start < genes[j].Start
			}
			return genes[i].End < genes[j].End
		})
	}
}

// binsOf returns the inclusive bin range covering [start, end). Ports
// RefBin.getBins (which divides the half-open end by binSize, inclusive).
func binsOf(start, end int) (int, int) {
	return start / binSize, end / binSize
}

// FindGenes returns the genes overlapping [start, end) on ref (0-based
// half-open), de-duplicated and sorted by position. Ports findAnnotation
// (overlap mode).
func (s *AnnotationSource) FindGenes(ref string, start, end int) []*Gene {
	q := span{ref: ref, start: start, end: end, strand: bed.StrandNone}
	seen := make(map[*Gene]bool)
	var out []*Gene
	lo, hi := binsOf(start, end)
	for b := lo; b <= hi; b++ {
		for _, g := range s.bins[refBin{ref, b}] {
			if g.Start > end {
				break // bin slice is start-sorted; nothing further can overlap
			}
			if seen[g] {
				continue
			}
			if g.coord().overlaps(q) {
				seen[g] = true
				out = append(out, g)
			}
		}
	}
	sort.Slice(out, func(i, j int) bool {
		if out[i].Start != out[j].Start {
			return out[i].Start < out[j].Start
		}
		return out[i].GeneID < out[j].GeneID
	})
	return out
}

// Genes returns all genes, ordered by chromosome (first-seen order) then
// position. Useful for whole-genome iteration (gene/exon emission, counting).
func (s *AnnotationSource) Genes() []*Gene {
	refOrder := make(map[string]int, len(s.refNames))
	for i, r := range s.refNames {
		refOrder[r] = i
	}
	out := make([]*Gene, len(s.genes))
	copy(out, s.genes)
	sort.SliceStable(out, func(i, j int) bool {
		if out[i].Ref != out[j].Ref {
			return refOrder[out[i].Ref] < refOrder[out[j].Ref]
		}
		if out[i].Start != out[j].Start {
			return out[i].Start < out[j].Start
		}
		return out[i].GeneID < out[j].GeneID
	})
	return out
}

// RefNames returns the chromosomes present in the GTF, in first-seen order.
// Suitable for building a vcf.ContigConverter.
func (s *AnnotationSource) RefNames() []string { return s.refNames }

// HasRef reports whether the GTF contains any gene on ref.
func (s *AnnotationSource) HasRef(ref string) bool {
	for _, r := range s.refNames {
		if r == ref {
			return true
		}
	}
	return false
}

// Size returns the number of genes parsed.
func (s *AnnotationSource) Size() int { return len(s.genes) }

// Provides reports whether the GTF supplied a given annotation field. Recognized
// keys: "biotype", "status", "gene_id", "gene_name", "start", "end", "strand".
// Ports AbstractAnnotationSource.provides + getAnnotationNames.
func (s *AnnotationSource) Provides(key string) bool {
	switch key {
	case "biotype":
		return s.hasBioType
	case "status":
		return s.hasStatus
	case "gene_id", "gene_name", "start", "end", "strand":
		return true
	}
	return false
}
