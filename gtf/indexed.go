package gtf

import (
	"container/list"
	"fmt"
	"math"
	"sort"
	"strconv"
	"strings"

	"github.com/compgenlab/cghts/bed"
	"github.com/compgenlab/cghts/htsio/tabix"
)

// defaultGeneCacheCap bounds the per-gene model cache. A position-sorted VCF hits
// many consecutive variants in the same gene, so a small cache serves almost every
// query from memory; 100 comfortably covers dense overlapping-gene clusters (HLA).
const defaultGeneCacheCap = 100

// IndexedAnnotationSource is a tabix-backed drop-in for AnnotationSource that
// queries a bgzipped + tabix-indexed GTF per position instead of loading the whole
// file into memory. For each variant it fetches only the overlapping gene(s), builds
// their full gene models on demand (reusing parse), and caches them in an LRU keyed
// by (ref, gene_id). Peak memory is "a handful of genes," independent of GTF size.
//
// It exposes the same FindGenes / FindGenicRegionForPos / RefNames / HasRef surface
// as AnnotationSource, so callers can hold either behind an interface. Region
// classification is correct because FindGenicRegionForPos is gene-local and each
// cached model contains the complete gene (all its exon/CDS rows, read from the
// gene's full span).
//
// Not safe for concurrent use by multiple goroutines (the LRU is mutated on every
// query); build one per annotation pass, as the pipeline already does.
type IndexedAnnotationSource struct {
	tr   *tabix.Reader
	tags []string
	cap  int

	ll    *list.List // MRU at front; values are *geneCacheEntry
	cache map[geneKey]*list.Element

	// What a sample of the file says it carries, and whether it has been taken.
	// See Provides.
	probed     bool
	hasBioType bool
	hasStatus  bool
}

type geneKey struct{ ref, geneID string }

type geneCacheEntry struct {
	key geneKey
	src *AnnotationSource // a single-gene model (complete)
}

// NewIndexedAnnotationSource opens a bgzipped + tabix-indexed GTF (a .tbi/.csi must
// sit beside gtfPath — build it with a GFF-preset tabix index). requiredTags is the
// same --gtf-tag / gtf_tags filter as NewAnnotationSource, applied when a queried
// gene's rows are parsed.
func NewIndexedAnnotationSource(gtfPath string, requiredTags []string) (*IndexedAnnotationSource, error) {
	tr, err := tabix.NewReader(gtfPath)
	if err != nil {
		return nil, fmt.Errorf("gtf: open indexed %s: %w", gtfPath, err)
	}
	return NewIndexedAnnotationSourceFrom(tr, requiredTags), nil
}

// NewIndexedAnnotationSourceFrom wraps an already-open tabix reader, so the gene
// model can be queried from a source that is not a local file — see
// [tabix.NewReaderFromSource]. Ownership of tr transfers: Close releases it.
func NewIndexedAnnotationSourceFrom(tr *tabix.Reader, requiredTags []string) *IndexedAnnotationSource {
	return &IndexedAnnotationSource{
		tr:    tr,
		tags:  requiredTags,
		cap:   defaultGeneCacheCap,
		ll:    list.New(),
		cache: make(map[geneKey]*list.Element),
	}
}

// RefNames returns the chromosomes present in the index, mirroring
// AnnotationSource.RefNames (suitable for building a vcf.ContigConverter).
func (s *IndexedAnnotationSource) RefNames() []string { return s.tr.RefNames() }

// HasRef reports whether the index contains ref.
func (s *IndexedAnnotationSource) HasRef(ref string) bool { return s.tr.HasRef(ref) }

// Close releases the underlying tabix reader.
func (s *IndexedAnnotationSource) Close() error { return s.tr.Close() }

// FindGenes returns the genes overlapping [start, end) on ref, deduplicated and
// sorted by (Start, GeneID) — matching AnnotationSource.FindGenes.
func (s *IndexedAnnotationSource) FindGenes(ref string, start, end int) []*Gene {
	spans := s.overlappingSpans(ref, start, end)
	q := span{ref: ref, start: start, end: end, strand: bed.StrandNone}
	seen := make(map[string]bool)
	var out []*Gene
	for _, sp := range spans {
		g := s.gene(ref, sp)
		if g == nil || seen[g.GeneID] || !g.coord().overlaps(q) {
			continue
		}
		seen[g.GeneID] = true
		out = append(out, g)
	}
	sort.Slice(out, func(i, j int) bool {
		if out[i].Start != out[j].Start {
			return out[i].Start < out[j].Start
		}
		return out[i].GeneID < out[j].GeneID
	})
	return out
}

// FindGenicRegionForPos classifies pos within gene geneID, mirroring
// AnnotationSource.FindGenicRegionForPos. It delegates to the gene's cached
// single-gene model (self-building it from a point query if not already warm).
func (s *IndexedAnnotationSource) FindGenicRegionForPos(ref string, pos int, strand bed.Strand, geneID string) GenicRegion {
	if _, ok := s.cache[geneKey{ref, geneID}]; !ok {
		for _, sp := range s.overlappingSpans(ref, pos, pos+1) {
			if sp.geneID == geneID {
				s.gene(ref, sp) // build + cache
				break
			}
		}
	}
	el, ok := s.cache[geneKey{ref, geneID}]
	if !ok {
		return Intergenic
	}
	s.ll.MoveToFront(el)
	return el.Value.(*geneCacheEntry).src.FindGenicRegionForPos(ref, pos, strand, geneID)
}

// geneSpan is a gene's id and its full 0-based half-open extent, read from the
// gene's `gene` row (or, as a fallback, the union of its overlapping `transcript`
// rows) returned by a tabix query.
type geneSpan struct {
	geneID     string
	start, end int
}

// overlappingSpans queries [start, end) and returns, per overlapping gene_id, the
// union of every returned row's extent. Because a `gene` (or `transcript`) row
// spans the whole gene (introns included) and is returned by any point inside it,
// the union yields the gene's full span whenever such a row exists — the common
// GENCODE/Ensembl case. For GTFs with only exon/CDS rows the union covers the
// exon(s) overlapping the query, which is enough to classify an exonic position;
// intronic positions in such gene-row-less GTFs can't be resolved (a known limit).
func (s *IndexedAnnotationSource) overlappingSpans(ref string, start, end int) []geneSpan {
	seq, err := s.tr.Query(ref, start, end)
	if err != nil {
		return nil // unknown ref → no overlap (same outcome as an empty bin lookup)
	}
	byGene := make(map[string]geneSpan)
	var order []string
	for rec, err := range seq {
		if err != nil {
			return nil
		}
		cols := strings.Split(rec.Line, "\t")
		if len(cols) < 9 {
			continue
		}
		gs, ge, ok := spanOf(cols)
		if !ok {
			continue
		}
		gid := parseAttributes(cols[8]).geneID
		if gid == "" {
			continue
		}
		if u, ok := byGene[gid]; ok {
			if gs < u.start {
				u.start = gs
			}
			if ge > u.end {
				u.end = ge
			}
			byGene[gid] = u
		} else {
			byGene[gid] = geneSpan{gid, gs, ge}
			order = append(order, gid)
		}
	}
	out := make([]geneSpan, 0, len(order))
	for _, gid := range order {
		out = append(out, byGene[gid])
	}
	return out
}

// gene returns the complete model for one gene, from the LRU cache or by querying
// the gene's full span, parsing it, and caching a single-gene source.
func (s *IndexedAnnotationSource) gene(ref string, sp geneSpan) *Gene {
	key := geneKey{ref, sp.geneID}
	if el, ok := s.cache[key]; ok {
		s.ll.MoveToFront(el)
		return soleGene(el.Value.(*geneCacheEntry).src, sp.geneID)
	}

	seq, err := s.tr.Query(ref, sp.start, sp.end)
	if err != nil {
		return nil
	}
	var b strings.Builder
	for rec, err := range seq {
		if err != nil {
			return nil
		}
		b.WriteString(rec.Line)
		b.WriteByte('\n')
	}
	full, err := parse(strings.NewReader(b.String()), s.tags)
	if err != nil {
		return nil
	}
	g := soleGene(full, sp.geneID)
	if g == nil {
		return nil
	}
	// Cache just this gene (a single-gene source) — the span query may also pull
	// neighboring genes, but only this gene's full extent was queried, so only it is
	// guaranteed complete.
	s.putGene(key, newSourceFromGenes([]*Gene{g}))
	return g
}

func (s *IndexedAnnotationSource) putGene(key geneKey, src *AnnotationSource) {
	el := s.ll.PushFront(&geneCacheEntry{key: key, src: src})
	s.cache[key] = el
	for s.ll.Len() > s.cap {
		back := s.ll.Back()
		if back == nil {
			break
		}
		s.ll.Remove(back)
		delete(s.cache, back.Value.(*geneCacheEntry).key)
	}
}

// soleGene returns the gene with the given id from a source, or nil.
func soleGene(src *AnnotationSource, geneID string) *Gene {
	for _, g := range src.genes {
		if g.GeneID == geneID {
			return g
		}
	}
	return nil
}

// newSourceFromGenes builds a minimal AnnotationSource holding the given genes with
// a fresh bin index — enough for FindGenes / FindGenicRegionForPos on a per-gene model.
func newSourceFromGenes(genes []*Gene) *AnnotationSource {
	s := &AnnotationSource{bins: make(map[refBin][]*Gene), genes: genes}
	s.buildIndex()
	return s
}

// spanOf extracts a GTF row's 0-based half-open [start, end) from cols 4/5.
func spanOf(cols []string) (int, int, bool) {
	start, err := strconv.Atoi(cols[3])
	if err != nil {
		return 0, 0, false
	}
	end, err := strconv.Atoi(cols[4])
	if err != nil {
		return 0, 0, false
	}
	return start - 1, end, true // 1-based inclusive → 0-based half-open
}

// Provides reports whether the GTF supplies a given annotation field, matching
// [AnnotationSource.Provides].
//
// The whole-file model knows because it has read everything. This one has read
// nothing until it is queried, which is the point of it — so for the two fields
// that vary between GTFs it samples the file instead: a bounded read of the
// first records it can find, checking whether they carry the attribute.
//
// A sample rather than a scan because the answer is a property of the file's
// attribute schema, not of any record. A GTF is emitted by one pipeline with one
// set of attributes, so gene_biotype is on every gene row or on none; reading a
// few hundred settles it, and reading the rest would cost exactly what the index
// was for.
//
// Unknown resolves to true. When the sample finds no records at all — an empty
// GTF, or one whose first contig has none — the honest answer is "cannot tell",
// and declaring a field that is never written costs a header line nobody reads
// while writing one never declared is a record a strict parser may reject.
func (s *IndexedAnnotationSource) Provides(key string) bool {
	switch key {
	case "gene_id", "gene_name", "start", "end", "strand":
		return true
	case "biotype":
		s.probe()
		return s.hasBioType
	case "status":
		s.probe()
		return s.hasStatus
	}
	return false
}

// probeLines bounds the sample, in lines read rather than gene rows found.
//
// Lines, because that is what the read costs. A budget spent on findings does
// not bound anything: a file whose rows mostly lack a gene_id would be read
// until enough gene rows turned up, which could be all of it.
//
// A thousand is a couple of bgzf blocks and several hundred gene rows in any
// ordinary GTF — far past the point where a file carrying biotypes would have
// shown one.
const probeLines = 1000

// probe samples the file once, recording which optional attributes it carries.
//
// The budget is spent across contigs, not per contig, so a GTF whose first
// contigs hold nothing parseable cannot make this read the whole file a
// thousand lines at a time.
func (s *IndexedAnnotationSource) probe() {
	if s.probed {
		return
	}
	s.probed = true
	// Unknown until proven otherwise, so a file this cannot read reports both
	// present — see Provides on why that is the safe direction.
	s.hasBioType, s.hasStatus = true, true

	lines, genes := 0, 0
	bio, status := false, false
	for _, ref := range s.tr.RefNames() {
		if lines >= probeLines {
			break
		}
		seq, err := s.tr.Query(ref, 0, math.MaxInt32)
		if err != nil {
			continue
		}
		for rec, qErr := range seq {
			if qErr != nil {
				return // cannot sample; leave the unknown answer in place
			}
			if lines++; lines > probeLines {
				break
			}
			cols := strings.Split(rec.Line, "\t")
			if len(cols) < 9 {
				continue
			}
			a := parseAttributes(cols[8])
			if a.geneID == "" {
				// Not a gene row as far as anything else here is concerned —
				// overlappingSpans skips these too. Counting one as a finding
				// would let a file of nothing but malformed rows report that it
				// carries no biotypes, which is a different claim from "cannot
				// tell".
				continue
			}
			genes++
			if a.bioType != "" {
				bio = true
			}
			if a.status != "" {
				status = true
			}
			if bio && status {
				break // nothing left to learn
			}
		}
		if bio && status {
			break
		}
	}
	if genes == 0 {
		return // nothing to go on; the unknown answer stands
	}
	s.hasBioType, s.hasStatus = bio, status
}
