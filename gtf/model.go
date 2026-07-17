package gtf

import (
	"sort"

	"github.com/compgenlab/cghts/bed"
)

// Exon is a single feature row of a transcript — an exon, a CDS segment, or a
// start/stop codon. Coordinates are 0-based half-open. Attributes holds the GTF
// attributes that are not recognized gene/transcript keys, in file order. Ports
// GTFAnnotationSource.GTFExon.
type Exon struct {
	Start      int
	End        int
	Attributes [][2]string
}

// Attribute returns the first value for key, or "" if the exon has no such
// attribute.
func (e *Exon) Attribute(key string) string {
	for _, kv := range e.Attributes {
		if kv[0] == key {
			return kv[1]
		}
	}
	return ""
}

func (e *Exon) region(ref string) span {
	return span{ref: ref, start: e.Start, end: e.End, strand: bed.StrandNone}
}

// Transcript is one transcript of a gene: its exons and, for coding
// transcripts, CDS segments and start/stop codons. CdsStart/CdsEnd bound the
// coding sequence; CdsEnd is the half-open end of the last CDS segment — i.e.
// one past the last coding base, and (per GTF convention, where the stop codon
// is a separate feature) it excludes the stop codon. Ports
// GTFAnnotationSource.GTFTranscript.
type Transcript struct {
	ID          string
	Start       int // -1 until the first exon is added
	End         int
	CdsStart    int // -1 until the first CDS is added
	CdsEnd      int
	Exons       []*Exon
	CDS         []*Exon
	StartCodons []*Exon
	StopCodons  []*Exon
}

func newTranscript(id string) *Transcript {
	return &Transcript{ID: id, Start: -1, End: -1, CdsStart: -1, CdsEnd: -1}
}

// HasCDS reports whether the transcript is coding.
func (t *Transcript) HasCDS() bool { return t.CdsStart > -1 && t.CdsEnd > -1 }

func (t *Transcript) addExon(start, end int, attrs [][2]string) {
	if t.Start == -1 || t.Start > start {
		t.Start = start
	}
	if t.End == -1 || t.End < end {
		t.End = end
	}
	t.Exons = append(t.Exons, &Exon{Start: start, End: end, Attributes: attrs})
	sortExons(t.Exons)
}

func (t *Transcript) addCDS(start, end int, attrs [][2]string) {
	if t.CdsStart == -1 || t.CdsStart > start {
		t.CdsStart = start
	}
	if t.CdsEnd == -1 || t.CdsEnd < end {
		t.CdsEnd = end
	}
	t.CDS = append(t.CDS, &Exon{Start: start, End: end, Attributes: attrs})
	sortExons(t.CDS)
}

func (t *Transcript) addStartCodon(start, end int, attrs [][2]string) {
	t.StartCodons = append(t.StartCodons, &Exon{Start: start, End: end, Attributes: attrs})
	sortExons(t.StartCodons)
}

func (t *Transcript) addStopCodon(start, end int, attrs [][2]string) {
	t.StopCodons = append(t.StopCodons, &Exon{Start: start, End: end, Attributes: attrs})
	sortExons(t.StopCodons)
}

func sortExons(xs []*Exon) {
	sort.Slice(xs, func(i, j int) bool {
		if xs[i].Start != xs[j].Start {
			return xs[i].Start < xs[j].Start
		}
		return xs[i].End < xs[j].End
	})
}

// Gene aggregates the transcripts that share a gene_id on one chromosome.
// Coordinates are 0-based half-open. Ports GTFAnnotationSource.GTFGene.
type Gene struct {
	Ref         string
	GeneID      string
	GeneName    string
	BioType     string
	Status      string
	Start       int
	End         int
	Strand      bed.Strand
	Transcripts map[string]*Transcript
}

func newGene(geneID, geneName, ref string, start, end int, strand bed.Strand, bioType, status string) *Gene {
	return &Gene{
		Ref: ref, GeneID: geneID, GeneName: geneName, BioType: bioType, Status: status,
		Start: start, End: end, Strand: strand,
		Transcripts: map[string]*Transcript{},
	}
}

func (g *Gene) transcript(id string) *Transcript {
	t, ok := g.Transcripts[id]
	if !ok {
		t = newTranscript(id)
		g.Transcripts[id] = t
	}
	return t
}

func (g *Gene) addExon(txID string, start, end int, attrs [][2]string) {
	g.transcript(txID).addExon(start, end, attrs)
	// Only exons extend the gene span (matching GTFGene.addExon); CDS/codon
	// features do not.
	if g.Start == -1 || g.Start > start {
		g.Start = start
	}
	if g.End == -1 || g.End < end {
		g.End = end
	}
}

func (g *Gene) addCDS(txID string, start, end int, attrs [][2]string) {
	g.transcript(txID).addCDS(start, end, attrs)
}

func (g *Gene) addStartCodon(txID string, start, end int, attrs [][2]string) {
	g.transcript(txID).addStartCodon(start, end, attrs)
}

func (g *Gene) addStopCodon(txID string, start, end int, attrs [][2]string) {
	g.transcript(txID).addStopCodon(start, end, attrs)
}

func (g *Gene) coord() span {
	return span{ref: g.Ref, start: g.Start, end: g.End, strand: g.Strand}
}

// IsCoding reports whether the gene has any coding transcript (one with a CDS).
// This is the structural definition used by the region classifier, independent
// of whether a biotype attribute is present.
func (g *Gene) IsCoding() bool {
	for _, t := range g.Transcripts {
		if t.HasCDS() {
			return true
		}
	}
	return false
}

// SortedTranscripts returns the gene's transcripts ordered by transcript ID.
// (The classifier iterates in this stable order; the Java original iterates an
// unordered HashMap, so this also makes region results deterministic.)
func (g *Gene) SortedTranscripts() []*Transcript {
	out := make([]*Transcript, 0, len(g.Transcripts))
	for _, t := range g.Transcripts {
		out = append(out, t)
	}
	sort.Slice(out, func(i, j int) bool { return out[i].ID < out[j].ID })
	return out
}

// Exons returns the gene's exons de-duplicated across transcripts (by
// coordinate) and sorted by position. Ports GTFGene.getExons().
func (g *Gene) Exons() []*Exon {
	seen := make(map[[2]int]*Exon)
	for _, t := range g.Transcripts {
		for _, e := range t.Exons {
			seen[[2]int{e.Start, e.End}] = e
		}
	}
	out := make([]*Exon, 0, len(seen))
	for _, e := range seen {
		out = append(out, e)
	}
	sortExons(out)
	return out
}
