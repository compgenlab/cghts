package gtf

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/compgenlab/cghts/bed"
)

// testGTF is a small synthetic GTF (1-based inclusive coordinates) exercising:
//   - G1: a + strand coding gene with 5'UTR / CDS / 3'UTR across 4 exons and
//     introns of each flavor (tagged "basic")
//   - G2: a non-coding (lncRNA) gene, 2 exons
//   - G3: a - strand coding gene (to exercise the strand-flipped UTR logic)
//   - G4: a RefSeq-style gene using `gene` (not gene_name) and `gene_biotype`
//
// 0-based half-open spans (what the parser stores) are noted per line.
const testGTF = `#!genome-build test
chr1	test	gene	101	800	.	+	.	gene_id "G1"; gene_name "GeneOne"; gene_type "protein_coding";
chr1	test	exon	101	200	.	+	.	gene_id "G1"; gene_name "GeneOne"; transcript_id "T1"; gene_type "protein_coding"; tag "basic";
chr1	test	exon	301	400	.	+	.	gene_id "G1"; gene_name "GeneOne"; transcript_id "T1"; gene_type "protein_coding"; tag "basic";
chr1	test	exon	501	600	.	+	.	gene_id "G1"; gene_name "GeneOne"; transcript_id "T1"; gene_type "protein_coding"; tag "basic";
chr1	test	exon	701	800	.	+	.	gene_id "G1"; gene_name "GeneOne"; transcript_id "T1"; gene_type "protein_coding"; tag "basic";
chr1	test	CDS	351	400	.	+	0	gene_id "G1"; gene_name "GeneOne"; transcript_id "T1"; gene_type "protein_coding"; tag "basic";
chr1	test	CDS	501	550	.	+	0	gene_id "G1"; gene_name "GeneOne"; transcript_id "T1"; gene_type "protein_coding"; tag "basic";
chr1	test	exon	1001	1100	.	+	.	gene_id "G2"; gene_name "GeneTwo"; transcript_id "T2"; gene_type "lncRNA";
chr1	test	exon	1301	1400	.	+	.	gene_id "G2"; gene_name "GeneTwo"; transcript_id "T2"; gene_type "lncRNA";
chr1	test	exon	2001	2100	.	-	.	gene_id "G3"; gene_name "GeneThree"; transcript_id "T3"; gene_type "protein_coding";
chr1	test	exon	2301	2400	.	-	.	gene_id "G3"; gene_name "GeneThree"; transcript_id "T3"; gene_type "protein_coding";
chr1	test	CDS	2051	2100	.	-	0	gene_id "G3"; gene_name "GeneThree"; transcript_id "T3"; gene_type "protein_coding";
chr1	test	CDS	2301	2350	.	-	0	gene_id "G3"; gene_name "GeneThree"; transcript_id "T3"; gene_type "protein_coding";
chr2	test	exon	11	20	.	+	.	gene_id "G4"; gene "RefGene"; transcript_id "T4"; gene_biotype "protein_coding";
chr2	test	CDS	13	18	.	+	0	gene_id "G4"; gene "RefGene"; transcript_id "T4"; gene_biotype "protein_coding";
`

func writeGTF(t *testing.T, content string) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), "test.gtf")
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("write gtf: %v", err)
	}
	return path
}

func loadTestGTF(t *testing.T, requiredTags ...string) *AnnotationSource {
	t.Helper()
	src, err := NewAnnotationSource(writeGTF(t, testGTF), requiredTags)
	if err != nil {
		t.Fatalf("NewAnnotationSource: %v", err)
	}
	return src
}

func geneByID(s *AnnotationSource, id string) *Gene {
	for _, g := range s.Genes() {
		if g.GeneID == id {
			return g
		}
	}
	return nil
}

func TestParseStructure(t *testing.T) {
	s := loadTestGTF(t)
	if s.Size() != 4 {
		t.Fatalf("gene count = %d, want 4", s.Size())
	}
	if !s.Provides("biotype") {
		t.Errorf("Provides(biotype) = false, want true")
	}
	if s.Provides("status") {
		t.Errorf("Provides(status) = true, want false")
	}
	if got := s.RefNames(); len(got) != 2 || got[0] != "chr1" || got[1] != "chr2" {
		t.Errorf("RefNames = %v, want [chr1 chr2]", got)
	}

	g1 := geneByID(s, "G1")
	if g1 == nil {
		t.Fatal("G1 not found")
	}
	// Gene span comes from the exons: [100, 800).
	if g1.Start != 100 || g1.End != 800 {
		t.Errorf("G1 span = [%d,%d), want [100,800)", g1.Start, g1.End)
	}
	if g1.Strand != bed.StrandPlus {
		t.Errorf("G1 strand = %v, want +", g1.Strand)
	}
	if g1.BioType != "protein_coding" {
		t.Errorf("G1 biotype = %q, want protein_coding", g1.BioType)
	}
	if len(g1.Transcripts) != 1 {
		t.Fatalf("G1 transcripts = %d, want 1", len(g1.Transcripts))
	}
	t1 := g1.Transcripts["T1"]
	if len(t1.Exons) != 4 {
		t.Errorf("T1 exons = %d, want 4", len(t1.Exons))
	}
	if !t1.HasCDS() {
		t.Error("T1 should be coding")
	}
	// CDS [350,400) + [500,550): cdsStart=350, cdsEnd=550.
	if t1.CdsStart != 350 || t1.CdsEnd != 550 {
		t.Errorf("T1 cds = [%d,%d), want [350,550)", t1.CdsStart, t1.CdsEnd)
	}
}

func TestParseRefSeqGeneNameFallback(t *testing.T) {
	s := loadTestGTF(t)
	g4 := geneByID(s, "G4")
	if g4 == nil {
		t.Fatal("G4 not found")
	}
	if g4.GeneName != "RefGene" {
		t.Errorf("G4 gene_name = %q, want RefGene (from `gene`)", g4.GeneName)
	}
	if g4.BioType != "protein_coding" {
		t.Errorf("G4 biotype = %q, want protein_coding (from gene_biotype)", g4.BioType)
	}
}

// TestParseRefSeqBiotypeBackfill reproduces the RefSeq layout after coordinate
// sorting: the gene_biotype lives only on the "gene" feature line, and a
// transcript/exon row for the same gene_id sorts ahead of it. The gene must still
// pick up its biotype from the later "gene" row rather than being frozen with the
// empty biotype of whichever row seeded it first.
func TestParseRefSeqBiotypeBackfill(t *testing.T) {
	const refseq = `chr3	BestRefSeq	transcript	101	800	.	+	.	gene_id "G5"; transcript_id "NM_1"; gene "RefFive"; transcript_biotype "mRNA";
chr3	BestRefSeq	exon	101	200	.	+	.	gene_id "G5"; transcript_id "NM_1"; gene "RefFive"; transcript_biotype "mRNA";
chr3	BestRefSeq	CDS	101	200	.	+	0	gene_id "G5"; transcript_id "NM_1"; gene "RefFive"; transcript_biotype "mRNA";
chr3	BestRefSeq,Gnomon	gene	101	800	.	+	.	gene_id "G5"; gene "RefFive"; gene_biotype "protein_coding";
`
	s, err := NewAnnotationSource(writeGTF(t, refseq), nil)
	if err != nil {
		t.Fatalf("NewAnnotationSource: %v", err)
	}
	g5 := geneByID(s, "G5")
	if g5 == nil {
		t.Fatal("G5 not found")
	}
	if g5.BioType != "protein_coding" {
		t.Errorf("G5 biotype = %q, want protein_coding (backfilled from the gene row)", g5.BioType)
	}
	if !s.hasBioType {
		t.Error("source hasBioType = false, want true")
	}
}

func TestRequiredTags(t *testing.T) {
	// Only G1's rows carry tag "basic"; with the filter, only G1 survives.
	s := loadTestGTF(t, "basic")
	if s.Size() != 1 {
		t.Fatalf("gene count with tag filter = %d, want 1", s.Size())
	}
	if geneByID(s, "G1") == nil {
		t.Error("G1 should survive the tag filter")
	}
}

func TestFindGenes(t *testing.T) {
	s := loadTestGTF(t)
	// A position inside G1 exon1.
	genes := s.FindGenes("chr1", 150, 151)
	if len(genes) != 1 || genes[0].GeneID != "G1" {
		t.Errorf("FindGenes@150 = %v, want [G1]", geneIDs(genes))
	}
	// Intergenic gap.
	if genes := s.FindGenes("chr1", 5000, 5001); len(genes) != 0 {
		t.Errorf("FindGenes@5000 = %v, want []", geneIDs(genes))
	}
	// Unknown contig.
	if s.HasRef("chrZ") {
		t.Error("HasRef(chrZ) = true, want false")
	}
}

func geneIDs(genes []*Gene) []string {
	out := make([]string, len(genes))
	for i, g := range genes {
		out[i] = g.GeneID
	}
	return out
}
