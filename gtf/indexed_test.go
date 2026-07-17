package gtf

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/compgenlab/cghts/bed"
	"github.com/compgenlab/cghts/htsio/tabix"
)

// fixtureGTF is a small gene model: a coding gene A (5'UTR/coding/intron/3'UTR), a
// non-coding gene B, and a non-coding gene C overlapping A's 3' end. Rows are
// intentionally out of position order — the tabix writer sorts them.
var fixtureGTF = []string{
	`chr1	t	exon	1900	2200	.	+	.	gene_id "C"; transcript_id "C1";`,
	`chr1	t	gene	5000	6000	.	-	.	gene_id "B"; gene_name "GENEB"; gene_type "lincRNA";`,
	`chr1	t	transcript	5000	6000	.	-	.	gene_id "B"; transcript_id "B1";`,
	`chr1	t	exon	5000	5200	.	-	.	gene_id "B"; transcript_id "B1";`,
	`chr1	t	exon	5800	6000	.	-	.	gene_id "B"; transcript_id "B1";`,
	`chr1	t	gene	1000	2000	.	+	.	gene_id "A"; gene_name "GENEA"; gene_type "protein_coding";`,
	`chr1	t	transcript	1000	2000	.	+	.	gene_id "A"; transcript_id "A1";`,
	`chr1	t	exon	1000	1200	.	+	.	gene_id "A"; transcript_id "A1";`,
	`chr1	t	CDS	1100	1200	.	+	0	gene_id "A"; transcript_id "A1";`,
	`chr1	t	exon	1500	2000	.	+	.	gene_id "A"; transcript_id "A1";`,
	`chr1	t	CDS	1500	1800	.	+	1	gene_id "A"; transcript_id "A1";`,
	`chr1	t	gene	1900	2200	.	+	.	gene_id "C"; gene_name "GENEC"; gene_type "lincRNA";`,
	`chr1	t	transcript	1900	2200	.	+	.	gene_id "C"; transcript_id "C1";`,
}

// writeIndexedFixture writes fixtureGTF as a bgzipped, GFF-tabix-indexed file and
// returns its path (with a sibling .tbi).
func writeIndexedFixture(t *testing.T, dir string) string {
	t.Helper()
	path := filepath.Join(dir, "genes.gtf.gz")
	w := tabix.NewWriter(path, tabix.NewWriterOpts().GFF().AutoIndex())
	for _, line := range fixtureGTF {
		if err := w.Write(line); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	return path
}

// writePlainFixture writes fixtureGTF as a plain .gtf and returns its path.
func writePlainFixture(t *testing.T, dir string) string {
	t.Helper()
	path := filepath.Join(dir, "genes.gtf")
	var b []byte
	for _, line := range fixtureGTF {
		b = append(b, line...)
		b = append(b, '\n')
	}
	if err := os.WriteFile(path, b, 0o644); err != nil {
		t.Fatal(err)
	}
	return path
}

// TestIndexedMatchesInMemory: the tabix-backed reader returns the same overlapping
// genes and the same genic-region code as the full in-memory model, for a battery
// of positions (exon/intron/UTR/overlap/intergenic).
func TestIndexedMatchesInMemory(t *testing.T) {
	dir := t.TempDir()
	mem, err := NewAnnotationSource(writePlainFixture(t, dir), nil)
	if err != nil {
		t.Fatalf("in-memory: %v", err)
	}
	idx, err := NewIndexedAnnotationSource(writeIndexedFixture(t, dir), nil)
	if err != nil {
		t.Fatalf("indexed: %v", err)
	}
	defer idx.Close()

	// 0-based positions (file coords minus 1).
	positions := []int{1049, 1149, 1299, 1599, 1899, 1949, 5099, 5499, 8999}
	sawCode := map[string]bool{}
	for _, pos := range positions {
		memGenes := idsOf(mem.FindGenes("chr1", pos, pos+1))
		idxGenes := idsOf(idx.FindGenes("chr1", pos, pos+1))
		if !eqStrings(memGenes, idxGenes) {
			t.Errorf("pos %d: FindGenes indexed=%v want %v", pos, idxGenes, memGenes)
		}
		for _, gid := range memGenes {
			m := mem.FindGenicRegionForPos("chr1", pos, bed.StrandNone, gid).Code
			i := idx.FindGenicRegionForPos("chr1", pos, bed.StrandNone, gid).Code
			if m != i {
				t.Errorf("pos %d gene %s: region indexed=%q want %q", pos, gid, i, m)
			}
			sawCode[m] = true
		}
	}
	// Sanity: the fixture exercises real classifications, not a vacuous all-empty match.
	for _, want := range []string{"5_utr", "coding_exon", "coding_intron", "3_utr", "nc_exon", "nc_intron"} {
		if !sawCode[want] {
			t.Errorf("expected to observe region code %q across the fixture; got %v", want, keys(sawCode))
		}
	}
	// The overlap position must return both A and C.
	if g := idsOf(idx.FindGenes("chr1", 1949, 1950)); len(g) != 2 {
		t.Errorf("overlap pos: want 2 genes (A,C), got %v", g)
	}
}

// TestIndexedCacheBounded: the per-gene LRU never exceeds its cap even when many
// distinct genes are queried.
func TestIndexedCacheBounded(t *testing.T) {
	dir := t.TempDir()
	idx, err := NewIndexedAnnotationSource(writeIndexedFixture(t, dir), nil)
	if err != nil {
		t.Fatal(err)
	}
	defer idx.Close()
	idx.cap = 2 // force eviction with only 3 genes
	for i := 0; i < 50; i++ {
		idx.FindGenes("chr1", 1149, 1150) // A
		idx.FindGenes("chr1", 5099, 5100) // B
		idx.FindGenes("chr1", 1949, 1950) // A + C
	}
	if got := idx.ll.Len(); got > idx.cap {
		t.Errorf("cache holds %d entries, cap is %d", got, idx.cap)
	}
	if len(idx.cache) != idx.ll.Len() {
		t.Errorf("cache map (%d) and list (%d) out of sync", len(idx.cache), idx.ll.Len())
	}
}

func idsOf(genes []*Gene) []string {
	out := make([]string, len(genes))
	for i, g := range genes {
		out[i] = g.GeneID
	}
	return out
}

func eqStrings(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func keys(m map[string]bool) []string {
	out := make([]string, 0, len(m))
	for k := range m {
		out = append(out, k)
	}
	return out
}
