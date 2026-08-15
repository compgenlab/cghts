package annotate

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/gtf"
	"github.com/compgenlab/cghts/htsio/tabix"
)

// writeIndexedGTF writes the same fixture bgzipped with a tabix index beside it.
func writeIndexedGTF(t *testing.T) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), "ann.gtf.gz")
	w := tabix.NewWriter(path, tabix.NewWriterOpts().Columns(1, 4, 5).Meta('#').AutoIndex())
	for _, line := range strings.Split(strings.TrimRight(annGTF, "\n"), "\n") {
		if line == "" {
			continue
		}
		if err := w.Write(line); err != nil {
			t.Fatalf("write %q: %v", line, err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatalf("close: %v", err)
	}
	return path
}

// An indexed GTF is queried through its index rather than read into memory.
//
// This is the whole reason the annotator holds an interface. Reading a GTF
// wholly in is fine for a small annotation set and is not fine for a human
// GENCODE — a caller annotating a genome against two of them ran out of memory,
// and worked around it by writing a second annotator rather than fixing this
// one. The workaround is the thing being removed.
func TestAnIndexedGtfIsQueriedThroughItsIndex(t *testing.T) {
	path := writeIndexedGTF(t)
	if _, err := os.Stat(path + ".tbi"); err != nil {
		t.Fatalf("the fixture has no index, so this proves nothing: %v", err)
	}

	a, err := NewGtfAnnotator(GtfOptions{Filename: path})
	if err != nil {
		t.Fatalf("NewGtfAnnotator: %v", err)
	}
	defer a.Close()

	// The index was chosen, not the whole-file reader.
	if _, ok := a.src.(*gtf.IndexedAnnotationSource); !ok {
		t.Fatalf("gene model is %T; an indexed GTF should be queried through its "+
			"index, which is what bounds memory to the genes actually looked at", a.src)
	}
}

// An unindexed GTF still works, read wholly into memory.
//
// Falling back rather than refusing: an unindexed GTF is a working GTF, it just
// costs memory proportional to the file instead of to the query. Refusing one
// would break every caller that has been passing plain files.
func TestAnUnindexedGtfStillWorks(t *testing.T) {
	a, err := NewGtfAnnotator(GtfOptions{Filename: writeAnnGTF(t)})
	if err != nil {
		t.Fatalf("NewGtfAnnotator: %v", err)
	}
	defer a.Close()
	if _, ok := a.src.(*gtf.AnnotationSource); !ok {
		t.Errorf("gene model is %T; an unindexed GTF has to be read in", a.src)
	}
}

// Both backends annotate identically.
//
// The point of the interface is that the choice is invisible above it. A
// difference here would mean the memory saving had quietly changed the answer,
// which is the failure the second implementation was risking all along.
func TestBothGeneModelsAnnotateTheSame(t *testing.T) {
	plainA, err := NewGtfAnnotator(GtfOptions{Filename: writeAnnGTF(t)})
	if err != nil {
		t.Fatal(err)
	}
	defer plainA.Close()
	indexedA, err := NewGtfAnnotator(GtfOptions{Filename: writeIndexedGTF(t)})
	if err != nil {
		t.Fatal(err)
	}
	defer indexedA.Close()

	const line = "chr1\t170\t.\tA\tG\t.\tPASS\t.\tGT\t0/1" // GeneA + GeneB overlap
	for _, key := range []string{"GTF_GENE", "GTF_GENEID", "GTF_STRAND", "GTF_REGION"} {
		var got [2]string
		for i, a := range []*GtfAnnotator{plainA, indexedA} {
			h, recs := bedRecs(t, line)
			if err := a.SetupHeader(h); err != nil {
				t.Fatal(err)
			}
			if err := a.Annotate(recs[0]); err != nil {
				t.Fatal(err)
			}
			got[i], _ = info(t, recs[0], key)
		}
		if got[0] != got[1] {
			t.Errorf("%s: in-memory %q, indexed %q — the two gene models disagree",
				key, got[0], got[1])
		}
		if got[0] == "" {
			t.Errorf("%s came back empty from both; the comparison is vacuous", key)
		}
	}
}
