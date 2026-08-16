package annotate

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/gtf"
	"github.com/compgenlab/cghts/htsio/tabix"
)

// writeIndexedGTF writes the shared fixture bgzipped with a tabix index beside it.
func writeIndexedGTF(t *testing.T) string {
	t.Helper()
	return writeIndexedFrom(t, annGTF)
}

// writeIndexedFrom indexes arbitrary GTF text, so a test can supply its own.
func writeIndexedFrom(t *testing.T, body string) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), "ann.gtf.gz")
	w := tabix.NewWriter(path, tabix.NewWriterOpts().Columns(1, 4, 5).Meta('#').AutoIndex())
	for _, line := range strings.Split(strings.TrimRight(body, "\n"), "\n") {
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

// A caller gets the fields it selected and no others.
//
// The seven are one annotator's output but not one annotation. A configuration
// naming two of them should produce two INFO fields, not seven with five nobody
// asked for on every record of a whole-genome VCF — which is the other reason a
// caller went and wrote its own GTF annotator instead of using this one.
func TestOnlyTheSelectedGtfFieldsAreWritten(t *testing.T) {
	a, err := NewGtfAnnotator(GtfOptions{
		Filename: writeAnnGTF(t),
		Fields:   []string{GtfGeneSymbol, GtfRegion},
	})
	if err != nil {
		t.Fatal(err)
	}
	defer a.Close()

	h, recs := bedRecs(t, "chr1\t170\t.\tA\tG\t.\tPASS\t.\tGT\t0/1")
	if err := a.SetupHeader(h); err != nil {
		t.Fatal(err)
	}
	if err := a.Annotate(recs[0]); err != nil {
		t.Fatal(err)
	}

	for _, want := range []string{"GTF_GENE", "GTF_REGION"} {
		if _, ok := info(t, recs[0], want); !ok {
			t.Errorf("%s was selected but not written", want)
		}
		if _, ok := h.InfoDef(want); !ok {
			t.Errorf("%s was selected but not declared", want)
		}
	}
	for _, unwanted := range []string{"GTF_GENEID", "GTF_STRAND", "GTF_BIOTYPE", "GTF_CODING", "GTF_NONCODING"} {
		if _, ok := info(t, recs[0], unwanted); ok {
			t.Errorf("%s was written though it was not selected", unwanted)
		}
		if _, ok := h.InfoDef(unwanted); ok {
			t.Errorf("%s was declared though it was not selected", unwanted)
		}
	}
}

// Selecting nothing means all of them, which is what "the GTF annotations"
// means to a caller that has not asked for a subset.
func TestNoSelectionMeansEveryGtfField(t *testing.T) {
	a, err := NewGtfAnnotator(GtfOptions{Filename: writeAnnGTF(t)})
	if err != nil {
		t.Fatal(err)
	}
	defer a.Close()
	h, recs := bedRecs(t, "chr1\t170\t.\tA\tG\t.\tPASS\t.\tGT\t0/1")
	if err := a.SetupHeader(h); err != nil {
		t.Fatal(err)
	}
	if err := a.Annotate(recs[0]); err != nil {
		t.Fatal(err)
	}
	for _, want := range []string{"GTF_GENE", "GTF_GENEID", "GTF_STRAND", "GTF_REGION"} {
		if _, ok := info(t, recs[0], want); !ok {
			t.Errorf("%s missing when nothing was selected", want)
		}
	}
}

// The biotype field is declared when the GTF carries one and not when it does
// not — the same answer from either gene model.
//
// Until the indexed model could answer, the annotator assumed yes, so every
// indexed GTF got a header line for a field that might never be written. Both
// models sample or know now, so both headers describe the file in front of them.
func TestTheBiotypeFieldFollowsTheFile(t *testing.T) {
	// The shared fixture carries gene_type; this one carries nothing optional.
	const bare = "chr1\tt\texon\t101\t200\t.\t+\t.\t" +
		`gene_id "GeneA"; gene_name "GeneA"; transcript_id "TA";` + "\n" +
		"chr1\tt\tCDS\t101\t200\t.\t+\t0\t" +
		`gene_id "GeneA"; gene_name "GeneA"; transcript_id "TA";` + "\n"
	bareDir := t.TempDir()
	barePlain := filepath.Join(bareDir, "bare.gtf")
	if err := os.WriteFile(barePlain, []byte(bare), 0o644); err != nil {
		t.Fatal(err)
	}
	bareIndexed := writeIndexedFrom(t, bare)

	for _, tc := range []struct {
		name string
		path string
		want bool
	}{
		{"in-memory, carries biotypes", writeAnnGTF(t), true},
		{"indexed, carries biotypes", writeIndexedGTF(t), true},
		{"in-memory, carries none", barePlain, false},
		{"indexed, carries none", bareIndexed, false},
	} {
		a, err := NewGtfAnnotator(GtfOptions{Filename: tc.path})
		if err != nil {
			t.Fatalf("%s: %v", tc.name, err)
		}
		h, _ := bedRecs(t, "chr1\t170\t.\tA\tG\t.\tPASS\t.\tGT\t0/1")
		if err := a.SetupHeader(h); err != nil {
			t.Fatal(err)
		}
		a.Close()
		if _, ok := h.InfoDef("GTF_BIOTYPE"); ok != tc.want {
			t.Errorf("%s: biotype declared = %v, want %v", tc.name, ok, tc.want)
		}
		// Not vacuous: the fields it does carry are declared either way.
		if _, ok := h.InfoDef("GTF_GENE"); !ok {
			t.Errorf("%s: the header is empty, so this proves nothing", tc.name)
		}
	}
}
