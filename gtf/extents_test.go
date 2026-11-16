package gtf

import (
	"errors"
	"os"
	"path/filepath"
	"testing"

	"github.com/compgenlab/cghts/bed"
	"github.com/compgenlab/cghts/htsio/tabix"
)

func writeIndexed(t *testing.T, lines ...string) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), "g.gtf.gz")
	w := tabix.NewWriter(path, tabix.NewWriterOpts().GFF().AutoIndex())
	for _, l := range lines {
		if err := w.Write(l); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	return path
}

// exon and CDS rows only, as UCSC genePredToGtf writes: a two-exon coding
// transcript with an intron at [200,500).
var noExtentRows = []string{
	"chr1\tt\texon\t101\t200\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\";",
	"chr1\tt\tCDS\t121\t200\t.\t+\t0\tgene_id \"G\"; transcript_id \"T\";",
	"chr1\tt\texon\t501\t600\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\";",
	"chr1\tt\tCDS\t501\t550\t.\t+\t0\tgene_id \"G\"; transcript_id \"T\";",
}

// Through the index, a position in the intron overlaps none of G's rows, so it
// would read as intergenic -- the whole-file model, which sees every row, says
// coding intron. Rather than answer that wrongly, the indexed reader refuses
// such a GTF.
func TestIndexedRefusesAGtfWithoutGeneExtents(t *testing.T) {
	_, err := NewIndexedAnnotationSource(writeIndexed(t, noExtentRows...), nil)
	if !errors.Is(err, ErrNoGeneExtents) {
		t.Fatalf("err = %v, want ErrNoGeneExtents", err)
	}

	// The same rows read whole are fine, and show what the index would miss.
	path := filepath.Join(t.TempDir(), "g.gtf")
	body := ""
	for _, l := range noExtentRows {
		body += l + "\n"
	}
	if err := os.WriteFile(path, []byte(body), 0o644); err != nil {
		t.Fatal(err)
	}
	mem, err := NewAnnotationSource(path, nil)
	if err != nil {
		t.Fatal(err)
	}
	if got := mem.FindGenicRegionForPos("chr1", 300, bed.StrandNone, ""); got != CodingIntron {
		t.Errorf("whole-file model at the intron: %s, want coding_intron", got.Code)
	}
}

// Either a gene row or transcript rows give a gene its extent.
func TestIndexedAcceptsGeneOrTranscriptRows(t *testing.T) {
	for name, extra := range map[string]string{
		"gene row":       "chr1\tt\tgene\t101\t600\t.\t+\t.\tgene_id \"G\";",
		"transcript row": "chr1\tt\ttranscript\t101\t600\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\";",
	} {
		src, err := NewIndexedAnnotationSource(writeIndexed(t, append([]string{extra}, noExtentRows...)...), nil)
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		if got := src.FindGenicRegionForPos("chr1", 300, bed.StrandNone, ""); got != CodingIntron {
			t.Errorf("%s: intron reads as %s, want coding_intron", name, got.Code)
		}
		src.Close()
	}
}

// A sample with nothing to judge by is "cannot tell", not "has none".
func TestIndexedDoesNotRefuseAnEmptyGtf(t *testing.T) {
	src, err := NewIndexedAnnotationSource(writeIndexed(t), nil)
	if err != nil {
		t.Fatalf("an empty GTF should open: %v", err)
	}
	src.Close()
}
