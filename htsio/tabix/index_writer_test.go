package tabix

import (
	"io"
	"path/filepath"
	"testing"

	"github.com/compgenlab/cghts/htsio/bgzf"
)

// writeBGZF writes the given lines to a BGZF file with no index.
func writeBGZF(t *testing.T, path string, lines ...string) {
	t.Helper()
	w, err := bgzf.NewBGZipFile(path)
	if err != nil {
		t.Fatal(err)
	}
	for _, l := range lines {
		if _, err := io.WriteString(w, l+"\n"); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
}

func TestIndexWriterBED(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "data.bed.gz")
	writeBGZF(t, path,
		"chr1\t90\t110\tgeneA",
		"chr1\t145\t155\tenhB",
		"chr2\t400\t600\tgeneC") // single chr2 record (and the last line)

	if err := NewIndexWriter(NewWriterOpts().BED()).WriteIndex(path); err != nil {
		t.Fatalf("WriteIndex: %v", err)
	}

	r, err := NewReader(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()

	cases := []struct {
		ref        string
		start, end int
		want       string
	}{
		{"chr1", 99, 100, "chr1\t90\t110\tgeneA"},
		{"chr1", 149, 150, "chr1\t145\t155\tenhB"},
		{"chr2", 499, 500, "chr2\t400\t600\tgeneC"},
	}
	for _, c := range cases {
		got := queryLines(t, r, c.ref, c.start, c.end)
		if len(got) != 1 || got[0] != c.want {
			t.Errorf("query %s:%d-%d = %v, want [%q]", c.ref, c.start, c.end, got, c.want)
		}
	}
}

func TestIndexWriterVCFWithHeader(t *testing.T) {
	// A '#'-prefixed header line must be skipped (meta='#'); 1-based positions.
	dir := t.TempDir()
	path := filepath.Join(dir, "v.vcf.gz")
	writeBGZF(t, path,
		"##fileformat=VCFv4.2",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
		"chr1\t100\t.\tA\tG\t.\tPASS\t.",
		"chr2\t500\t.\tC\tT\t.\tPASS\t.")

	if err := NewIndexWriter(NewWriterOpts().VCF()).WriteIndex(path); err != nil {
		t.Fatalf("WriteIndex: %v", err)
	}
	r, err := NewReader(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	if got := queryLines(t, r, "chr2", 499, 500); len(got) != 1 {
		t.Errorf("chr2 query = %v, want 1 row", got)
	}
	if got := queryLines(t, r, "chr1", 99, 100); len(got) != 1 {
		t.Errorf("chr1 query = %v, want 1 row", got)
	}
}
