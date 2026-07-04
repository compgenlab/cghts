package vcf

import (
	"io"
	"os"
	"path/filepath"
	"testing"

	"github.com/compgenlab/hts/htsio/bgzf"
)

// TestOpenVcfWriterBGZF: a ".gz" (or ".bgz") output is BGZF (block-gzip), not
// plain gzip — so it round-trips through the BGZF reader and is tabix-indexable.
// BGZF sets the gzip FEXTRA flag (0x04) with a "BC" subfield; stdlib gzip does
// not, so that flag distinguishes a regression back to plain gzip.
func TestOpenVcfWriterBGZF(t *testing.T) {
	for _, ext := range []string{".vcf.gz", ".vcf.bgz"} {
		path := filepath.Join(t.TempDir(), "out"+ext)
		w, err := OpenVcfWriter(path)
		if err != nil {
			t.Fatalf("OpenVcfWriter(%s): %v", path, err)
		}
		if err := w.WriteLine("##fileformat=VCFv4.2"); err != nil {
			t.Fatal(err)
		}
		if err := w.WriteLine("chr1\t100\t.\tA\tG\t.\t.\t."); err != nil {
			t.Fatal(err)
		}
		if err := w.Close(); err != nil {
			t.Fatalf("close: %v", err)
		}

		raw, err := os.ReadFile(path)
		if err != nil {
			t.Fatal(err)
		}
		// gzip magic + FEXTRA flag set (BGZF), which stdlib gzip.NewWriter omits.
		if len(raw) < 4 || raw[0] != 0x1f || raw[1] != 0x8b {
			t.Fatalf("%s: not a gzip stream", ext)
		}
		if raw[3]&0x04 == 0 {
			t.Errorf("%s: FEXTRA flag not set — output is plain gzip, not BGZF", ext)
		}

		// And it decodes back through the BGZF reader.
		f, err := os.Open(path)
		if err != nil {
			t.Fatal(err)
		}
		got, err := io.ReadAll(bgzf.NewReader(f))
		f.Close()
		if err != nil {
			t.Fatalf("%s: BGZF read: %v", ext, err)
		}
		want := "##fileformat=VCFv4.2\nchr1\t100\t.\tA\tG\t.\t.\t.\n"
		if string(got) != want {
			t.Errorf("%s: round-trip = %q, want %q", ext, got, want)
		}
	}
}
