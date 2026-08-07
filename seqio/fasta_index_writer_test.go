package seqio

import (
	"bytes"
	"compress/gzip"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/htsio/bgzf"
)

const idxFasta = ">chr1 with a description\n" +
	"ACGTACGTACGTACGTACGT\nACGTACGTACGTACGTACGT\nACGTACGT\n" +
	">chr2\n" +
	"TTTTAAAACCCCGGGG\nTTTTAAAACCCCGGGG\nTTTT\n" +
	">chr3\n" +
	"NNNNNNNNNN\n"

// bigFasta spans several BGZF blocks, so a .gzi carries real entries. A small
// fixture only ever proves the count is zero.
func bigFasta() string {
	var b strings.Builder
	line := strings.Repeat("ACGT", 15) // 60 bases
	for c := 1; c <= 3; c++ {
		fmt.Fprintf(&b, ">chr%d some description\n", c)
		for i := 0; i < 1500; i++ {
			b.WriteString(line)
			b.WriteByte('\n')
		}
		b.WriteString("ACGTACGT\n") // short final line
	}
	return b.String()
}

// bgzipTo compresses body into path using this package's own BGZF writer.
func bgzipTo(t *testing.T, path, body string) {
	t.Helper()
	w, err := bgzf.NewBGZipFile(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := w.Write([]byte(body)); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
}

// samtools is the oracle. Our indexes are only correct if htslib's readers
// accept them, and the cheapest proof is producing byte-identical output to
// htslib's own writer.
func samtools(t *testing.T) string {
	t.Helper()
	p, err := exec.LookPath("samtools")
	if err != nil {
		t.Skip("samtools not installed; skipping the htslib byte comparison")
	}
	return p
}

// The index has to be usable, which means the reader in this package must read
// correct bases through it — including across a line boundary, where a wrong
// LineWidth shows up.
func TestBuildFastaIndexRoundTrips(t *testing.T) {
	dir := t.TempDir()
	p := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(p, []byte(idxFasta), 0o644); err != nil {
		t.Fatal(err)
	}
	entries, err := BuildFastaIndex(p)
	if err != nil {
		t.Fatal(err)
	}
	if len(entries) != 3 {
		t.Fatalf("indexed %d sequences, want 3", len(entries))
	}
	if entries[0].Name != "chr1" || entries[0].Length != 48 {
		t.Errorf("chr1 = %+v, want name chr1 length 48", entries[0])
	}

	r, err := NewIndexedFastaReader(p)
	if err != nil {
		t.Fatalf("our own reader rejected the index: %v", err)
	}
	defer r.Close()
	got, err := r.GetSequenceRange("chr2", 14, 18)
	if err != nil {
		t.Fatal(err)
	}
	if strings.ToUpper(string(got)) != "GGTT" {
		t.Errorf("chr2:14-18 = %q, want GGTT (spans a line boundary)", got)
	}
}

func TestFaiMatchesHtslib(t *testing.T) {
	sam := samtools(t)
	// Separate directories: both writers put the index beside the file.
	ours, theirs := t.TempDir(), t.TempDir()
	for _, d := range []string{ours, theirs} {
		if err := os.WriteFile(filepath.Join(d, "ref.fa"), []byte(idxFasta), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	if out, err := exec.Command(sam, "faidx", filepath.Join(theirs, "ref.fa")).CombinedOutput(); err != nil {
		t.Fatalf("samtools faidx: %v\n%s", err, out)
	}
	if _, err := BuildFastaIndex(filepath.Join(ours, "ref.fa")); err != nil {
		t.Fatal(err)
	}
	want, _ := os.ReadFile(filepath.Join(theirs, "ref.fa.fai"))
	got, _ := os.ReadFile(filepath.Join(ours, "ref.fa.fai"))
	if string(got) != string(want) {
		t.Errorf(".fai differs from htslib.\n--- ours ---\n%s--- samtools ---\n%s", got, want)
	}
}

// Both index the same compressed bytes, so the block boundaries are fixed and
// the .gzi must agree exactly — this catches an off-by-one in the block walk
// that a "does it load" check would not. The EOF marker is the specific trap:
// counting it yields one entry too many.
func TestFaiAndGziMatchHtslibForBGZF(t *testing.T) {
	sam := samtools(t)

	src := t.TempDir()
	gz := filepath.Join(src, "ref.fa.gz")
	bgzipTo(t, gz, bigFasta())
	body, err := os.ReadFile(gz)
	if err != nil {
		t.Fatal(err)
	}

	ours, theirs := t.TempDir(), t.TempDir()
	for _, d := range []string{ours, theirs} {
		if err := os.WriteFile(filepath.Join(d, "ref.fa.gz"), body, 0o644); err != nil {
			t.Fatal(err)
		}
	}
	if o, err := exec.Command(sam, "faidx", filepath.Join(theirs, "ref.fa.gz")).CombinedOutput(); err != nil {
		t.Fatalf("samtools faidx (bgzf): %v\n%s", err, o)
	}
	if _, err := BuildFastaIndex(filepath.Join(ours, "ref.fa.gz")); err != nil {
		t.Fatal(err)
	}
	for _, ext := range []string{".fai", ".gzi"} {
		want, err := os.ReadFile(filepath.Join(theirs, "ref.fa.gz"+ext))
		if err != nil {
			t.Fatalf("samtools wrote no %s: %v", ext, err)
		}
		got, err := os.ReadFile(filepath.Join(ours, "ref.fa.gz"+ext))
		if err != nil {
			t.Fatalf("we wrote no %s: %v", ext, err)
		}
		if !bytes.Equal(got, want) {
			t.Errorf("%s differs from htslib (ours %d bytes, samtools %d)", ext, len(got), len(want))
		}
	}
}

// A .fai addresses uncompressed offsets, so it is only usable with a container
// that can seek to one. Plain gzip cannot, and an unusable index is worse than
// a refusal.
func TestBuildFastaIndexRefusesPlainGzip(t *testing.T) {
	dir := t.TempDir()
	p := filepath.Join(dir, "ref.fa.gz")
	f, err := os.Create(p)
	if err != nil {
		t.Fatal(err)
	}
	zw := gzip.NewWriter(f) // no BGZF extra field
	if _, err := zw.Write([]byte(idxFasta)); err != nil {
		t.Fatal(err)
	}
	zw.Close()
	f.Close()

	if _, err := BuildFastaIndex(p); err == nil {
		t.Fatal("a plain-gzip FASTA was accepted; its index would be unusable")
	} else if !strings.Contains(err.Error(), "BGZF") {
		t.Errorf("error does not name the problem: %v", err)
	}
}

// A ragged record cannot be described by a .fai, and producing one anyway would
// give wrong coordinates rather than an error.
func TestBuildFastaIndexRejectsRaggedLines(t *testing.T) {
	dir := t.TempDir()
	p := filepath.Join(dir, "ragged.fa")
	if err := os.WriteFile(p, []byte(">chr1\nACGTACGT\nACG\nACGTACGT\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if _, err := BuildFastaIndex(p); err == nil {
		t.Error("uneven line lengths were accepted")
	}
}
