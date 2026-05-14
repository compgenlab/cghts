package seqio

import (
	"bytes"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

// createTestFasta creates a plain FASTA + .fai in the given directory.
// Returns the FASTA file path.
func createTestFasta(t *testing.T, dir string) string {
	t.Helper()

	// Create a FASTA with two sequences. Use 80-char lines.
	var buf bytes.Buffer
	buf.WriteString(">seq1\n")
	seq1 := strings.Repeat("ACGTACGTNN", 200) // 2000 bases
	for i := 0; i < len(seq1); i += 80 {
		end := i + 80
		if end > len(seq1) {
			end = len(seq1)
		}
		buf.WriteString(seq1[i:end])
		buf.WriteByte('\n')
	}

	buf.WriteString(">seq2\n")
	seq2 := strings.Repeat("TTGGCCAA", 125) // 1000 bases
	for i := 0; i < len(seq2); i += 80 {
		end := i + 80
		if end > len(seq2) {
			end = len(seq2)
		}
		buf.WriteString(seq2[i:end])
		buf.WriteByte('\n')
	}

	fastaPath := filepath.Join(dir, "test.fa")
	if err := os.WriteFile(fastaPath, buf.Bytes(), 0644); err != nil {
		t.Fatal(err)
	}

	// Compute .fai entries.
	// seq1: 2000 bases, offset = len(">seq1\n") = 6, lineBases=80, lineWidth=81
	// seq2: 1000 bases, offset = 6 + 25*81 + len(">seq2\n") = 6 + 2025 + 6 = 2037
	seq1Offset := 6
	seq1Lines := (2000 + 79) / 80 // 25 lines
	seq2Offset := seq1Offset + seq1Lines*81 + 6

	faiContent := strings.Join([]string{
		"seq1\t2000\t" + itoa(seq1Offset) + "\t80\t81",
		"seq2\t1000\t" + itoa(seq2Offset) + "\t80\t81",
	}, "\n") + "\n"

	faiPath := fastaPath + ".fai"
	if err := os.WriteFile(faiPath, []byte(faiContent), 0644); err != nil {
		t.Fatal(err)
	}

	return fastaPath
}

func itoa(n int) string {
	return strconv.Itoa(n)
}

func TestIndexedFastaReaderGetSequence(t *testing.T) {
	dir := t.TempDir()
	fastaPath := createTestFasta(t, dir)

	r, err := NewIndexedFastaReader(fastaPath)
	if err != nil {
		t.Fatalf("NewIndexedFastaReader: %v", err)
	}
	defer r.Close()

	// Verify names.
	names := r.Names()
	if len(names) != 2 || names[0] != "seq1" || names[1] != "seq2" {
		t.Errorf("Names() = %v, want [seq1, seq2]", names)
	}

	// Verify lengths.
	if l, ok := r.SequenceLength("seq1"); !ok || l != 2000 {
		t.Errorf("SequenceLength(seq1) = %d, %v, want 2000, true", l, ok)
	}
	if l, ok := r.SequenceLength("seq2"); !ok || l != 1000 {
		t.Errorf("SequenceLength(seq2) = %d, %v, want 1000, true", l, ok)
	}

	// Get full seq1.
	seq1, err := r.GetSequence("seq1")
	if err != nil {
		t.Fatalf("GetSequence(seq1): %v", err)
	}
	expected1 := strings.Repeat("ACGTACGTNN", 200)
	if string(seq1) != expected1 {
		t.Errorf("GetSequence(seq1): got %d bytes, want %d", len(seq1), len(expected1))
		if len(seq1) > 0 && len(seq1) < 50 {
			t.Errorf("  got: %s", string(seq1[:min(50, len(seq1))]))
		}
	}

	// Get full seq2.
	seq2, err := r.GetSequence("seq2")
	if err != nil {
		t.Fatalf("GetSequence(seq2): %v", err)
	}
	expected2 := strings.Repeat("TTGGCCAA", 125)
	if string(seq2) != expected2 {
		t.Errorf("GetSequence(seq2): got %d bytes, want %d", len(seq2), len(expected2))
	}
}

func TestIndexedFastaReaderGetSequenceRange(t *testing.T) {
	dir := t.TempDir()
	fastaPath := createTestFasta(t, dir)

	r, err := NewIndexedFastaReader(fastaPath)
	if err != nil {
		t.Fatalf("NewIndexedFastaReader: %v", err)
	}
	defer r.Close()

	expected1 := strings.Repeat("ACGTACGTNN", 200)

	tests := []struct {
		name       string
		seq        string
		start, end int
		want       string
	}{
		{"first 10", "seq1", 0, 10, expected1[:10]},
		{"mid range", "seq1", 100, 200, expected1[100:200]},
		{"cross line boundary", "seq1", 75, 85, expected1[75:85]},
		{"last 10", "seq1", 1990, 2000, expected1[1990:2000]},
		{"full seq2", "seq2", 0, 1000, strings.Repeat("TTGGCCAA", 125)},
		{"clamped end", "seq1", 1990, 9999, expected1[1990:]},
		{"empty range", "seq1", 100, 100, ""},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := r.GetSequenceRange(tt.seq, tt.start, tt.end)
			if err != nil {
				t.Fatalf("GetSequenceRange(%s, %d, %d): %v", tt.seq, tt.start, tt.end, err)
			}
			if tt.want == "" {
				if len(got) != 0 {
					t.Errorf("expected empty, got %d bytes", len(got))
				}
				return
			}
			if string(got) != tt.want {
				t.Errorf("got %d bytes, want %d", len(got), len(tt.want))
				if len(got) < 50 {
					t.Errorf("  got:  %q", string(got))
					t.Errorf("  want: %q", tt.want[:min(50, len(tt.want))])
				}
			}
		})
	}
}

func TestIndexedFastaReaderErrors(t *testing.T) {
	dir := t.TempDir()
	fastaPath := createTestFasta(t, dir)

	r, err := NewIndexedFastaReader(fastaPath)
	if err != nil {
		t.Fatalf("NewIndexedFastaReader: %v", err)
	}
	defer r.Close()

	// Missing sequence.
	_, err = r.GetSequence("nonexistent")
	if err == nil {
		t.Error("expected error for missing sequence")
	}

	// Missing .fai.
	_, err = NewIndexedFastaReader(filepath.Join(dir, "nofai.fa"))
	if err == nil {
		t.Error("expected error when .fai is missing")
	}
}

func TestIndexedFastaReaderWithCRAMTestdata(t *testing.T) {
	// Use the existing CRAM test reference if available.
	fastaPath := "../htsio/cram/testdata/ref.fa"
	if _, err := os.Stat(fastaPath); err != nil {
		t.Skip("CRAM testdata not available")
	}
	faiPath := fastaPath + ".fai"
	if _, err := os.Stat(faiPath); err != nil {
		t.Skip("CRAM testdata .fai not available")
	}

	r, err := NewIndexedFastaReader(fastaPath)
	if err != nil {
		t.Fatalf("NewIndexedFastaReader: %v", err)
	}
	defer r.Close()

	// Verify chr1 length.
	l, ok := r.SequenceLength("chr1")
	if !ok {
		t.Fatal("chr1 not found")
	}
	if l != 100000 {
		t.Errorf("chr1 length = %d, want 100000", l)
	}

	// Read a small range and verify it's uppercase ACGTN.
	seq, err := r.GetSequenceRange("chr1", 0, 100)
	if err != nil {
		t.Fatalf("GetSequenceRange: %v", err)
	}
	if len(seq) != 100 {
		t.Fatalf("got %d bases, want 100", len(seq))
	}
	for i, b := range seq {
		if b != 'A' && b != 'C' && b != 'G' && b != 'T' && b != 'N' {
			t.Errorf("base %d: got %c, want ACGTN", i, b)
			break
		}
	}

	// Verify a range matches full sequence slice.
	full, err := r.GetSequence("chr1")
	if err != nil {
		t.Fatalf("GetSequence: %v", err)
	}
	partial, err := r.GetSequenceRange("chr1", 5000, 6000)
	if err != nil {
		t.Fatalf("GetSequenceRange: %v", err)
	}
	if !bytes.Equal(partial, full[5000:6000]) {
		t.Error("range [5000:6000] does not match full sequence slice")
	}
}
