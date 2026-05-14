package seqio

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestOpenReferenceIndexed(t *testing.T) {
	dir := t.TempDir()
	fastaPath := createTestFasta(t, dir) // creates .fa + .fai

	ref, err := OpenReference(fastaPath)
	if err != nil {
		t.Fatalf("OpenReference: %v", err)
	}
	defer ref.Close()

	// Should use IndexedFastaReader.
	if _, ok := ref.(*IndexedFastaReader); !ok {
		t.Errorf("expected IndexedFastaReader, got %T", ref)
	}

	// Verify interface works.
	names := ref.Names()
	if len(names) != 2 || names[0] != "seq1" || names[1] != "seq2" {
		t.Errorf("Names() = %v", names)
	}

	l, ok := ref.SequenceLength("seq1")
	if !ok || l != 2000 {
		t.Errorf("SequenceLength(seq1) = %d, %v", l, ok)
	}

	seq, err := ref.GetSequenceRange("seq1", 0, 10)
	if err != nil {
		t.Fatalf("GetSequenceRange: %v", err)
	}
	if string(seq) != "ACGTACGTNN" {
		t.Errorf("got %q, want ACGTACGTNN", string(seq))
	}
}

func TestOpenReferenceFallback(t *testing.T) {
	dir := t.TempDir()

	// Create a FASTA without .fai.
	var buf strings.Builder
	buf.WriteString(">chrX\n")
	buf.WriteString("ACGTACGT\n")
	buf.WriteString(">chrY\n")
	buf.WriteString("TTGGCCAA\n")

	fastaPath := filepath.Join(dir, "noidx.fa")
	if err := os.WriteFile(fastaPath, []byte(buf.String()), 0644); err != nil {
		t.Fatal(err)
	}

	ref, err := OpenReference(fastaPath)
	if err != nil {
		t.Fatalf("OpenReference: %v", err)
	}
	defer ref.Close()

	// Should use inMemoryFasta fallback.
	if _, ok := ref.(*inMemoryFasta); !ok {
		t.Errorf("expected inMemoryFasta, got %T", ref)
	}

	names := ref.Names()
	if len(names) != 2 || names[0] != "chrX" || names[1] != "chrY" {
		t.Errorf("Names() = %v", names)
	}

	l, ok := ref.SequenceLength("chrX")
	if !ok || l != 8 {
		t.Errorf("SequenceLength(chrX) = %d, %v", l, ok)
	}

	seq, err := ref.GetSequence("chrY")
	if err != nil {
		t.Fatalf("GetSequence: %v", err)
	}
	if string(seq) != "TTGGCCAA" {
		t.Errorf("got %q, want TTGGCCAA", string(seq))
	}

	// Range query.
	partial, err := ref.GetSequenceRange("chrX", 2, 6)
	if err != nil {
		t.Fatalf("GetSequenceRange: %v", err)
	}
	if string(partial) != "GTAC" {
		t.Errorf("got %q, want GTAC", string(partial))
	}
}

func TestOpenReferenceNotFound(t *testing.T) {
	_, err := OpenReference("/nonexistent/ref.fa")
	if err == nil {
		t.Error("expected error for missing file")
	}
}
