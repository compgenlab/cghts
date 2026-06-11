package cram

import (
	"bytes"
	"testing"

	"github.com/compgen-io/cgkit/htsio"
)

func minimalHeader() *htsio.SamHeader {
	h := htsio.NewSamHeader()
	h.AddLine("@HD\tVN:1.6")
	h.AddLine("@SQ\tSN:chr1\tLN:1000")
	return h
}

// TestWriterCloseIdempotent verifies a second Close() does not write a second
// EOF container or error — it returns the first call's result.
func TestWriterCloseIdempotent(t *testing.T) {
	var buf bytes.Buffer
	w, err := NewWriterFromWriter(&buf, minimalHeader())
	if err != nil {
		t.Fatalf("NewWriterFromWriter: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("first Close: %v", err)
	}
	n := buf.Len()
	if err := w.Close(); err != nil {
		t.Fatalf("second Close: %v", err)
	}
	if buf.Len() != n {
		t.Fatalf("second Close wrote %d extra bytes (double EOF / double write)", buf.Len()-n)
	}
}

// TestWriterWriteAfterClose verifies writes after Close are rejected.
func TestWriterWriteAfterClose(t *testing.T) {
	var buf bytes.Buffer
	w, err := NewWriterFromWriter(&buf, minimalHeader())
	if err != nil {
		t.Fatalf("NewWriterFromWriter: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	rec := &htsio.SamRecord{ReadName: "r", Flag: 4, RefName: "*", Cigar: "*", RefNext: "*", Seq: "ACGT", Qual: "IIII"}
	if err := w.Write(rec); err == nil {
		t.Fatal("expected error writing to a closed writer, got nil")
	}
}

// TestWriterRejectsCigarSeqMismatch verifies a record whose CIGAR query length
// disagrees with its SEQ length is rejected instead of silently losing bases.
func TestWriterRejectsCigarSeqMismatch(t *testing.T) {
	var buf bytes.Buffer
	w, err := NewWriterFromWriter(&buf, minimalHeader())
	if err != nil {
		t.Fatalf("NewWriterFromWriter: %v", err)
	}
	defer w.Close()
	bad := &htsio.SamRecord{
		ReadName: "r", Flag: 0, RefName: "chr1", Pos: 1, MapQ: 60,
		Cigar: "10M", RefNext: "*", Seq: "ACGT", Qual: "IIII", // 10M needs 10 bases, has 4
	}
	if err := w.Write(bad); err == nil {
		t.Fatal("expected CIGAR/SEQ mismatch error, got nil")
	}
}
