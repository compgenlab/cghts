package bam

import (
	"bytes"
	"fmt"
	"testing"

	"github.com/compgenlab/hts/htsio"
)

func bamTestHeader() *htsio.SamHeader {
	h := htsio.NewSamHeader()
	h.AddLine("@HD\tVN:1.6\tSO:coordinate")
	h.AddLine("@SQ\tSN:chr1\tLN:1000")
	return h
}

// TestBamWriterRejectsCigarSeqMismatch verifies the BAM writer rejects a record
// whose CIGAR query length disagrees with its SEQ length (malformed per spec),
// synchronously at Write time.
func TestBamWriterRejectsCigarSeqMismatch(t *testing.T) {
	var buf bytes.Buffer
	w := NewWriterFromWriter(&buf, bamTestHeader())
	defer w.Close()
	bad := &htsio.SamRecord{
		ReadName: "r", Flag: 0, RefName: "chr1", Pos: 1, MapQ: 60,
		Cigar: "10M", RefNext: "*", Seq: "ACGT", Qual: "IIII", // 10M needs 10 bases
	}
	if err := w.Write(bad); err == nil {
		t.Fatal("expected CIGAR/SEQ mismatch error, got nil")
	}
}

// failWriter fails every write, simulating a broken sink (disk full / pipe).
type failWriter struct{}

func (failWriter) Write(p []byte) (int, error) { return 0, fmt.Errorf("injected write error") }

// TestBamWriterPropagatesSinkError verifies that a failing underlying writer
// surfaces an error through Write or Close rather than being silently
// swallowed, so a broken sink never looks like a successful write.
func TestBamWriterPropagatesSinkError(t *testing.T) {
	w := NewWriterFromWriter(failWriter{}, bamTestHeader())
	rec := &htsio.SamRecord{
		ReadName: "r", Flag: 0, RefName: "chr1", Pos: 1, MapQ: 60,
		Cigar: "4M", RefNext: "*", Seq: "ACGT", Qual: "IIII",
	}
	werr := w.Write(rec)
	cerr := w.Close()
	if werr == nil && cerr == nil {
		t.Fatal("expected a write or close error from a failing sink, got nil from both")
	}
}
