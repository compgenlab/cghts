package bam

import (
	"bytes"
	"fmt"
	"io"
	"strings"
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

// TestBamWriterRejectsMalformedCigar verifies the BAM writer rejects a record
// whose CIGAR is malformed, synchronously at Write time. encodeCigar runs on the
// async encode path where an error cannot reach the caller, and it used to encode
// an unrecognized operation as M, silently corrupting the record.
//
// Every case here is one that TestBamWriterRejectsCigarSeqMismatch's check lets
// through, so the subtest asserts that precondition rather than assuming it:
// CigarQueryLen skips what it does not recognize, so "10Z" implies query length
// 0, and ValidateCigarSeq is a no-op when SEQ is "*".
func TestBamWriterRejectsMalformedCigar(t *testing.T) {
	cases := []struct{ name, cigar, seq string }{
		{"invalid operation", "10Z", "*"},
		{"trailing length with no operation", "10M5", "ACGTACGTAC"},
		{"B is not a SAM operation", "10B", "*"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if err := htsio.ValidateCigarSeq(c.cigar, c.seq); err != nil {
				t.Fatalf("precondition: ValidateCigarSeq already rejects %q (%v), so this case "+
					"does not exercise the CIGAR parse check", c.cigar, err)
			}
			var buf bytes.Buffer
			w := NewWriterFromWriter(&buf, bamTestHeader())
			defer w.Close()
			bad := &htsio.SamRecord{
				ReadName: "r", Flag: 0, RefName: "chr1", Pos: 1, MapQ: 60,
				Cigar: c.cigar, RefNext: "*", Seq: c.seq, Qual: "*",
			}
			if err := w.Write(bad); err == nil {
				t.Errorf("Write accepted malformed CIGAR %q, want an error", c.cigar)
			}
		})
	}
}

// TestBamCigarRoundTrip verifies that every CIGAR operation survives a write and
// read back unchanged, pinning the operation -> opcode mapping in encodeCigar
// against the decode table in reader.go.
//
// Without this, transposing two opcodes in cigarOpCode (say I and D) passes the
// entire suite: the malformed-CIGAR tests only exercise rejection, and nothing
// else asserts an encoded CIGAR's value.
func TestBamCigarRoundTrip(t *testing.T) {
	cigars := []string{
		"10M",
		"5M1I4M",             // I
		"5M1D5M",             // D
		"5M100N5M",           // N
		"3S4M3S",             // S
		"3H4M3H",             // H
		"4M2P4M",             // P
		"5=5X",               // = and X
		"2H3S4M1I1D2X2=3S2H", // every operation at once
		"*",
	}
	for _, cigar := range cigars {
		t.Run(cigar, func(t *testing.T) {
			var buf bytes.Buffer
			w := NewWriterFromWriter(&buf, bamTestHeader())
			seq := "*"
			if n := htsio.CigarQueryLen(cigar); n > 0 {
				seq = strings.Repeat("A", n)
			}
			want := &htsio.SamRecord{
				ReadName: "r", Flag: 0, RefName: "chr1", Pos: 1, MapQ: 60,
				Cigar: cigar, RefNext: "*", Seq: seq, Qual: "*",
			}
			if err := w.Write(want); err != nil {
				t.Fatalf("Write(%q): %v", cigar, err)
			}
			if err := w.Close(); err != nil {
				t.Fatalf("Close: %v", err)
			}

			r, err := NewReader(io.NopCloser(&buf), "test.bam", nil)
			if err != nil {
				t.Fatalf("NewReader: %v", err)
			}
			defer r.Close()

			n := 0
			for got, err := range r.Records() {
				if err != nil {
					t.Fatalf("Records: %v", err)
				}
				if got.Cigar != cigar {
					t.Errorf("round trip of %q gave %q", cigar, got.Cigar)
				}
				n++
			}
			if n != 1 {
				t.Errorf("read %d records, want 1", n)
			}
		})
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
