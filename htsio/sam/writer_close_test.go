package sam

import (
	"bytes"
	"testing"

	"github.com/compgen-io/cgkit/htsio"
)

// TestSamWriterCloseIdempotent verifies a second Close() is a no-op that returns
// the first call's result, aligning the SAM writer with the BAM/CRAM writers.
func TestSamWriterCloseIdempotent(t *testing.T) {
	h := htsio.NewSamHeader()
	h.AddLine("@HD\tVN:1.6")
	var buf bytes.Buffer
	w := NewWriterFromWriter(&buf, h)

	if err := w.Close(); err != nil {
		t.Fatalf("first Close: %v", err)
	}
	n := buf.Len()
	if err := w.Close(); err != nil {
		t.Fatalf("second Close: %v", err)
	}
	if buf.Len() != n {
		t.Fatalf("second Close wrote %d extra bytes", buf.Len()-n)
	}
}
