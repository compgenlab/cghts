package bam

import (
	"bytes"
	"encoding/binary"
	"testing"

	"github.com/compgen-io/cgkit/htsio/bgzf"
)

// bgzfReaderFrom wraps raw bytes in a single BGZF block and returns a reader
// over it, matching how a real BAM stream is framed.
func bgzfReaderFrom(t *testing.T, raw []byte) *bgzf.Reader {
	t.Helper()
	var buf bytes.Buffer
	w := bgzf.NewWriter(&buf)
	if _, err := w.Write(raw); err != nil {
		t.Fatalf("bgzf write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("bgzf close: %v", err)
	}
	return bgzf.NewReader(bytes.NewReader(buf.Bytes()))
}

// le32 encodes int32 values little-endian (BAM header field encoding).
func le32(vals ...int32) []byte {
	b := &bytes.Buffer{}
	for _, v := range vals {
		_ = binary.Write(b, binary.LittleEndian, v)
	}
	return b.Bytes()
}

// TestReadHeaderMalformed asserts the BAM header parser rejects a zero-length
// reference name (which would otherwise panic on nameBuf[:nameLen-1]) and an
// implausibly large reference count (which would otherwise pre-allocate a huge
// slice), returning an error in both cases rather than panicking or OOMing.
func TestReadHeaderMalformed(t *testing.T) {
	// magic + l_text=0 + n_ref=1 + ref0 l_name=0
	zeroName := append([]byte("BAM\x01"), le32(0, 1, 0)...)
	b := &Reader{r: bgzfReaderFrom(t, zeroName)}
	if err := b.readHeader(); err == nil {
		t.Fatal("expected error for l_name=0, got nil")
	}

	// magic + l_text=0 + huge n_ref, then EOF
	hugeNRef := append([]byte("BAM\x01"), le32(0, 0x7FFFFFFF)...)
	b2 := &Reader{r: bgzfReaderFrom(t, hugeNRef)}
	if err := b2.readHeader(); err == nil {
		t.Fatal("expected error for huge n_ref, got nil")
	}

	// magic + huge l_text, then EOF
	hugeText := append([]byte("BAM\x01"), le32(0x7FFFFFFF)...)
	b3 := &Reader{r: bgzfReaderFrom(t, hugeText)}
	if err := b3.readHeader(); err == nil {
		t.Fatal("expected error for huge l_text, got nil")
	}
}

// TestReadHeaderValid confirms a minimal well-formed header (one reference)
// still parses after the hardening.
func TestReadHeaderValid(t *testing.T) {
	raw := []byte("BAM\x01")
	raw = append(raw, le32(0, 1)...)          // l_text=0, n_ref=1
	raw = append(raw, le32(int32(len("c1\x00")))...) // l_name (includes NUL)
	raw = append(raw, []byte("c1\x00")...)    // name
	raw = append(raw, le32(1000)...)          // l_ref

	b := &Reader{r: bgzfReaderFrom(t, raw)}
	if err := b.readHeader(); err != nil {
		t.Fatalf("valid header rejected: %v", err)
	}
	if len(b.refs) != 1 || b.refs[0].name != "c1" || b.refs[0].length != 1000 {
		t.Fatalf("unexpected refs: %+v", b.refs)
	}
}
