package bgzf

import (
	"bytes"
	"io"
	"testing"
)

// A full block holds 65536 bytes, one more than a virtual offset's 16-bit
// within-block field can hold, so "the end of this block" has to be reported
// as "the start of the next one" -- as htslib does. It used to wrap to byte 0
// of the same block, and a tabix index built from it sent readers into the
// middle of a line.
func TestVirtualTellAtTheEndOfAFullBlock(t *testing.T) {
	full := bytes.Repeat([]byte("x"), MaxUncompressedSize)
	var buf bytes.Buffer
	w := NewWriter(&buf)
	if _, err := w.Write(full); err != nil {
		t.Fatal(err)
	}
	// The block is full: the writer flushes it at once, so its position is
	// already the next block's start rather than byte 65536 of this one.
	wt := w.VirtualTell()
	if wt.WithinBlock() != 0 || wt.BlockOffset() == 0 {
		t.Fatalf("writer after a full block: %d:%d, want next-block:0", wt.BlockOffset(), wt.WithinBlock())
	}
	if _, err := w.Write([]byte("y\n")); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	next := wt.BlockOffset() // where the second block starts

	r := NewReader(bytes.NewReader(buf.Bytes()))
	if _, err := io.ReadFull(r, make([]byte, MaxUncompressedSize)); err != nil {
		t.Fatal(err)
	}
	if got := r.VirtualTell(); got != NewVirtualOffset(next, 0) {
		t.Errorf("Reader after a full block: %d:%d, want %d:0", got.BlockOffset(), got.WithinBlock(), next)
	}

	ir := NewIndexedReader(bytes.NewReader(buf.Bytes()))
	if err := ir.SeekToVirtualOffset(0); err != nil {
		t.Fatal(err)
	}
	if _, err := io.ReadFull(ir, make([]byte, MaxUncompressedSize)); err != nil {
		t.Fatal(err)
	}
	if got := ir.VirtualTell(); got != NewVirtualOffset(next, 0) {
		t.Errorf("IndexedReader after a full block: %d:%d, want %d:0", got.BlockOffset(), got.WithinBlock(), next)
	}
	// And that offset reads what follows.
	if err := ir.SeekToVirtualOffset(NewVirtualOffset(next, 0)); err != nil {
		t.Fatal(err)
	}
	rest, _ := io.ReadAll(ir)
	if string(rest) != "y\n" {
		t.Errorf("reading from the reported offset: %q, want %q", rest, "y\n")
	}
}
