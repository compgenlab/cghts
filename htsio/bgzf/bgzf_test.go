package bgzf

import (
	"bytes"
	"crypto/rand"
	"io"
	"testing"
)

func TestRoundTrip(t *testing.T) {
	original := []byte("Hello, BGZF! This is a round-trip test of the bgzf package.")

	var buf bytes.Buffer
	w := NewWriter(&buf)
	if _, err := w.Write(original); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	r := NewReader(&buf)
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if !bytes.Equal(got, original) {
		t.Errorf("round-trip mismatch: got %q, want %q", got, original)
	}
}

func TestRoundTripLarge(t *testing.T) {
	// Write enough data to span multiple BGZF blocks.
	original := make([]byte, MaxUncompressedSize*3+1234)
	for i := range original {
		original[i] = byte(i % 251) // deterministic pattern
	}

	var buf bytes.Buffer
	w := NewWriter(&buf)
	if _, err := w.Write(original); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	r := NewReader(&buf)
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if !bytes.Equal(got, original) {
		t.Errorf("round-trip mismatch: lengths got=%d want=%d", len(got), len(original))
	}
}

func TestRoundTripRandom(t *testing.T) {
	original := make([]byte, 100000)
	if _, err := rand.Read(original); err != nil {
		t.Fatalf("rand.Read: %v", err)
	}

	var buf bytes.Buffer
	w := NewWriter(&buf)
	if _, err := w.Write(original); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	r := NewReader(&buf)
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if !bytes.Equal(got, original) {
		t.Errorf("round-trip mismatch for random data")
	}
}

func TestParallelRoundTrip(t *testing.T) {
	original := make([]byte, MaxUncompressedSize*5+999)
	for i := range original {
		original[i] = byte(i % 251)
	}

	var buf bytes.Buffer
	w := NewParallelWriter(&buf, 4)
	if _, err := w.Write(original); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	r := NewReader(&buf)
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if !bytes.Equal(got, original) {
		t.Errorf("parallel round-trip mismatch: got %d bytes, want %d", len(got), len(original))
	}
}

func TestParallelRoundTripRandom(t *testing.T) {
	original := make([]byte, 200000)
	if _, err := rand.Read(original); err != nil {
		t.Fatalf("rand.Read: %v", err)
	}

	var buf bytes.Buffer
	w := NewParallelWriter(&buf, 4)
	if _, err := w.Write(original); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	r := NewReader(&buf)
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if !bytes.Equal(got, original) {
		t.Errorf("parallel round-trip mismatch for random data")
	}
}

func TestParallelEmpty(t *testing.T) {
	var buf bytes.Buffer
	w := NewParallelWriter(&buf, 4)
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	r := NewReader(&buf)
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(got) != 0 {
		t.Errorf("expected empty output, got %d bytes", len(got))
	}
}

func TestParallelSmallWrites(t *testing.T) {
	original := []byte("ACGTACGTACGTACGT")

	var buf bytes.Buffer
	w := NewParallelWriter(&buf, 2)
	for _, b := range original {
		if _, err := w.Write([]byte{b}); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	r := NewReader(&buf)
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if !bytes.Equal(got, original) {
		t.Errorf("got %q, want %q", got, original)
	}
}

func TestParallelThreadsZero(t *testing.T) {
	// threads=0 should use runtime.NumCPU() and not panic.
	var buf bytes.Buffer
	w := NewParallelWriter(&buf, 0)
	if _, err := w.Write([]byte("test")); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	r := NewReader(&buf)
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if string(got) != "test" {
		t.Errorf("got %q, want %q", got, "test")
	}
}

func TestParallelThreadsOne(t *testing.T) {
	// threads=1 should fall back to single-threaded mode.
	var buf bytes.Buffer
	w := NewParallelWriter(&buf, 1)
	if _, err := w.Write([]byte("single")); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	r := NewReader(&buf)
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if string(got) != "single" {
		t.Errorf("got %q, want %q", got, "single")
	}
}

func TestEmptyFile(t *testing.T) {
	// A valid BGZF file with no data is just the EOF block.
	var buf bytes.Buffer
	w := NewWriter(&buf)
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	r := NewReader(&buf)
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(got) != 0 {
		t.Errorf("expected empty output, got %d bytes", len(got))
	}
}

func TestReadByte(t *testing.T) {
	original := []byte{0xDE, 0xAD, 0xBE, 0xEF}

	var buf bytes.Buffer
	w := NewWriter(&buf)
	if _, err := w.Write(original); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	r := NewReader(&buf)
	for i, want := range original {
		got, err := r.ReadByte()
		if err != nil {
			t.Fatalf("ReadByte[%d]: %v", i, err)
		}
		if got != want {
			t.Errorf("ReadByte[%d]: got %#x, want %#x", i, got, want)
		}
	}
	// Next read should return EOF.
	_, err := r.ReadByte()
	if err != io.EOF {
		t.Errorf("expected io.EOF, got %v", err)
	}
}

func TestVirtualOffset(t *testing.T) {
	// Write two blocks: first is exactly MaxUncompressedSize, second has remainder.
	block1 := make([]byte, MaxUncompressedSize)
	for i := range block1 {
		block1[i] = 'A'
	}
	block2 := []byte("second block")

	var buf bytes.Buffer
	w := NewWriter(&buf)
	if _, err := w.Write(block1); err != nil {
		t.Fatalf("Write block1: %v", err)
	}
	if _, err := w.Write(block2); err != nil {
		t.Fatalf("Write block2: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	r := NewReader(&buf)

	// Before reading, virtual offset should be 0:0.
	vo := r.VirtualTell()
	if vo.BlockOffset() != 0 {
		t.Errorf("initial block offset: got %d, want 0", vo.BlockOffset())
	}
	if vo.WithinBlock() != 0 {
		t.Errorf("initial within-block: got %d, want 0", vo.WithinBlock())
	}

	// Read the first block entirely.
	discard := make([]byte, MaxUncompressedSize)
	if _, err := io.ReadFull(r, discard); err != nil {
		t.Fatalf("ReadFull block1: %v", err)
	}

	// Read one byte from block2 — this triggers loading the next block.
	b, err := r.ReadByte()
	if err != nil {
		t.Fatalf("ReadByte: %v", err)
	}
	if b != 's' {
		t.Errorf("first byte of block2: got %c, want s", b)
	}

	// Virtual offset should now point to block2, within-block=1.
	vo = r.VirtualTell()
	if vo.BlockOffset() == 0 {
		t.Error("after reading into block2, block offset should be non-zero")
	}
	if vo.WithinBlock() != 1 {
		t.Errorf("within-block after 1 byte: got %d, want 1", vo.WithinBlock())
	}
}

func TestVirtualOffsetEncoding(t *testing.T) {
	vo := NewVirtualOffset(12345, 678)
	if vo.BlockOffset() != 12345 {
		t.Errorf("BlockOffset: got %d, want 12345", vo.BlockOffset())
	}
	if vo.WithinBlock() != 678 {
		t.Errorf("WithinBlock: got %d, want 678", vo.WithinBlock())
	}
}

func TestInvalidMagic(t *testing.T) {
	data := []byte("this is not a bgzf file")
	r := NewReader(bytes.NewReader(data))
	_, err := r.ReadByte()
	if err == nil {
		t.Fatal("expected error for invalid magic")
	}
}

func TestEOFBlockOnly(t *testing.T) {
	r := NewReader(bytes.NewReader(bgzfEOFBlock))
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(got) != 0 {
		t.Errorf("expected empty output, got %d bytes", len(got))
	}
}

func TestSmallWrites(t *testing.T) {
	// Write one byte at a time to test buffering.
	original := []byte("ACGTACGTACGTACGT")

	var buf bytes.Buffer
	w := NewWriter(&buf)
	for _, b := range original {
		if _, err := w.Write([]byte{b}); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	r := NewReader(&buf)
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if !bytes.Equal(got, original) {
		t.Errorf("got %q, want %q", got, original)
	}
}

func TestFlush(t *testing.T) {
	var buf bytes.Buffer
	w := NewWriter(&buf)

	if _, err := w.Write([]byte("block1")); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}

	if _, err := w.Write([]byte("block2")); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	r := NewReader(&buf)
	got, err := io.ReadAll(r)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if string(got) != "block1block2" {
		t.Errorf("got %q, want %q", got, "block1block2")
	}
}
