package cram

import (
	"bytes"
	"testing"
)

// TestReadSizedBytes checks the shared sized-block reader rejects a negative
// size, reports truncation (rather than OOM-allocating) for an implausibly huge
// size backed by little data, and still reads an exact-length block.
func TestReadSizedBytes(t *testing.T) {
	if _, err := readSizedBytes(bytes.NewReader(nil), -1); err == nil {
		t.Fatal("expected error for negative size")
	}
	// A 1 GiB declared size with 3 bytes available must fail as a truncation,
	// and must do so without allocating 1 GiB (LimitReader bounds the buffer).
	if _, err := readSizedBytes(bytes.NewReader([]byte("abc")), 1<<30); err == nil {
		t.Fatal("expected truncation error for huge size with short input")
	}
	got, err := readSizedBytes(bytes.NewReader([]byte("hello")), 5)
	if err != nil || string(got) != "hello" {
		t.Fatalf("readSizedBytes exact: got %q err %v", got, err)
	}
}

// TestReadITF8ArrayMalformedNoOOM encodes a huge array length followed by no
// element data and asserts the reader errors (growing by append) instead of
// pre-allocating from the attacker-controlled count.
func TestReadITF8ArrayMalformedNoOOM(t *testing.T) {
	var buf bytes.Buffer
	if err := writeITF8(&buf, 0x0FFFFFFF); err != nil {
		t.Fatalf("writeITF8: %v", err)
	}
	if _, err := readITF8Array(bytes.NewReader(buf.Bytes())); err == nil {
		t.Fatal("expected error reading huge ITF8 array with no element data")
	}
}

// TestReadITF8ArrayValid confirms a small, well-formed array still round-trips.
func TestReadITF8ArrayValid(t *testing.T) {
	var buf bytes.Buffer
	if err := writeITF8Array(&buf, []int32{10, 20, 30}); err != nil {
		t.Fatalf("writeITF8Array: %v", err)
	}
	vals, err := readITF8Array(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("readITF8Array: %v", err)
	}
	if len(vals) != 3 || vals[0] != 10 || vals[2] != 30 {
		t.Fatalf("readITF8Array round-trip = %v", vals)
	}
}
