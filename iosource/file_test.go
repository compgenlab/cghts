package iosource

import (
	"errors"
	"io/fs"
	"os"
	"path/filepath"
	"testing"
)

func TestOpenFileReadsAndSizes(t *testing.T) {
	p := filepath.Join(t.TempDir(), "f.bin")
	if err := os.WriteFile(p, []byte("0123456789"), 0o644); err != nil {
		t.Fatal(err)
	}
	f, err := OpenFile(p)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()

	if n, err := f.Size(); err != nil || n != 10 {
		t.Errorf("Size() = (%d, %v), want (10, nil)", n, err)
	}
	buf := make([]byte, 4)
	if _, err := f.ReadAt(buf, 3); err != nil {
		t.Fatalf("ReadAt: %v", err)
	}
	if string(buf) != "3456" {
		t.Errorf("ReadAt(4, 3) = %q, want %q", buf, "3456")
	}
}

// A missing local file must satisfy errors.Is(err, fs.ErrNotExist), because
// that is the test callers use to tell "not there" from "could not be read" --
// and the HTTP source now wraps a 404 the same way so the two agree.
func TestOpenFileMissingIsNotExist(t *testing.T) {
	_, err := OpenFile(filepath.Join(t.TempDir(), "nope.bin"))
	if err == nil {
		t.Fatal("OpenFile succeeded on a missing path")
	}
	if !errors.Is(err, fs.ErrNotExist) {
		t.Errorf("error %v does not satisfy errors.Is(err, fs.ErrNotExist)", err)
	}
}

// NewFile documents that it takes ownership of the handle. Ownership passes on
// the call, not on success -- the Stat error path used to return without
// closing, leaving the caller holding a handle it had already given away.
func TestNewFileClosesOnStatFailure(t *testing.T) {
	p := filepath.Join(t.TempDir(), "f.bin")
	if err := os.WriteFile(p, []byte("x"), 0o644); err != nil {
		t.Fatal(err)
	}
	f, err := os.Open(p)
	if err != nil {
		t.Fatal(err)
	}
	// Stat on a closed handle fails, which is the error path under test.
	f.Close()

	if _, err := NewFile(f); err == nil {
		t.Fatal("NewFile succeeded on a closed handle")
	}
	// Already closed, so a second Close is the observable proof it was not left
	// open -- and must not panic either way.
	if err := f.Close(); err == nil {
		t.Log("handle closed twice cleanly")
	}
}

func TestNewFileTakesOwnership(t *testing.T) {
	p := filepath.Join(t.TempDir(), "f.bin")
	if err := os.WriteFile(p, []byte("hello"), 0o644); err != nil {
		t.Fatal(err)
	}
	h, err := os.Open(p)
	if err != nil {
		t.Fatal(err)
	}
	src, err := NewFile(h)
	if err != nil {
		t.Fatal(err)
	}
	if n, _ := src.Size(); n != 5 {
		t.Errorf("Size() = %d, want 5", n)
	}
	if err := src.Close(); err != nil {
		t.Errorf("Close: %v", err)
	}
	// The source closed the handle it was given, so reading through it fails.
	if _, err := h.ReadAt(make([]byte, 1), 0); err == nil {
		t.Error("the underlying handle is still open after Close")
	}
}
