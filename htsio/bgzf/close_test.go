package bgzf

import (
	"io"
	"os"
	"path/filepath"
	"runtime"
	"testing"
)

// openFDs counts this process's open file descriptors, or skips where that
// cannot be read (it needs /proc).
func openFDs(t *testing.T) int {
	t.Helper()
	if runtime.GOOS != "linux" {
		t.Skip("counting descriptors needs /proc")
	}
	ents, err := os.ReadDir("/proc/self/fd")
	if err != nil {
		t.Skipf("cannot read /proc/self/fd: %v", err)
	}
	return len(ents)
}

func writeTestBgzf(t *testing.T) string {
	t.Helper()
	data, _ := writeBgzfBlocks([][]byte{[]byte("hello, "), []byte("world\n")})
	path := filepath.Join(t.TempDir(), "x.gz")
	if err := os.WriteFile(path, data, 0o644); err != nil {
		t.Fatal(err)
	}
	return path
}

func TestIndexedReaderCloseClosesTheFileItOpened(t *testing.T) {
	r, err := OpenIndexedReader(writeTestBgzf(t))
	if err != nil {
		t.Fatal(err)
	}
	f := r.rs.(*os.File)
	if err := r.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	if _, err := f.Read(make([]byte, 1)); err == nil {
		t.Error("the file OpenIndexedReader opened is still open after Close")
	}
	if err := r.Close(); err != nil {
		t.Errorf("a second Close: %v", err)
	}
	if _, err := r.Read(make([]byte, 1)); err != errReaderClosed {
		t.Errorf("Read after Close: %v, want errReaderClosed", err)
	}
	// A seek clears the sticky error, so it must check for Close itself.
	if err := r.SeekToVirtualOffset(0); err != errReaderClosed {
		t.Errorf("SeekToVirtualOffset after Close: %v, want errReaderClosed", err)
	}
}

func TestIndexedReaderCloseLeavesACallersSourceOpen(t *testing.T) {
	f, err := os.Open(writeTestBgzf(t))
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	r := NewIndexedReader(f)
	if err := r.SeekToVirtualOffset(0); err != nil {
		t.Fatal(err)
	}
	if err := r.Close(); err != nil {
		t.Fatal(err)
	}
	if _, err := f.Seek(0, io.SeekStart); err != nil {
		t.Fatalf("Close closed a ReadSeeker the caller owns: %v", err)
	}
	if err := r.SeekToVirtualOffset(0); err != errReaderClosed {
		t.Errorf("the reader still works after Close: %v", err)
	}
}

// The leak itself: without Close, every OpenIndexedReader held a descriptor
// until the garbage collector finalized its file.
func TestOpenIndexedReaderDoesNotLeakDescriptors(t *testing.T) {
	path := writeTestBgzf(t)
	before := openFDs(t)
	for i := 0; i < 100; i++ {
		r, err := OpenIndexedReader(path)
		if err != nil {
			t.Fatal(err)
		}
		if _, err := io.ReadAll(r); err != nil {
			t.Fatal(err)
		}
		if err := r.Close(); err != nil {
			t.Fatal(err)
		}
	}
	if after := openFDs(t); after > before {
		t.Errorf("open descriptors grew from %d to %d over 100 open/close cycles", before, after)
	}
}
