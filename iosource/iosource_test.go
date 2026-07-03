package iosource

import (
	"bytes"
	"io"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"testing"
	"time"
)

var zeroTime time.Time

func makeData(n int) []byte {
	b := make([]byte, n)
	for i := range b {
		b[i] = byte(i % 251) // 251 is prime → no short period aligned to reads
	}
	return b
}

func TestFileByteSource(t *testing.T) {
	data := makeData(4096)
	path := filepath.Join(t.TempDir(), "data.bin")
	if err := os.WriteFile(path, data, 0o644); err != nil {
		t.Fatal(err)
	}

	src, err := OpenFile(path)
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()

	if sz, err := src.Size(); err != nil || sz != int64(len(data)) {
		t.Fatalf("Size() = %d, %v; want %d", sz, err, len(data))
	}

	buf := make([]byte, 100)
	n, err := src.ReadAt(buf, 1000)
	if err != nil {
		t.Fatalf("ReadAt: %v", err)
	}
	if !bytes.Equal(buf[:n], data[1000:1100]) {
		t.Fatalf("ReadAt returned wrong bytes")
	}

	// Reading past EOF returns io.EOF with the available bytes.
	tail := make([]byte, 50)
	n, err = src.ReadAt(tail, int64(len(data)-10))
	if err != io.EOF {
		t.Fatalf("ReadAt at EOF: err = %v, want io.EOF", err)
	}
	if n != 10 || !bytes.Equal(tail[:n], data[len(data)-10:]) {
		t.Fatalf("ReadAt at EOF returned %d bytes, wrong content", n)
	}
}

func TestHTTPRange(t *testing.T) {
	data := makeData(1 << 16)
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		// http.ServeContent handles Range + HEAD + Content-Range correctly.
		http.ServeContent(w, r, "data.bin", zeroTime, bytes.NewReader(data))
	}))
	defer srv.Close()

	src, err := NewHTTPRange(srv.URL)
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()

	if sz, err := src.Size(); err != nil || sz != int64(len(data)) {
		t.Fatalf("Size() = %d, %v; want %d", sz, err, len(data))
	}

	for _, tc := range []struct{ off, n int }{{0, 10}, {12345, 4096}, {len(data) - 5, 5}} {
		buf := make([]byte, tc.n)
		got, err := src.ReadAt(buf, int64(tc.off))
		if err != nil && err != io.EOF {
			t.Fatalf("ReadAt(%d,%d): %v", tc.off, tc.n, err)
		}
		if !bytes.Equal(buf[:got], data[tc.off:tc.off+got]) {
			t.Fatalf("ReadAt(%d,%d) wrong bytes", tc.off, tc.n)
		}
	}

	// Reading past the end returns io.EOF.
	buf := make([]byte, 100)
	n, err := src.ReadAt(buf, int64(len(data)-10))
	if err != io.EOF {
		t.Fatalf("ReadAt past end: err = %v, want io.EOF", err)
	}
	if n != 10 {
		t.Fatalf("ReadAt past end: n = %d, want 10", n)
	}
}

// TestHTTPRangeIgnoredRange verifies correctness against a server that ignores
// Range headers and always returns 200 with the whole body.
func TestHTTPRangeIgnoredRange(t *testing.T) {
	data := makeData(8192)
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.Header().Set("Content-Type", "application/octet-stream")
		w.WriteHeader(http.StatusOK)
		w.Write(data)
	}))
	defer srv.Close()

	src, err := NewHTTPRange(srv.URL)
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()

	buf := make([]byte, 200)
	n, err := src.ReadAt(buf, 5000)
	if err != nil && err != io.EOF {
		t.Fatalf("ReadAt: %v", err)
	}
	if !bytes.Equal(buf[:n], data[5000:5000+n]) {
		t.Fatalf("ReadAt against range-ignoring server returned wrong bytes")
	}
}

func TestResolveSibling(t *testing.T) {
	dir := t.TempDir()
	base := filepath.Join(dir, "reads.bam")
	if err := os.WriteFile(base+".csi", []byte("csi"), 0o644); err != nil {
		t.Fatal(err)
	}

	rc, matched, err := ResolveSibling(base, []string{".bai", ".csi"}, FileSibling)
	if err != nil {
		t.Fatalf("ResolveSibling: %v", err)
	}
	defer rc.Close()
	if matched != ".csi" {
		t.Fatalf("matched = %q, want .csi", matched)
	}

	if _, _, err := ResolveSibling(base, []string{".bai"}, FileSibling); err == nil {
		t.Fatalf("expected error when no sibling exists")
	}
}
