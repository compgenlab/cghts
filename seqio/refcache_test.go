package seqio

import (
	"os"
	"path/filepath"
	"testing"
)

func TestExpandCachePattern(t *testing.T) {
	tests := []struct {
		pattern string
		md5     string
		want    string
	}{
		{"/cache/%2s/%2s/%s", "abcdef1234567890", "/cache/ab/cd/abcdef1234567890"},
		{"/cache/%s", "abcdef", "/cache/abcdef"},
		{"/cache/%3s/%s", "abcdef", "/cache/abc/abcdef"},
		{"/simple/path", "abcdef", "/simple/path"},
	}
	for _, tc := range tests {
		got := expandCachePattern(tc.pattern, tc.md5)
		if got != tc.want {
			t.Errorf("expandCachePattern(%q, %q) = %q, want %q", tc.pattern, tc.md5, got, tc.want)
		}
	}
}

func TestRefCacheReaderLocalCache(t *testing.T) {
	// Set up a temp REF_CACHE directory with a sequence file.
	dir := t.TempDir()
	md5 := "d41d8cd98f00b204e9800998ecf8427e"

	// Create cache hierarchy: %2s/%2s/%s
	subdir := filepath.Join(dir, md5[:2], md5[2:4])
	if err := os.MkdirAll(subdir, 0755); err != nil {
		t.Fatal(err)
	}
	seqData := "ACGTACGTNNTTGGCC"
	if err := os.WriteFile(filepath.Join(subdir, md5), []byte(seqData), 0644); err != nil {
		t.Fatal(err)
	}

	// Set env vars.
	pattern := filepath.Join(dir, "%2s/%2s/%s")
	t.Setenv("REF_CACHE", pattern)
	t.Setenv("REF_PATH", "")

	md5s := map[string]string{"chr1": md5}
	lengths := map[string]int{"chr1": 16}
	names := []string{"chr1"}

	rr, err := NewRefCacheReader(md5s, lengths, names)
	if err != nil {
		t.Fatalf("NewRefCacheReader: %v", err)
	}
	defer rr.Close()

	// Full sequence.
	seq, err := rr.GetSequence("chr1")
	if err != nil {
		t.Fatalf("GetSequence: %v", err)
	}
	if string(seq) != "ACGTACGTNNTTGGCC" {
		t.Errorf("GetSequence = %q, want ACGTACGTNNTTGGCC", string(seq))
	}

	// Range.
	sub, err := rr.GetSequenceRange("chr1", 0, 4)
	if err != nil {
		t.Fatalf("GetSequenceRange: %v", err)
	}
	if string(sub) != "ACGT" {
		t.Errorf("GetSequenceRange(0,4) = %q, want ACGT", string(sub))
	}

	// Missing sequence.
	_, err = rr.GetSequenceRange("chrX", 0, 10)
	if err == nil {
		t.Error("expected error for missing sequence")
	}
}

func TestRefCacheReaderRefPath(t *testing.T) {
	// Test REF_PATH with a plain directory.
	dir := t.TempDir()
	md5 := "abc123def456"
	seqData := "GGCCTTAA"
	if err := os.WriteFile(filepath.Join(dir, md5), []byte(seqData), 0644); err != nil {
		t.Fatal(err)
	}

	t.Setenv("REF_PATH", dir)
	t.Setenv("REF_CACHE", "")

	md5s := map[string]string{"seq1": md5}
	lengths := map[string]int{"seq1": 8}
	names := []string{"seq1"}

	rr, err := NewRefCacheReader(md5s, lengths, names)
	if err != nil {
		t.Fatalf("NewRefCacheReader: %v", err)
	}
	defer rr.Close()

	seq, err := rr.GetSequence("seq1")
	if err != nil {
		t.Fatalf("GetSequence: %v", err)
	}
	if string(seq) != "GGCCTTAA" {
		t.Errorf("GetSequence = %q", string(seq))
	}
}

func TestRefCacheReaderNoEnvVars(t *testing.T) {
	t.Setenv("REF_PATH", "")
	t.Setenv("REF_CACHE", "")

	_, err := NewRefCacheReader(
		map[string]string{"chr1": "abc"},
		map[string]int{"chr1": 10},
		[]string{"chr1"},
	)
	if err == nil {
		t.Error("expected error when env vars not set")
	}
}
