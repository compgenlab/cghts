package tabix

import (
	"fmt"
	"path/filepath"
	"testing"
)

// writeBigBED builds a BGZF+tabix BED spanning many 16 kb windows on two chroms,
// with enough records to fill multiple bgzf blocks (so some lines straddle block
// boundaries). Returns the path.
func writeBigBED(t *testing.T) string {
	t.Helper()
	dir := t.TempDir()
	path := filepath.Join(dir, "big.bed.gz")
	w := NewWriter(path, NewWriterOpts().BED().AutoIndex())
	for _, chrom := range []string{"chr1", "chr2"} {
		for i := 0; i < 4000; i++ {
			start := i * 50 // 0..200000 → ~12 windows of 16384
			if err := w.Write(fmt.Sprintf("%s\t%d\t%d\t%s_feat%d", chrom, start, start+20, chrom, i)); err != nil {
				t.Fatal(err)
			}
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	return path
}

// TestCacheEquivalence: the cache must be transparent — a cache-on reader returns
// byte-identical record sets to a cache-off reader for every query, including point
// queries, small ranges, window-crossing ranges, empty regions, and unknown refs.
// The multi-block fixture also exercises records straddling bgzf block boundaries
// (hydrate reuses the same scanner as the direct path, so straddlers reassemble).
func TestCacheEquivalence(t *testing.T) {
	path := writeBigBED(t)

	on, err := NewReader(path) // cache enabled (default)
	if err != nil {
		t.Fatal(err)
	}
	defer on.Close()
	off, err := NewReaderSize(path, 0) // cache disabled
	if err != nil {
		t.Fatal(err)
	}
	defer off.Close()

	// Unknown ref must error on both paths (parity preserved through the cache).
	if _, err := on.Query("chrX", 0, 1); err == nil {
		t.Error("cache-on: unknown ref should error")
	}
	if _, err := off.Query("chrX", 0, 1); err == nil {
		t.Error("cache-off: unknown ref should error")
	}

	type q struct {
		ref        string
		start, end int
	}
	var queries []q
	for _, ref := range []string{"chr1", "chr2"} {
		for _, s := range []int{0, 1, 49, 50, 100, 16380, 16384, 16390, 99950, 100000, 200000} {
			queries = append(queries,
				q{ref, s, s + 1},     // point
				q{ref, s, s + 100},   // small range
				q{ref, s, s + 20000}) // spans >1 window (bypasses cache)
		}
	}

	for _, qq := range queries {
		want := queryLines(t, off, qq.ref, qq.start, qq.end)
		got := queryLines(t, on, qq.ref, qq.start, qq.end)
		if len(want) != len(got) {
			t.Fatalf("query %v: len got=%d want=%d", qq, len(got), len(want))
		}
		for i := range want {
			if got[i] != want[i] {
				t.Fatalf("query %v row %d: got %q want %q", qq, i, got[i], want[i])
			}
		}
	}
}

// TestCacheReuseAndEviction: repeated point queries in one window hydrate once and
// then reuse (one cache entry); querying more distinct windows than the cap evicts
// the oldest, and an evicted window re-hydrates correctly on the next query.
func TestCacheReuseAndEviction(t *testing.T) {
	path := writeBigBED(t)
	r, err := NewReaderSize(path, 2) // cap = 2 windows
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()

	// Two point queries in the same 16 kb window → one entry, served from cache.
	queryLines(t, r, "chr1", 100, 101)
	queryLines(t, r, "chr1", 200, 201)
	if n := r.recCache.ll.Len(); n != 1 {
		t.Fatalf("after two same-window queries, cache entries = %d, want 1", n)
	}

	// Touch three distinct windows → cap 2 must evict the oldest.
	queryLines(t, r, "chr1", 20000, 20001) // window 1
	queryLines(t, r, "chr1", 40000, 40001) // window 2 → evicts window 0
	if n := r.recCache.ll.Len(); n != 2 {
		t.Fatalf("after three windows with cap 2, entries = %d, want 2", n)
	}

	// The evicted window 0 re-hydrates and still returns the right record.
	got := queryLines(t, r, "chr1", 100, 101)
	if len(got) != 1 || got[0] != "chr1\t100\t120\tchr1_feat2" {
		t.Fatalf("re-query evicted window = %v, want [chr1_feat2]", got)
	}
}
