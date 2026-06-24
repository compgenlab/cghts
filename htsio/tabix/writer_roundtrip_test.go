package tabix

import (
	"path/filepath"
	"testing"
)

// queryLines returns the raw lines from a region query.
func queryLines(t *testing.T, r *Reader, ref string, start, end int) []string {
	t.Helper()
	seq, err := r.Query(ref, start, end)
	if err != nil {
		t.Fatalf("Query(%s,%d,%d): %v", ref, start, end, err)
	}
	var out []string
	for rec, err := range seq {
		if err != nil {
			t.Fatalf("query iter: %v", err)
		}
		out = append(out, rec.Line)
	}
	return out
}

// TestWriterLastRecordPerRefQueryable guards against a regression where the
// index closed each bin's chunk at the *start* of the final record, producing a
// zero-length chunk that made the last record of a reference (and any reference
// whose only record was last) unqueryable.
func TestWriterLastRecordPerRefQueryable(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "regions.bed.gz")

	w := NewWriter(path, NewWriterOpts().BED().AutoIndex())
	// chr2 has a single record and is the last reference written.
	for _, line := range []string{
		"chr1\t90\t110\tgeneA",
		"chr1\t145\t155\tenhB",
		"chr2\t400\t600\tgeneC",
	} {
		if err := w.Write(line); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	r, err := NewReader(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()

	// The single/last chr2 record must be found.
	got := queryLines(t, r, "chr2", 499, 500)
	if len(got) != 1 || got[0] != "chr2\t400\t600\tgeneC" {
		t.Errorf("chr2 query = %v, want [geneC]", got)
	}
	// chr1 records still queryable.
	if got := queryLines(t, r, "chr1", 99, 100); len(got) != 1 || got[0] != "chr1\t90\t110\tgeneA" {
		t.Errorf("chr1:100 query = %v, want [geneA]", got)
	}
	if got := queryLines(t, r, "chr1", 149, 150); len(got) != 1 || got[0] != "chr1\t145\t155\tenhB" {
		t.Errorf("chr1:150 query = %v, want [enhB]", got)
	}
}

// TestWriterSingleRecordQueryable covers the simplest case: one record total.
func TestWriterSingleRecordQueryable(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "one.bed.gz")
	w := NewWriter(path, NewWriterOpts().BED().AutoIndex())
	if err := w.Write("chr1\t1000\t2000\tonly"); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	r, err := NewReader(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	if got := queryLines(t, r, "chr1", 1500, 1501); len(got) != 1 {
		t.Errorf("single-record query = %v, want 1 row", got)
	}
}
