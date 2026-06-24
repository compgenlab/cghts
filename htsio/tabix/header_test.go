package tabix

import (
	"path/filepath"
	"testing"
)

func TestColumnNamesFromSkippedHeader(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "scores.tab.gz")
	w := NewWriter(path, NewWriterOpts().Columns(1, 2, 0).Meta('#').Skip(1).AutoIndex())
	w.WriteHeader("#chrom\tpos\tref\talt\tscore")
	for _, l := range []string{
		"chr1\t100\tA\tG\t0.9",
		"chr2\t500\tC\tA\t0.3",
	} {
		if err := w.Write(l); err != nil {
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

	names, err := r.ColumnNames()
	if err != nil {
		t.Fatal(err)
	}
	want := []string{"chrom", "pos", "ref", "alt", "score"}
	if len(names) != len(want) {
		t.Fatalf("ColumnNames = %v, want %v", names, want)
	}
	for i := range want {
		if names[i] != want[i] {
			t.Errorf("ColumnNames[%d] = %q, want %q", i, names[i], want[i])
		}
	}
	// 1-based lookups.
	for name, col := range map[string]int{"chrom": 1, "alt": 4, "score": 5} {
		if c, err := r.ColumnByName(name); err != nil || c != col {
			t.Errorf("ColumnByName(%q) = %d (err=%v), want %d", name, c, err, col)
		}
	}
	if _, err := r.ColumnByName("missing"); err == nil {
		t.Error("expected error for missing column")
	}
	// A data row is still queryable (header is skipped).
	if got := queryLines(t, r, "chr2", 499, 500); len(got) != 1 {
		t.Errorf("chr2 query = %v, want 1 row", got)
	}
}

func TestColumnNamesMetaCommentButNoSkip(t *testing.T) {
	// A '#'-comment header line with Skip=0 is NOT a column header: without a
	// skipped line there is no header, even though the comment line exists.
	dir := t.TempDir()
	path := filepath.Join(dir, "meta.tab.gz")
	w := NewWriter(path, NewWriterOpts().Columns(1, 2, 0).Meta('#').AutoIndex()) // Skip defaults to 0
	w.WriteHeader("#chrom\tpos\tref\talt\tscore")
	if err := w.Write("chr1\t100\tA\tG\t0.9"); err != nil {
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

	names, err := r.ColumnNames()
	if err != nil {
		t.Fatal(err)
	}
	if len(names) != 0 {
		t.Errorf("ColumnNames = %v, want none (Skip=0 => no header)", names)
	}
	if _, err := r.ColumnByName("alt"); err == nil {
		t.Error("expected error: a meta comment without a skipped line is not a header")
	}
	// The comment line must not break data queries.
	if got := queryLines(t, r, "chr1", 99, 100); len(got) != 1 {
		t.Errorf("query rows = %v, want 1", got)
	}
}

func TestColumnNamesNoHeader(t *testing.T) {
	// No skipped line => no column header, even with a meta character.
	dir := t.TempDir()
	path := filepath.Join(dir, "noheader.bed.gz")
	w := NewWriter(path, NewWriterOpts().BED().AutoIndex())
	if err := w.Write("chr1\t90\t110\tgeneA"); err != nil {
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
	names, err := r.ColumnNames()
	if err != nil {
		t.Fatal(err)
	}
	if names != nil {
		t.Errorf("ColumnNames = %v, want nil (no header)", names)
	}
	if _, err := r.ColumnByName("anything"); err == nil {
		t.Error("expected error resolving a name with no header")
	}
}
