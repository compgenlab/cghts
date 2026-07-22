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

// TestWriterLinearIndexCoversFullSpan verifies that the linear index records
// the file offset of a large feature for every 16kb window it spans, not just
// the window containing its start. Without this, a query at a distant intronic
// position can miss the gene row because the linear-index guard filters out
// the chunk that contains it.
func TestWriterLinearIndexCoversFullSpan(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "genes.gtf.gz")

	w := NewWriter(path, NewWriterOpts().GFF().AutoIndex())

	// A small record before the gene so the gene's file offset is non-zero
	// (distinguishable from "unset" 0 entries in the linear index).
	if err := w.Write("chr22\tRefSeq\texon\t10000\t10100\t.\t+\t.\tgene_id \"EARLY\"; transcript_id \"NM_000\";"); err != nil {
		t.Fatal(err)
	}
	// A gene row spanning ~200kb (like CECR2), covering many 16kb windows.
	geneLine := "chr22\tRefSeq\tgene\t17360000\t17558000\t.\t+\t.\tgene_id \"TESTGENE\";"
	if err := w.Write(geneLine); err != nil {
		t.Fatal(err)
	}
	// A small exon far from the gene start, in a different 16kb window.
	if err := w.Write("chr22\tRefSeq\texon\t17500000\t17500100\t.\t+\t.\tgene_id \"TESTGENE\"; transcript_id \"NM_001\";"); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	// Read the index and check the linear index entries.
	r, err := NewReader(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()

	idx, ok := r.idx.(*BinIndex)
	if !ok {
		t.Skip("index is not a BinIndex (CSI), skipping linear index check")
	}

	refID := idx.RefID("chr22")
	if refID < 0 {
		t.Fatal("chr22 not found in index")
	}
	ref := &idx.refs[refID]

	// The gene spans 0-based [17359999, 17558000), covering 16kb windows
	// 17359999>>14 = 1059 through (17558000-1)>>14 = 1071.
	geneStartWin := 17359999 >> 14 // 1059
	_ = (17558000 - 1) >> 14 // geneEndWin = 1071

	// Every window the gene spans must have a non-zero linear index entry
	// (or the same offset as the gene start window). A zero entry in windows
	// past the gene start means the linear index wasn't updated for those
	// windows, so the gene row's chunk could be filtered out during queries.
	//
	// The gene row is the first record, so its offset (at the gene start
	// window) may be 0. But windows between the gene start and the exon
	// (17500000, window 1068) that have no records should still have the
	// gene row's offset propagated, not 0 (which means "unset" for non-first
	// windows with no records starting in them).
	geneRowOffset := ref.linearIdx[geneStartWin]

	// The exon is at 0-based 17499999, window 17499999>>14 = 1068.
	// Check windows between gene start+1 and exon-1 (1060..1067) — these
	// have no small features starting in them.
	exonWin := 17499999 >> 14 // 1068
	for w := geneStartWin + 1; w < exonWin; w++ {
		if w >= len(ref.linearIdx) {
			t.Errorf("linear index too short: window %d beyond length %d", w, len(ref.linearIdx))
			continue
		}
		vo := ref.linearIdx[w]
		// Without the fix, these windows have 0 (unset) because only the
		// gene start window was updated. With the fix, they should have
		// the gene row's offset.
		if vo != geneRowOffset {
			t.Errorf("window %d: linearIdx=%d, want %d (gene row offset); gene row spans this window but linear index wasn't updated",
				w, vo, geneRowOffset)
		}
	}
}
