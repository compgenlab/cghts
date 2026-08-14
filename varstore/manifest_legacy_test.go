package varstore

import (
	"bytes"
	"compress/gzip"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// renameTablesKeyToMembers rewrites a manifest the way it was written before
// the tables/volumes rename, so a test can read a store nobody can produce any
// more.
func renameTablesKeyToMembers(t *testing.T, base string) {
	t.Helper()
	path := VolumeManifestPath(base)
	raw, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	zr, err := gzip.NewReader(bytes.NewReader(raw))
	if err != nil {
		t.Fatal(err)
	}
	doc, err := io.ReadAll(zr)
	if err != nil {
		t.Fatal(err)
	}
	old := bytes.Replace(doc, []byte(`"tables":`), []byte(`"members":`), 1)
	if bytes.Equal(old, doc) {
		t.Fatal(`the manifest carries no "tables" key, so this test rewrote nothing`)
	}
	var buf bytes.Buffer
	zw := gzip.NewWriter(&buf)
	if _, err := zw.Write(old); err != nil {
		t.Fatal(err)
	}
	if err := zw.Close(); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(path, buf.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}
}

// A store written before the rename still opens.
//
// The table index moved from "members" to "tables" when "member" stopped
// meaning two things at once. Every store converted before that carries the old
// key, and reconverting a whole-genome callset to rename a JSON field is not a
// migration anyone should be asked to run.
func TestOpenReadsThePreRenameTablesKey(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	buildCensusStore(t, base, WriterOpts{MinDP: 10})
	renameTablesKeyToMembers(t, base)

	man, err := ReadVolumeManifest(base)
	if err != nil {
		t.Fatal(err)
	}
	if len(man.Tables) == 0 {
		t.Fatal("the table index did not survive the old key")
	}
	if man.Tables[CallsTable].Rows == 0 {
		t.Errorf("calls table read with no rows: %+v", man.Tables)
	}

	s, err := OpenParquet(base)
	if err != nil {
		t.Fatalf("a pre-rename store no longer opens: %v", err)
	}
	defer s.Close()
	if _, err := s.Samples(); err != nil {
		t.Fatal(err)
	}
}

// And the check the index exists for still runs against it.
//
// THIS IS THE HALF THAT MATTERS. Reading the old key wrong does not fail
// loudly: an empty index means the row counts every table is checked against
// are absent, and a foreign table slipped into the store would be accepted in
// silence. Opening successfully proves nothing on its own -- this proves the
// index that came back is the one doing the work.
func TestThePreRenameTablesKeyStillCatchesAForeignTable(t *testing.T) {
	dir := t.TempDir()
	mine := filepath.Join(dir, "mine")
	theirs := filepath.Join(dir, "theirs")
	buildCensusStore(t, mine, WriterOpts{MinDP: 10})

	w, err := NewWriter(theirs, WriterOpts{Samples: []string{"S1", "S2"}, MinDP: 10, RowGroupSize: 4})
	if err != nil {
		t.Fatal(err)
	}
	if err := w.WriteSite(Site{Chrom: "chr1", Pos: 100, Ref: "A", Alt: "T", AC: 1, AN: 4}); err != nil {
		t.Fatal(err)
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}

	swapped, err := os.ReadFile(SitesPath(theirs))
	if err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(SitesPath(mine), swapped, 0o644); err != nil {
		t.Fatal(err)
	}
	renameTablesKeyToMembers(t, mine)

	_, err = OpenParquet(mine)
	if err == nil {
		t.Fatal("a pre-rename store accepted a sites table from a different conversion")
	}
	if !strings.Contains(err.Error(), "does not belong") {
		t.Errorf("error does not explain the mismatch: %v", err)
	}
}
