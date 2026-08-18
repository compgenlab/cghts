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
// the tables rename, so a test can produce a store nobody can write any more.
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

// A store written before the tables rename is REFUSED, and told why.
//
// It used to be read through a compatibility path, which is gone: one consumer,
// cheap re-conversion, and a second format to keep correct forever is a poor
// trade for both.
//
// REFUSING IS THE POINT, and not merely the consequence. Without the old key
// the table index comes back empty, and an empty index does not fail on its
// own -- it makes every row-count check pass over nothing. The store would
// open, answer, and have silently skipped the one verification that catches a
// table belonging to another conversion. So the absence is caught where it can
// still be explained, rather than surfacing as a row-count mismatch or as
// nothing at all.
func TestAPreRenameStoreIsRefusedWithAReason(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	buildCensusStore(t, base, WriterOpts{MinDP: 10})
	renameTablesKeyToMembers(t, base)

	_, err := OpenParquet(base)
	if err == nil {
		t.Fatal("a store predating the tables rename opened, so its tables were never verified")
	}
	// The remedy has to be in the message: there is nothing to retry and no
	// flag to pass, only a re-conversion.
	for _, want := range []string{"re-convert", "rename"} {
		if !strings.Contains(err.Error(), want) {
			t.Errorf("the refusal does not mention %q: %v", want, err)
		}
	}
}
