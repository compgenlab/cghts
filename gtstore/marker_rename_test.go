package gtstore

import (
	"context"
	"os"
	"path/filepath"
	"testing"
)

// An archive bound before the gtset rename carries its manifest as
// gtstore.json.gz. Only the file name changed, so it must still open, and the
// next write moves it to the current name.
func TestArchiveUnderThePreviousMarkerNameStillOpens(t *testing.T) {
	ctx := context.Background()
	dir := buildSet(t, []string{"chr1", "chr2"}, []string{"S1", "S2"})
	if _, err := os.Stat(filepath.Join(dir, StoreManifestFile)); err != nil {
		t.Fatalf("BuildStore did not write %s: %v", StoreManifestFile, err)
	}
	prev := filepath.Join(dir, PreviousStoreManifestFile)
	if err := os.Rename(StoreManifestPath(dir), prev); err != nil {
		t.Fatal(err)
	}

	if !IsStore(ctx, dir) {
		t.Fatal("an archive under the previous marker name is not recognised")
	}
	if _, err := Open(ctx, dir, ""); err != nil {
		t.Fatalf("Open: %v", err)
	}
	st, err := OpenStore(ctx, dir)
	if err != nil {
		t.Fatalf("OpenStore: %v", err)
	}
	if s, _ := st.Samples(); len(s) != 2 {
		t.Fatalf("samples = %v", s)
	}

	// Writing the manifest again moves it to the current name.
	man, err := readStoreManifest(ctx, dir)
	if err != nil {
		t.Fatal(err)
	}
	if err := WriteStoreManifest(dir, *man); err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(StoreManifestPath(dir)); err != nil {
		t.Errorf("current marker not written: %v", err)
	}
	if _, err := os.Stat(prev); !os.IsNotExist(err) {
		t.Errorf("previous marker not removed: %v", err)
	}
}

// When both names are present the current one is the archive's manifest.
func TestCurrentMarkerWinsOverThePreviousName(t *testing.T) {
	ctx := context.Background()
	dir := buildSet(t, []string{"chr1", "chr2"}, []string{"S1", "S2"})
	man, err := readStoreManifest(ctx, dir)
	if err != nil {
		t.Fatal(err)
	}
	stale := *man
	stale.Program = "stale"
	// Write a stale copy under the previous name by hand.
	if err := WriteStoreManifest(dir, stale); err != nil {
		t.Fatal(err)
	}
	if err := os.Rename(StoreManifestPath(dir), filepath.Join(dir, PreviousStoreManifestFile)); err != nil {
		t.Fatal(err)
	}
	if err := WriteStoreManifest(dir, *man); err != nil {
		t.Fatal(err)
	}
	// WriteStoreManifest removed the stale copy; put one back to test precedence.
	if err := os.WriteFile(filepath.Join(dir, PreviousStoreManifestFile), []byte("not gzip"), 0o644); err != nil {
		t.Fatal(err)
	}
	got, err := readStoreManifest(ctx, dir)
	if err != nil {
		t.Fatalf("read with both names present: %v", err)
	}
	if got.Program == "stale" {
		t.Error("the previous name was read although the current one exists")
	}
}
