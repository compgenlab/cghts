package gtstore

import (
	"context"
	"encoding/json"
	"errors"
	"os"
	"path/filepath"
	"testing"
)

// Data from before the gtstore rename is refused by name, at every level a
// reader can meet it: a volume, a store and a fragment. Reading it instead
// would look up gtstore.* metadata keys in a file that carries cgkit.* ones,
// and absent metadata is a wrong answer, not an error.

func TestOpenRefusesAPreRenameVolume(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	buildCensusStore(t, base, WriterOpts{MinDP: 10})
	m, err := ReadVolumeManifest(base)
	if err != nil {
		t.Fatal(err)
	}
	m.FormatVersion = 1
	if err := WriteVolumeManifest(base, *m); err != nil {
		t.Fatal(err)
	}
	if _, err := OpenParquet(base); !errors.Is(err, ErrLegacyStore) {
		t.Fatalf("OpenParquet = %v, want ErrLegacyStore", err)
	}
}

func TestOpenRefusesAPreRenameStoreManifest(t *testing.T) {
	dir := buildSet(t, []string{"chr1", "chr2"}, []string{"S1", "S2"})
	man, err := readStoreManifest(context.Background(), dir)
	if err != nil {
		t.Fatal(err)
	}
	man.FormatVersion = 1
	if err := WriteStoreManifest(dir, *man); err != nil {
		t.Fatal(err)
	}
	if _, err := OpenStore(context.Background(), dir); !errors.Is(err, ErrLegacyStore) {
		t.Fatalf("OpenStore = %v, want ErrLegacyStore", err)
	}
}

// The old archive marker is the only thing at the top of a pre-rename store,
// so without recognising it the error would be "cannot tell what kind of
// store", which sends the reader off to check their path.
func TestOpenNamesAPreRenameArchiveMarker(t *testing.T) {
	dir := buildSet(t, []string{"chr1", "chr2"}, []string{"S1", "S2"})
	if err := os.Rename(StoreManifestPath(dir), filepath.Join(dir, LegacyStoreMarker)); err != nil {
		t.Fatal(err)
	}
	ctx := context.Background()
	if _, err := Open(ctx, dir, ""); !errors.Is(err, ErrLegacyStore) {
		t.Fatalf("Open = %v, want ErrLegacyStore", err)
	}
	if _, err := OpenStore(ctx, dir); !errors.Is(err, ErrLegacyStore) {
		t.Fatalf("OpenStore = %v, want ErrLegacyStore", err)
	}
	if _, err := StoreKind(ctx, dir); !errors.Is(err, ErrLegacyStore) {
		t.Fatalf("StoreKind = %v, want ErrLegacyStore", err)
	}
}

func TestReadFragmentMetaRefusesAPreRenameFragment(t *testing.T) {
	dir := buildFragment(t, "S1", RefBlock{SampleID: "S1", Chrom: "chr1", Start: 100, End: 200, MinDP: 20})
	p := filepath.Join(dir, FragmentMetaFile)
	raw, err := os.ReadFile(p)
	if err != nil {
		t.Fatal(err)
	}
	var m map[string]any
	if err := json.Unmarshal(raw, &m); err != nil {
		t.Fatal(err)
	}
	m["format_version"] = 1
	raw, _ = json.Marshal(m)
	if err := os.WriteFile(p, raw, 0o644); err != nil {
		t.Fatal(err)
	}
	if _, err := ReadFragmentMeta(context.Background(), dir); !errors.Is(err, ErrLegacyStore) {
		t.Fatalf("ReadFragmentMeta = %v, want ErrLegacyStore", err)
	}
}
