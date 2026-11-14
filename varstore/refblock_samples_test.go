package varstore

import (
	"context"
	"path/filepath"
	"reflect"
	"testing"
)

func attachedList(t *testing.T, base string) ([]string, bool) {
	t.Helper()
	s, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	return s.RefBlockSamples()
}

// The list records whose fragments were attached, not who has blocks: S2's
// fragment has nothing on chr1, so the table holds no S2 row, but S2's gVCF WAS
// read -- which is what separates "not covered" from "never attached" for S3.
func TestAttachRecordsWhoseFragmentsWereAttached(t *testing.T) {
	ctx := context.Background()
	base := buildAttachVolume(t)

	if _, ok := attachedList(t, base); ok {
		t.Fatal("a volume never attached claims a list of attached samples")
	}

	if _, err := AttachRefBlocks(ctx, base, loadFrags(t,
		buildFragment(t, "S2", 10, blk("chr2", 1, 500, 30)),
		buildFragment(t, "S1", 10, blk("chr1", 900, 2500, 30)),
	), AttachOpts{}); err != nil {
		t.Fatal(err)
	}
	got, ok := attachedList(t, base)
	if !ok || !reflect.DeepEqual(got, []string{"S1", "S2"}) {
		t.Errorf("after attaching S1,S2: %v (known %v), want [S1 S2] in roster order", got, ok)
	}

	// Replace adds a sample whose fragment is new, keeping the others.
	if _, err := AttachRefBlocks(ctx, base, loadFrags(t,
		buildFragment(t, "S3", 10, blk("chr1", 1500, 1600, 40)),
	), AttachOpts{Replace: []string{"S3"}}); err != nil {
		t.Fatal(err)
	}
	if got, _ := attachedList(t, base); !reflect.DeepEqual(got, []string{"S1", "S2", "S3"}) {
		t.Errorf("after replacing S3 in: %v", got)
	}

	// Replace with no fragment removes the sample from the list as from the table.
	if _, err := AttachRefBlocks(ctx, base, nil, AttachOpts{Replace: []string{"S1"}}); err != nil {
		t.Fatal(err)
	}
	if got, _ := attachedList(t, base); !reflect.DeepEqual(got, []string{"S2", "S3"}) {
		t.Errorf("after removing S1: %v", got)
	}

	// Removing everyone leaves an empty list, which is a claim, not "unknown".
	if _, err := AttachRefBlocks(ctx, base, nil, AttachOpts{Replace: []string{"S2", "S3"}}); err != nil {
		t.Fatal(err)
	}
	if got, ok := attachedList(t, base); !ok || len(got) != 0 {
		t.Errorf("after removing everyone: %v (known %v), want a known empty list", got, ok)
	}

	// A rebuild states exactly its fragments.
	if _, err := AttachRefBlocks(ctx, base, loadFrags(t,
		buildFragment(t, "S3", 10, blk("chr1", 1500, 1600, 40)),
	), AttachOpts{Rebuild: true}); err != nil {
		t.Fatal(err)
	}
	if got, _ := attachedList(t, base); !reflect.DeepEqual(got, []string{"S3"}) {
		t.Errorf("after rebuilding from S3: %v", got)
	}
}

// A table attached before the list existed has no list, and replacing one
// sample in it must not invent one: the replaced sample alone would claim
// everybody else had been removed.
func TestReplaceOverAnUnrecordedTableStaysUnknown(t *testing.T) {
	ctx := context.Background()
	base := buildAttachVolume(t)
	if _, err := AttachRefBlocks(ctx, base, loadFrags(t,
		buildFragment(t, "S1", 10, blk("chr1", 900, 2500, 30)),
	), AttachOpts{}); err != nil {
		t.Fatal(err)
	}
	man, err := ReadVolumeManifestContext(ctx, base)
	if err != nil {
		t.Fatal(err)
	}
	man.RefBlockSamples = nil // as an attach before this field would have written
	if err := WriteVolumeManifest(base, *man); err != nil {
		t.Fatal(err)
	}

	if _, err := AttachRefBlocks(ctx, base, loadFrags(t,
		buildFragment(t, "S2", 10, blk("chr1", 1500, 1600, 40)),
	), AttachOpts{Replace: []string{"S2"}}); err != nil {
		t.Fatal(err)
	}
	if got, ok := attachedList(t, base); ok {
		t.Errorf("replace over an unrecorded table recorded %v", got)
	}
}

// An archive states the list once when its volumes agree, and per volume when
// they do not, so a disagreement is visible rather than averaged away.
func TestArchiveHoistsAnAgreedListAndKeepsADisagreement(t *testing.T) {
	ctx := context.Background()
	base := filepath.Join(t.TempDir(), "cohort")
	buildVolumeOn(t, filepath.Join(base, "chr1"), "chr1")
	buildVolumeOn(t, filepath.Join(base, "chr2"), "chr2")
	man, err := BuildStore(ctx, base, []string{"chr1", "chr2"}, "test", "")
	if err != nil {
		t.Fatal(err)
	}
	if err := WriteStoreManifest(base, *man); err != nil {
		t.Fatal(err)
	}
	if _, err := AttachRefBlocksToStore(ctx, base, loadFrags(t,
		buildFragment(t, "S1", 10, blk("chr1", 1, 900, 30)),
	), AttachOpts{}, "test", ""); err != nil {
		t.Fatal(err)
	}

	s, err := OpenStore(ctx, base)
	if err != nil {
		t.Fatal(err)
	}
	got, ok := s.RefBlockSamples()
	sm := s.StoreManifest()
	s.Close()
	if !ok || !reflect.DeepEqual(got, []string{"S1"}) {
		t.Errorf("archive list %v (known %v), want [S1]", got, ok)
	}
	for _, v := range sm.Volumes {
		if v.RefBlockSamples != nil {
			t.Errorf("volume %s repeats the agreed list", v.Name)
		}
	}

	// Replace on one volume only, as a run dying partway would leave it.
	if _, err := AttachRefBlocks(ctx, filepath.Join(base, "chr2"), loadFrags(t,
		buildFragment(t, "S2", 10, blk("chr2", 1, 900, 30)),
	), AttachOpts{Replace: []string{"S2"}}); err != nil {
		t.Fatal(err)
	}
	man, err = BuildStore(ctx, base, []string{"chr1", "chr2"}, "test", "")
	if err != nil {
		t.Fatal(err)
	}
	if man.RefBlockSamples != nil {
		t.Errorf("disagreeing volumes hoisted %v", *man.RefBlockSamples)
	}
	want := map[string][]string{"chr1": {"S1"}, "chr2": {"S1", "S2"}}
	for _, v := range man.Volumes {
		if v.RefBlockSamples == nil || !reflect.DeepEqual(*v.RefBlockSamples, want[v.Name]) {
			t.Errorf("volume %s: %v, want %v", v.Name, v.RefBlockSamples, want[v.Name])
		}
	}
}
