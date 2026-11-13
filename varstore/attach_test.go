package varstore

import (
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"reflect"
	"sort"
	"testing"
)

// A sharded volume (one site per shard) with catalog sites at chr1:1000 and
// chr1:3000, samples S1-S3, no reference blocks yet. The catalog's second shard
// starts at 3000, so attached blocks are cut there.
func buildAttachVolume(t *testing.T) string {
	t.Helper()
	base := filepath.Join(t.TempDir(), "chr1")
	w, err := NewWriter(base, WriterOpts{
		Samples: []string{"S1", "S2", "S3"}, MinDP: 10, RowGroupSize: 4, ShardSites: 1,
	})
	if err != nil {
		t.Fatal(err)
	}
	for _, pos := range []int32{1000, 3000} {
		if err := w.WriteSite(Site{Chrom: "chr1", Pos: pos, Ref: "G", Alt: "A", AC: 1, AN: 6, NCarriers: 1, NCalled: 3}); err != nil {
			t.Fatal(err)
		}
		if err := w.WriteCall(Call{SampleID: "S3", Chrom: "chr1", Pos: pos, Ref: "G", Alt: "A",
			GT: "0/1", DP: 40, ADRef: 20, ADAlt: 20, GQ: 99}); err != nil {
			t.Fatal(err)
		}
		for _, s := range []string{"S1", "S2", "S3"} {
			if err := w.WriteRegion(CalledSiteRun{SampleID: s, Chrom: "chr1", Start: pos, End: pos, NSites: 1, MinDP: 30}); err != nil {
				t.Fatal(err)
			}
		}
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
	return base
}

func buildFragment(t *testing.T, sample string, maxGap int32, blocks ...RefBlock) string {
	t.Helper()
	dir := filepath.Join(t.TempDir(), sample)
	w, err := NewFragmentWriter(dir, FragmentMeta{SampleID: sample, MinDP: 10, MaxGap: maxGap})
	if err != nil {
		t.Fatal(err)
	}
	for _, b := range blocks {
		b.SampleID = sample
		if err := w.Write(b); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
	return dir
}

func loadFrags(t *testing.T, dirs ...string) []Fragment {
	t.Helper()
	f, err := LoadFragments(context.Background(), dirs)
	if err != nil {
		t.Fatal(err)
	}
	return f
}

func statesAt(t *testing.T, base string, pos int32) map[string]rbWant {
	t.Helper()
	s, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	return classifyBoth(t, s, Locus{Chrom: "chr1", Pos: pos, Ref: "C", Alt: "T"}, Gate{MinDP: 10})
}

func refBlockRows(t *testing.T, base string) []string {
	t.Helper()
	s, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	var out []string
	if err := s.RefBlocks(func(b RefBlock) bool {
		out = append(out, fmt.Sprintf("%s %s:%d-%d @%d", b.SampleID, b.Chrom, b.Start, b.End, b.MinDP))
		return true
	}); err != nil {
		t.Fatal(err)
	}
	sort.Strings(out)
	return out
}

func TestAttachBuildsRefBlocksFromFragments(t *testing.T) {
	base := buildAttachVolume(t)
	frags := loadFrags(t,
		buildFragment(t, "S1", 10, blk("chr1", 900, 2500, 30), blk("chr2", 100, 200, 30)),
		buildFragment(t, "S2", 10, blk("chr1", 1500, 3100, 40)),
	)
	// A batch of 2 over 2 inputs is one level; force the tree anyway with
	// three inputs by adding an empty-on-chr1 fragment.
	frags = append(frags, loadFrags(t, buildFragment(t, "S3", 10, blk("chr2", 5, 9, 50)))...)
	res, err := AttachRefBlocks(context.Background(), base, frags, AttachOpts{Batch: 2})
	if err != nil {
		t.Fatal(err)
	}

	// S2's block crosses the cut at 3000 and is split; chr2 is not this volume's.
	want := []string{
		"S1 chr1:900-2500 @30",
		"S2 chr1:1500-2999 @40",
		"S2 chr1:3000-3100 @40",
	}
	if got := refBlockRows(t, base); !reflect.DeepEqual(got, want) {
		t.Errorf("ref_blocks:\n got %v\nwant %v", got, want)
	}
	if res.Rows != 3 || res.Shards != 2 || res.Samples != 2 {
		t.Errorf("result = %+v, want 3 rows in 2 shards over 2 samples", res)
	}
	man, err := ReadVolumeManifest(base)
	if err != nil {
		t.Fatal(err)
	}
	if !man.Params.RefBlocks || man.Params.MaxGap != 10 || man.Tables[RefBlocksTable].Rows != 3 {
		t.Errorf("manifest params %+v, ref_blocks %+v", man.Params, man.Tables[RefBlocksTable])
	}

	nc := rbWant{StateNonCarrier, RefFromBlock}
	na := rbWant{StateNotAssayed, ""}
	if got := statesAt(t, base, 2000); !reflect.DeepEqual(got, map[string]rbWant{"S1": nc, "S2": nc, "S3": na}) {
		t.Errorf("off-catalog 2000: %v", got)
	}
	// The split-off piece answers past the cut.
	if got := statesAt(t, base, 3050); got["S2"] != nc {
		t.Errorf("off-catalog 3050: S2 = %v, want covered by its split block", got["S2"])
	}
}

func TestAttachReplacesOneSampleAndKeepsTheRest(t *testing.T) {
	base := buildAttachVolume(t)
	if _, err := AttachRefBlocks(context.Background(), base, loadFrags(t,
		buildFragment(t, "S1", 10, blk("chr1", 900, 2500, 30)),
		buildFragment(t, "S2", 10, blk("chr1", 1500, 3100, 40)),
	), AttachOpts{}); err != nil {
		t.Fatal(err)
	}
	before, _ := filepath.Glob(filepath.Join(base, RefBlocksTable, "*.parquet"))

	// S2's gVCF was redone: its new blocks stop well short of 2000.
	_, err := AttachRefBlocks(context.Background(), base, loadFrags(t,
		buildFragment(t, "S2", 10, blk("chr1", 1500, 1600, 40)),
	), AttachOpts{Replace: []string{"S2"}})
	if err != nil {
		t.Fatal(err)
	}
	want := []string{"S1 chr1:900-2500 @30", "S2 chr1:1500-1600 @40"}
	if got := refBlockRows(t, base); !reflect.DeepEqual(got, want) {
		t.Errorf("after replacing S2:\n got %v\nwant %v", got, want)
	}
	if got := statesAt(t, base, 2000); got["S1"].state != StateNonCarrier || got["S2"].state != StateNotAssayed {
		t.Errorf("off-catalog 2000 after replace: %v", got)
	}
	// The previous generation is gone; only the new one is referenced.
	for _, f := range before {
		if _, err := os.Stat(f); err == nil {
			t.Errorf("%s from the replaced generation is still on disk", f)
		}
	}

	// Replacing with no fragment removes the sample's blocks.
	if _, err := AttachRefBlocks(context.Background(), base, nil, AttachOpts{Replace: []string{"S1"}}); err != nil {
		t.Fatal(err)
	}
	if got := refBlockRows(t, base); !reflect.DeepEqual(got, []string{"S2 chr1:1500-1600 @40"}) {
		t.Errorf("after removing S1: %v", got)
	}
}

func TestAttachRefusals(t *testing.T) {
	ctx := context.Background()
	base := buildAttachVolume(t)
	s1 := loadFrags(t, buildFragment(t, "S1", 10, blk("chr1", 900, 2500, 30)))

	if _, err := AttachRefBlocks(ctx, base, s1, AttachOpts{Replace: []string{"S1"}}); err == nil {
		t.Error("replace on a volume with no ref_blocks was accepted")
	}
	if _, err := AttachRefBlocks(ctx, base, s1, AttachOpts{}); err != nil {
		t.Fatal(err)
	}
	if _, err := AttachRefBlocks(ctx, base, s1, AttachOpts{}); !errors.Is(err, ErrRefBlocksExist) {
		t.Errorf("a second attach without rebuild or replace: err = %v", err)
	}
	if _, err := AttachRefBlocks(ctx, base, s1, AttachOpts{Rebuild: true, Replace: []string{"S1"}}); err == nil {
		t.Error("rebuild and replace together were accepted")
	}
	s2 := loadFrags(t, buildFragment(t, "S2", 10, blk("chr1", 1, 5, 30)))
	if _, err := AttachRefBlocks(ctx, base, s2, AttachOpts{Replace: []string{"S1"}}); err == nil {
		t.Error("a fragment for a sample not being replaced was accepted")
	}
	if _, err := AttachRefBlocks(ctx, base, s1, AttachOpts{Rebuild: true}); err != nil {
		t.Errorf("rebuild: %v", err)
	}

	off := loadFrags(t, buildFragment(t, "S9", 10, blk("chr1", 1, 5, 30)))
	if _, err := AttachRefBlocks(ctx, base, off, AttachOpts{Rebuild: true}); err == nil {
		t.Error("a fragment for a sample off the roster was accepted")
	}

	dir := filepath.Join(t.TempDir(), "gated")
	w, _ := NewFragmentWriter(dir, FragmentMeta{SampleID: "S1", MinDP: 20, MaxGap: 10})
	_ = w.Finish()
	if _, err := AttachRefBlocks(ctx, base, loadFrags(t, dir), AttachOpts{Rebuild: true}); err == nil {
		t.Error("a fragment gated at a different --min-dp was accepted")
	}

	if _, err := LoadFragments(ctx, []string{
		buildFragment(t, "S1", 10), buildFragment(t, "S1", 10),
	}); err == nil {
		t.Error("two fragments for one sample were accepted")
	}
	if _, err := LoadFragments(ctx, []string{
		buildFragment(t, "S1", 10), buildFragment(t, "S2", 50),
	}); err == nil {
		t.Error("fragments merged at different gap tolerances were accepted")
	}
}

// The default volume is unsharded: the whole chromosome is one output file.
func TestAttachToAnUnshardedVolume(t *testing.T) {
	base := filepath.Join(t.TempDir(), "chr1")
	w, err := NewWriter(base, WriterOpts{Samples: []string{"S1"}, MinDP: 10})
	if err != nil {
		t.Fatal(err)
	}
	if err := w.WriteSite(Site{Chrom: "chr1", Pos: 1000, Ref: "G", Alt: "A", AN: 2, NCalled: 1}); err != nil {
		t.Fatal(err)
	}
	if err := w.WriteRegion(CalledSiteRun{SampleID: "S1", Chrom: "chr1", Start: 1000, End: 1000, NSites: 1, MinDP: 30}); err != nil {
		t.Fatal(err)
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
	res, err := AttachRefBlocks(context.Background(), base, loadFrags(t,
		buildFragment(t, "S1", 10, blk("chr1", 1, 999, 30), blk("chr1", 1001, 90000, 30)),
	), AttachOpts{})
	if err != nil {
		t.Fatal(err)
	}
	if res.Shards != 1 || res.Rows != 2 {
		t.Errorf("result = %+v, want 2 rows in one file", res)
	}
	if got := statesAt(t, base, 50000); got["S1"] != (rbWant{StateNonCarrier, RefFromBlock}) {
		t.Errorf("off-catalog 50000: %v", got["S1"])
	}
}

func buildVolumeOn(t *testing.T, base, chrom string) {
	t.Helper()
	w, err := NewWriter(base, WriterOpts{Samples: []string{"S1", "S2"}, MinDP: 10})
	if err != nil {
		t.Fatal(err)
	}
	if err := w.WriteSite(Site{Chrom: chrom, Pos: 1000, Ref: "G", Alt: "A", AN: 4, NCalled: 2}); err != nil {
		t.Fatal(err)
	}
	for _, s := range []string{"S1", "S2"} {
		if err := w.WriteRegion(CalledSiteRun{SampleID: s, Chrom: chrom, Start: 1000, End: 1000, NSites: 1, MinDP: 30}); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
}

func classifyStore(t *testing.T, s Store, chrom string, pos int32) map[string]State {
	t.Helper()
	states, err := s.Classify(Locus{Chrom: chrom, Pos: pos, Ref: "C", Alt: "T"}, Gate{MinDP: 10})
	if err != nil {
		t.Fatal(err)
	}
	return statesByName(states)
}

// A chromosome with no catalog site belonged to no volume, so its blocks had
// nowhere to be read from and every query there was not assayed. On a single
// volume, the volume now takes the chromosome.
func TestAttachToAVolumeTakesChromosomesWithoutSites(t *testing.T) {
	base := filepath.Join(t.TempDir(), "vol")
	buildVolumeOn(t, base, "chr1")
	frags := loadFrags(t, buildFragment(t, "S1", 10, blk("chr1", 1, 900, 30), blk("chr2", 100, 200, 30)))
	if _, err := AttachRefBlocksToStore(context.Background(), base, frags, AttachOpts{}, "test", ""); err != nil {
		t.Fatal(err)
	}
	s, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	if got := classifyStore(t, s, "chr2", 150); got["S1"] != StateNonCarrier || got["S2"] != StateNotAssayed {
		t.Errorf("chr2:150 on a volume whose catalog is chr1 only: %v", got)
	}
}

// On an archive, a chromosome no volume holds gets a blocks-only volume of its
// own, registered in the archive, so a query there is routed and answered.
func TestAttachToAnArchiveCreatesBlocksOnlyVolumes(t *testing.T) {
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

	frags := loadFrags(t,
		buildFragment(t, "S1", 10, blk("chr1", 1, 900, 30), blk("chr2", 1, 900, 30), blk("chrM", 1, 16000, 200)),
		buildFragment(t, "S2", 10, blk("chrM", 1, 500, 150)),
	)
	res, err := AttachRefBlocksToStore(ctx, base, frags, AttachOpts{}, "test", "")
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(res.Created, []string{"chrM"}) {
		t.Errorf("created %v, want a chrM volume", res.Created)
	}

	s, err := OpenStore(ctx, base)
	if err != nil {
		t.Fatal(err)
	}
	if got := classifyStore(t, s, "chrM", 1000); got["S1"] != StateNonCarrier || got["S2"] != StateNotAssayed {
		t.Errorf("chrM:1000: %v", got)
	}
	if got := classifyStore(t, s, "chr2", 500); got["S1"] != StateNonCarrier {
		t.Errorf("chr2:500 off the catalog: %v", got)
	}
	if got := classifyStore(t, s, "chrX", 500); got["S1"] != StateNotAssayed {
		t.Errorf("chrX, which no fragment covers: %v", got)
	}
	s.Close()

	// A rebuild now finds chrM held, and creates nothing.
	res, err = AttachRefBlocksToStore(ctx, base, frags, AttachOpts{Rebuild: true}, "test", "")
	if err != nil {
		t.Fatal(err)
	}
	if len(res.Created) != 0 || len(res.Volumes) != 3 {
		t.Errorf("rebuild: created %v over %d volumes, want none over 3", res.Created, len(res.Volumes))
	}
}
