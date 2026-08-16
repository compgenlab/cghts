package varstore

import (
	"fmt"
	"path/filepath"
	"strings"
	"testing"
)

// writeCoverageStore builds a small store, optionally carrying coverage blocks.
func writeCoverageStore(t *testing.T, base string, opts WriterOpts, blocks []CoverageBlock) {
	t.Helper()
	if opts.Samples == nil {
		opts.Samples = []string{"S1", "S2"}
	}
	if opts.RowGroupSize == 0 {
		opts.RowGroupSize = 8
	}
	w, err := NewWriter(base, opts)
	if err != nil {
		t.Fatal(err)
	}
	for pos := int32(100); pos < 106; pos++ {
		if err := w.WriteSite(Site{
			Chrom: "chr1", Pos: pos, Ref: "A", Alt: "T", AC: 1, AN: 4, NCarriers: 1, NCalled: 2,
		}); err != nil {
			t.Fatal(err)
		}
		if err := w.WriteCall(Call{
			SampleID: "S1", Chrom: "chr1", Pos: pos, Ref: "A", Alt: "T", GT: "0/1", DP: 30,
		}); err != nil {
			t.Fatal(err)
		}
	}
	for _, s := range opts.Samples {
		if err := w.WriteRegion(CalledSiteRun{
			SampleID: s, Chrom: "chr1", Start: 100, End: 105, NSites: 6, MinDP: 30,
		}); err != nil {
			t.Fatal(err)
		}
	}
	for _, b := range blocks {
		if err := w.WriteCoverage(b); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
}

// Coverage blocks survive a round trip.
func TestCoverageBlocksRoundTrip(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	want := []CoverageBlock{
		{SampleID: "S1", Chrom: "chr1", Start: 1, End: 5000, MinDP: 31, GQ: 60},
		{SampleID: "S2", Chrom: "chr1", Start: 1, End: 4200, MinDP: 12, GQ: 40},
		{SampleID: "S2", Chrom: "chr1", Start: 4300, End: 9000, MinDP: 28, GQ: 55},
	}
	writeCoverageStore(t, base, WriterOpts{MinDP: 10, Coverage: true, MaxGap: 10}, want)

	s, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	if !s.HasCoverage() {
		t.Fatal("a store written with coverage reports none")
	}
	if s.MaxGap() != 10 {
		t.Errorf("MaxGap is %d, want 10; a reader cannot tell what covered meant here", s.MaxGap())
	}

	var got []CoverageBlock
	if err := s.Coverage(func(b CoverageBlock) bool { got = append(got, b); return true }); err != nil {
		t.Fatal(err)
	}
	if len(got) != len(want) {
		t.Fatalf("read %d coverage blocks, wrote %d", len(got), len(want))
	}
	// Sorted by (chrom, start), so compare as a set keyed by sample and start.
	byKey := map[string]CoverageBlock{}
	for _, b := range got {
		byKey[fmt.Sprintf("%s:%d", b.SampleID, b.Start)] = b
	}
	for _, w := range want {
		g, ok := byKey[fmt.Sprintf("%s:%d", w.SampleID, w.Start)]
		if !ok {
			t.Errorf("block %s %s:%d-%d did not come back", w.SampleID, w.Chrom, w.Start, w.End)
			continue
		}
		if g.End != w.End || g.MinDP != w.MinDP || g.GQ != w.GQ {
			t.Errorf("block %s:%d read back as end=%d min_dp=%d gq=%d, wrote end=%d min_dp=%d gq=%d",
				w.SampleID, w.Start, g.End, g.MinDP, g.GQ, w.End, w.MinDP, w.GQ)
		}
	}
}

// A store without coverage says so, rather than reporting none.
//
// THE DISTINCTION IS THE WHOLE POINT. "No coverage recorded" and "covered
// nowhere" are different claims, and only the second would be an answer. A
// store converted from a plain pVCF has made no claim off its catalog, and
// anything that reads it as a negative turns an unknown into a wrong count.
func TestAStoreWithoutCoverageIsNotAStoreCoveringNothing(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	writeCoverageStore(t, base, WriterOpts{MinDP: 10}, nil)

	s, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	if s.HasCoverage() {
		t.Error("a store given no coverage reports that it has some")
	}
	n := 0
	if err := s.Coverage(func(CoverageBlock) bool { n++; return true }); err != nil {
		t.Errorf("reading coverage from a store without it should be silent, got %v", err)
	}
	if n != 0 {
		t.Errorf("streamed %d blocks from a store that has none", n)
	}

	// And the manifest must not carry a zero-row entry, which would read as a
	// present-but-empty table to anything inspecting it.
	man, err := ReadVolumeManifest(base)
	if err != nil {
		t.Fatal(err)
	}
	if _, listed := man.Tables[CoverageTable]; listed {
		t.Error("the manifest lists a coverage table for a store that has none, " +
			"which reads as covered-nowhere rather than nothing-said")
	}
}

// Writing coverage to a store that was not opened for it is refused.
//
// Silently dropping the block would be the expensive failure: the gVCF it came
// from is not read again, so the span is gone and nothing says so until a
// question years later comes back NotAssayed.
func TestCoverageIsRefusedUnlessTheStoreWasOpenedForIt(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	w, err := NewWriter(base, WriterOpts{Samples: []string{"S1"}, MinDP: 10})
	if err != nil {
		t.Fatal(err)
	}
	defer w.Discard()

	err = w.WriteCoverage(CoverageBlock{SampleID: "S1", Chrom: "chr1", Start: 1, End: 100, MinDP: 30})
	if err == nil {
		t.Fatal("a coverage block was accepted by a store that will not record one")
	}
	if !strings.Contains(err.Error(), "WriterOpts.Coverage") {
		t.Errorf("the error does not say how to fix it: %v", err)
	}
}

// A backwards block is refused rather than stored.
func TestABackwardsCoverageBlockIsRefused(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	w, err := NewWriter(base, WriterOpts{Samples: []string{"S1"}, MinDP: 10, Coverage: true})
	if err != nil {
		t.Fatal(err)
	}
	defer w.Discard()

	if err := w.WriteCoverage(CoverageBlock{
		SampleID: "S1", Chrom: "chr1", Start: 500, End: 100, MinDP: 30,
	}); err == nil {
		t.Fatal("a block ending before it starts was accepted")
	}
}
