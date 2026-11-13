package varstore

import (
	"context"
	"path/filepath"
	"reflect"
	"testing"
)

func coarsen(t *testing.T, c Coarsener, steps []any) []RefBlock {
	t.Helper()
	var out []RefBlock
	c.Emit = func(b RefBlock) error { out = append(out, b); return nil }
	for _, s := range steps {
		var err error
		switch v := s.(type) {
		case RefBlock:
			err = c.Block(v)
		case string: // "barrier"
			err = c.Barrier()
		}
		if err != nil {
			t.Fatal(err)
		}
	}
	if err := c.Flush(); err != nil {
		t.Fatal(err)
	}
	return out
}

func blk(chrom string, start, end, minDP int32) RefBlock {
	return RefBlock{SampleID: "S1", Chrom: chrom, Start: start, End: end, MinDP: minDP, GQ: Missing}
}

func TestCoarsenerMergesWithinABandAndBreaksAcrossOne(t *testing.T) {
	c := Coarsener{MinDP: 10, Bands: []int32{10, 20, 50}, MaxGap: 10}
	got := coarsen(t, c, []any{
		blk("chr1", 100, 199, 12), // band 10-19
		blk("chr1", 200, 299, 18), // same band, abutting: merges
		blk("chr1", 300, 399, 25), // band 20-49: breaks
		blk("chr1", 400, 499, 41), // same band: merges
	})
	want := []RefBlock{blk("chr1", 100, 299, 12), blk("chr1", 300, 499, 25)}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %+v\nwant %+v", got, want)
	}
}

// A gap of exactly MaxGap bases is bridged; one more is not. The gap counts the
// uncovered bases strictly between two blocks.
func TestCoarsenerGapToleranceIsInclusive(t *testing.T) {
	c := Coarsener{MinDP: 10, Bands: []int32{10, 20, 50}, MaxGap: 10}
	exact := coarsen(t, c, []any{blk("chr1", 100, 199, 30), blk("chr1", 210, 299, 30)}) // 200..209: 10 bases
	if want := []RefBlock{blk("chr1", 100, 299, 30)}; !reflect.DeepEqual(exact, want) {
		t.Errorf("a gap of exactly max-gap: got %+v, want one merged block", exact)
	}
	over := coarsen(t, c, []any{blk("chr1", 100, 199, 30), blk("chr1", 211, 299, 30)}) // 11 bases
	if len(over) != 2 {
		t.Errorf("a gap of max-gap+1 was bridged: %+v", over)
	}
}

// A block below the gate is dropped and its bases count as uncovered: short ones
// are bridged like any gap, long ones break the span.
func TestCoarsenerDropsSubGateBlocks(t *testing.T) {
	c := Coarsener{MinDP: 10, Bands: []int32{10, 20, 50}, MaxGap: 10}
	got := coarsen(t, c, []any{
		blk("chr1", 100, 199, 30),
		blk("chr1", 200, 204, 3), // 5 bases at 3x: dropped, gap of 5 bridged
		blk("chr1", 205, 299, 30),
		blk("chr1", 300, 499, 4), // 200 bases at 4x: dropped, gap too long
		blk("chr1", 500, 599, 30),
	})
	want := []RefBlock{blk("chr1", 100, 299, 30), blk("chr1", 500, 599, 30)}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %+v\nwant %+v", got, want)
	}
}

// THE RULE THAT MAKES COVERAGE SAFE. A variant record is where the sample may
// carry something, and its gap is what keeps a block from claiming reference
// there. A one-base SNV gap is well inside any tolerance, so this is exactly the
// case a gap merge would get wrong if it were allowed to.
func TestCoarsenerNeverBridgesAVariant(t *testing.T) {
	c := Coarsener{MinDP: 10, Bands: []int32{10, 20, 50}, MaxGap: 10}
	got := coarsen(t, c, []any{
		blk("chr1", 100, 199, 30),
		"barrier", // the sample's own SNV at 200
		blk("chr1", 201, 299, 30),
	})
	want := []RefBlock{blk("chr1", 100, 199, 30), blk("chr1", 201, 299, 30)}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("a block was merged across the sample's own variant: %+v", got)
	}
}

func TestCoarsenerBreaksAtChromosomesAndKeepsTheFloor(t *testing.T) {
	c := Coarsener{MinDP: 10, MaxGap: 1000} // unbanded: depth never breaks a span
	a := blk("chr1", 100, 199, 40)
	a.GQ = 60
	b := blk("chr1", 200, 299, 15)
	b.GQ = 30
	got := coarsen(t, c, []any{a, b, blk("chr2", 100, 199, 40)})
	if len(got) != 2 {
		t.Fatalf("got %+v, want one span per chromosome", got)
	}
	if got[0].MinDP != 15 || got[0].GQ != 30 || got[0].End != 299 {
		t.Errorf("merged span = %+v, want end 299, min_dp 15 (the floor of its parts), gq 30", got[0])
	}
	if got[1].Chrom != "chr2" {
		t.Errorf("second span on %s, want chr2", got[1].Chrom)
	}
}

func TestCoarsenerPassthroughKeepsSourceBlocks(t *testing.T) {
	c := Coarsener{MinDP: 10, Bands: []int32{10, 20, 50}, MaxGap: 10, Passthrough: true}
	in := []any{blk("chr1", 100, 199, 30), blk("chr1", 200, 299, 31), blk("chr1", 300, 399, 5)}
	got := coarsen(t, c, in)
	want := []RefBlock{blk("chr1", 100, 199, 30), blk("chr1", 200, 299, 31)}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("passthrough: got %+v, want the source's own blocks above the gate", got)
	}
}

func TestCoarsenerRefusesOverlappingBlocks(t *testing.T) {
	c := Coarsener{MinDP: 10, MaxGap: 10, Emit: func(RefBlock) error { return nil }}
	if err := c.Block(blk("chr1", 100, 199, 30)); err != nil {
		t.Fatal(err)
	}
	if err := c.Block(blk("chr1", 150, 250, 30)); err == nil {
		t.Error("an overlapping block was accepted")
	}
}

func TestFragmentRoundTrip(t *testing.T) {
	dir := filepath.Join(t.TempDir(), "S1.frag")
	w, err := NewFragmentWriter(dir, FragmentMeta{
		SampleID: "S1", Source: "s1.g.vcf.gz", MinDP: 10,
		DepthBands: []int32{10, 15, 22}, DepthRatio: 1.5, MaxGap: 10,
	})
	if err != nil {
		t.Fatal(err)
	}
	in := []RefBlock{blk("chr1", 100, 199, 30), blk("chr1", 300, 399, 12), blk("chr2", 50, 60, 20)}
	for _, b := range in {
		if err := w.Write(b); err != nil {
			t.Fatal(err)
		}
	}

	// Not finished, so not a fragment.
	if _, err := ReadFragmentMeta(context.Background(), dir); err == nil {
		t.Fatal("an unfinished fragment read as complete")
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}

	m, err := ReadFragmentMeta(context.Background(), dir)
	if err != nil {
		t.Fatal(err)
	}
	if m.SampleID != "S1" || len(m.Chroms) != 2 || m.Chroms[0].Rows != 2 || m.Chroms[1].Name != "chr2" {
		t.Errorf("meta = %+v", m)
	}
	if m.Chroms[0].First != 100 || m.Chroms[0].Last != 399 {
		t.Errorf("chr1 bounds %d-%d, want 100-399", m.Chroms[0].First, m.Chroms[0].Last)
	}
	var got []RefBlock
	for _, c := range m.Chroms {
		if err := FragmentBlocks(context.Background(), dir, c, func(b RefBlock) bool {
			got = append(got, b)
			return true
		}); err != nil {
			t.Fatal(err)
		}
	}
	if !reflect.DeepEqual(got, in) {
		t.Errorf("read back %+v\nwant %+v", got, in)
	}
}

func TestFragmentWriterRefusesDisorder(t *testing.T) {
	newW := func() *FragmentWriter {
		w, err := NewFragmentWriter(filepath.Join(t.TempDir(), "f"), FragmentMeta{SampleID: "S1", MinDP: 10})
		if err != nil {
			t.Fatal(err)
		}
		t.Cleanup(func() { w.Discard() })
		return w
	}
	w := newW()
	_ = w.Write(blk("chr1", 100, 199, 30))
	if err := w.Write(blk("chr1", 150, 160, 30)); err == nil {
		t.Error("an overlapping block was accepted")
	}
	w = newW()
	_ = w.Write(blk("chr1", 100, 199, 30))
	_ = w.Write(blk("chr2", 100, 199, 30))
	if err := w.Write(blk("chr1", 500, 599, 30)); err == nil {
		t.Error("a chromosome revisited after another was accepted")
	}
	w = newW()
	b := blk("chr1", 100, 199, 30)
	b.SampleID = "S2"
	if err := w.Write(b); err == nil {
		t.Error("another sample's block was accepted")
	}
}
