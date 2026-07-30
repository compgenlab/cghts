package varstore

import (
	"fmt"
	"path/filepath"
	"testing"

	"github.com/parquet-go/parquet-go"
)

// buildPruneStore writes a small two-chromosome store with deliberately tiny
// row groups, so there are many groups to prune between.
func buildPruneStore(t *testing.T, rowGroup int64) (base string, want []Call) {
	t.Helper()
	base = filepath.Join(t.TempDir(), "s")
	w, err := NewWriter(base, WriterOpts{
		RowGroupSize: rowGroup,
		Samples:      []string{"S1", "S2"},
		MinDP:        10,
	})
	if err != nil {
		t.Fatal(err)
	}
	// Coordinate-sorted, as a VCF would be: chr13 then chr17, ascending pos.
	for _, chrom := range []string{"chr13", "chr17"} {
		for pos := int32(100); pos <= 2000; pos += 100 {
			for _, s := range []string{"S1", "S2"} {
				c := Call{
					SampleID: s, Chrom: chrom, Pos: pos,
					Ref: "A", Alt: "G", GT: "0/1",
					DP: 30, ADRef: 15, ADAlt: 15, GQ: 99,
				}
				if err := w.WriteCall(c); err != nil {
					t.Fatal(err)
				}
				want = append(want, c)
			}
			if err := w.WriteSite(Site{
				Chrom: chrom, Pos: pos, Ref: "A", Alt: "G",
				AC: 2, AN: 4, NCarriers: 2, NCalled: 2,
			}); err != nil {
				t.Fatal(err)
			}
		}
		if err := w.WriteRegion(CalledSiteRun{
			SampleID: "S1", Chrom: chrom, Start: 100, End: 2000, NSites: 20,
		}); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	return base, want
}

// TestPruningNeverChangesTheAnswer is the invariant that makes pruning safe: it
// may only reduce how much is read. Every locus, including the first and last
// in the file where an off-by-one in the bounds check would bite, must resolve
// identically to a full scan.
func TestPruningNeverChangesTheAnswer(t *testing.T) {
	base, _ := buildPruneStore(t, 4)
	s, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	for _, chrom := range []string{"chr13", "chr17"} {
		for pos := int32(100); pos <= 2000; pos += 100 {
			l := Locus{Chrom: chrom, Pos: pos, Ref: "A", Alt: "G"}

			// Full scan, deliberately bypassing the pruning path.
			var full []Call
			if err := scanParquet(CallsPath(base), func(c Call) bool {
				if SameLocus(c.Locus(), l) {
					full = append(full, c)
				}
				return true
			}); err != nil {
				t.Fatal(err)
			}

			pruned, err := CollectCalls(s, Query{Loci: []Locus{l}})
			if err != nil {
				t.Fatal(err)
			}
			if len(pruned) != len(full) {
				t.Errorf("%s: pruned found %d carriers, full scan found %d",
					l, len(pruned), len(full))
			}
			if len(full) != 2 {
				t.Errorf("%s: fixture should have 2 carriers, got %d", l, len(full))
			}
		}
	}
}

// TestPruningActuallySkips confirms the filters reject groups, so the test above
// is not passing merely because nothing is ever pruned.
func TestPruningActuallySkips(t *testing.T) {
	base, _ := buildPruneStore(t, 4)

	groups, kept := 0, 0
	l := Locus{Chrom: "chr17", Pos: 1000, Ref: "A", Alt: "G"}
	keep := locusFilter(l)
	if err := eachRowGroup(CallsPath(base), func(rg parquet.RowGroup) {
		groups++
		if keep(rg) {
			kept++
		}
	}); err != nil {
		t.Fatal(err)
	}
	if groups < 4 {
		t.Fatalf("fixture produced only %d row groups; nothing to prune between", groups)
	}
	if kept >= groups {
		t.Errorf("no row groups were skipped (%d of %d kept)", kept, groups)
	}
	t.Logf("kept %d of %d row groups", kept, groups)
}

// TestPruningBoundsAreInclusive pins the edges: a locus exactly on a row group's
// min or max must be kept. An exclusive comparison here would silently drop
// real carriers.
func TestPruningBoundsAreInclusive(t *testing.T) {
	base, _ := buildPruneStore(t, 4)
	for _, pos := range []int32{100, 2000} { // first and last positions written
		for _, chrom := range []string{"chr13", "chr17"} {
			l := Locus{Chrom: chrom, Pos: pos, Ref: "A", Alt: "G"}
			s, err := OpenParquet(base)
			if err != nil {
				t.Fatal(err)
			}
			got, err := CollectCalls(s, Query{Loci: []Locus{l}})
			s.Close()
			if err != nil {
				t.Fatal(err)
			}
			if len(got) != 2 {
				t.Errorf("%s at a range boundary: got %d carriers, want 2", l, len(got))
			}
		}
	}
}

// TestSpanFilterKeepsOverlap pins the half-open span conversion, where the
// 0-based/1-based mismatch is easy to get wrong by one.
func TestSpanFilterKeepsOverlap(t *testing.T) {
	base, _ := buildPruneStore(t, 4)
	s, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	// Span covering exactly position 100 (1-based) is [99,100) 0-based.
	span := &Span{Chrom: "chr17", Start: 99, End: 100}
	got, err := CollectCalls(s, Query{Samples: []string{"S1"}, Spans: []Span{*span}})
	if err != nil {
		t.Fatal(err)
	}
	if len(got) != 1 || got[0].Pos != 100 {
		t.Errorf("span [99,100) should select exactly pos 100, got %v", summarise(got))
	}
}

func summarise(cs []Call) string {
	out := ""
	for _, c := range cs {
		out += fmt.Sprintf("%s:%d ", c.Chrom, c.Pos)
	}
	return out
}
