package bed

import (
	"fmt"
	"reflect"
	"sort"
	"strings"
	"testing"
)

func runSetOp(t *testing.T, aStr, bStr string, opts SetOpts) []*OutSegment {
	t.Helper()
	ra, err := NewBedReader(strings.NewReader(aStr))
	if err != nil {
		t.Fatal(err)
	}
	rb, err := NewBedReader(strings.NewReader(bStr))
	if err != nil {
		t.Fatal(err)
	}
	var segs []*OutSegment
	for s, err := range SetOperation(ra, rb, opts) {
		if err != nil {
			t.Fatalf("SetOperation: %v", err)
		}
		segs = append(segs, s)
	}
	return segs
}

func coords(segs []*OutSegment) []string {
	out := []string{}
	for _, s := range segs {
		out = append(out, fmt.Sprintf("%s:%d-%d", s.Ref, s.Start, s.End))
	}
	return out
}

func coordsStrand(segs []*OutSegment) []string {
	out := []string{}
	for _, s := range segs {
		out = append(out, fmt.Sprintf("%s:%d-%d:%s", s.Ref, s.Start, s.End, s.Strand))
	}
	return out
}

func TestSetOpCoordinates(t *testing.T) {
	cases := []struct {
		name string
		a, b string
		op   SetOp
		want []string
	}{
		// disjoint
		{"disjoint/inter", "chr1\t0\t10", "chr1\t20\t30", OpInter, []string{}},
		{"disjoint/union", "chr1\t0\t10", "chr1\t20\t30", OpUnion, []string{"chr1:0-10", "chr1:20-30"}},
		{"disjoint/sub", "chr1\t0\t10", "chr1\t20\t30", OpSub, []string{"chr1:0-10"}},
		{"disjoint/xor", "chr1\t0\t10", "chr1\t20\t30", OpXor, []string{"chr1:0-10", "chr1:20-30"}},
		// overlap
		{"overlap/inter", "chr1\t0\t20", "chr1\t10\t30", OpInter, []string{"chr1:10-20"}},
		{"overlap/union", "chr1\t0\t20", "chr1\t10\t30", OpUnion, []string{"chr1:0-30"}},
		{"overlap/sub", "chr1\t0\t20", "chr1\t10\t30", OpSub, []string{"chr1:0-10"}},
		{"overlap/xor", "chr1\t0\t20", "chr1\t10\t30", OpXor, []string{"chr1:0-10", "chr1:20-30"}},
		// abutting (half-open: must NOT intersect; union coalesces)
		{"abut/inter", "chr1\t0\t10", "chr1\t10\t20", OpInter, []string{}},
		{"abut/union", "chr1\t0\t10", "chr1\t10\t20", OpUnion, []string{"chr1:0-20"}},
		{"abut/xor", "chr1\t0\t10", "chr1\t10\t20", OpXor, []string{"chr1:0-20"}},
		// nested
		{"nested/inter", "chr1\t0\t30", "chr1\t10\t20", OpInter, []string{"chr1:10-20"}},
		{"nested/sub", "chr1\t0\t30", "chr1\t10\t20", OpSub, []string{"chr1:0-10", "chr1:20-30"}},
		{"nested/xor", "chr1\t0\t30", "chr1\t10\t20", OpXor, []string{"chr1:0-10", "chr1:20-30"}},
		// identical
		{"ident/inter", "chr1\t5\t15", "chr1\t5\t15", OpInter, []string{"chr1:5-15"}},
		{"ident/sub", "chr1\t5\t15", "chr1\t5\t15", OpSub, []string{}},
		{"ident/xor", "chr1\t5\t15", "chr1\t5\t15", OpXor, []string{}},
		// A self-overlap must flatten
		{"selfoverlap/inter", "chr1\t0\t15\nchr1\t10\t25", "chr1\t20\t30", OpInter, []string{"chr1:20-25"}},
		{"selfoverlap/union", "chr1\t0\t15\nchr1\t10\t25", "chr1\t20\t30", OpUnion, []string{"chr1:0-30"}},
		{"selfoverlap/sub", "chr1\t0\t15\nchr1\t10\t25", "chr1\t20\t30", OpSub, []string{"chr1:0-20"}},
		{"selfoverlap/xor", "chr1\t0\t15\nchr1\t10\t25", "chr1\t20\t30", OpXor, []string{"chr1:0-20", "chr1:25-30"}},
		// zero-length input is dropped
		{"zerolen/inter", "chr1\t10\t10", "chr1\t0\t20", OpInter, []string{}},
		{"zerolen/sub", "chr1\t10\t10", "chr1\t0\t20", OpSub, []string{}},
		{"zerolen/union", "chr1\t10\t10", "chr1\t0\t20", OpUnion, []string{"chr1:0-20"}},
		// multiple B inside an A gap
		{"multiB/sub", "chr1\t0\t100", "chr1\t10\t20\nchr1\t30\t40", OpSub, []string{"chr1:0-10", "chr1:20-30", "chr1:40-100"}},
		{"multiB/inter", "chr1\t0\t100", "chr1\t10\t20\nchr1\t30\t40", OpInter, []string{"chr1:10-20", "chr1:30-40"}},
		// empty inputs
		{"emptyA/sub", "", "chr1\t0\t10", OpSub, []string{}},
		{"emptyA/union", "", "chr1\t0\t10", OpUnion, []string{"chr1:0-10"}},
		{"emptyB/sub", "chr1\t0\t10", "", OpSub, []string{"chr1:0-10"}},
		// multi-chrom: chrom only in A
		{"chromOnlyA/inter", "chr1\t0\t10\nchr2\t0\t10", "chr1\t5\t15", OpInter, []string{"chr1:5-10"}},
		{"chromOnlyA/sub", "chr1\t0\t10\nchr2\t0\t10", "chr1\t5\t15", OpSub, []string{"chr1:0-5", "chr2:0-10"}},
		{"chromOnlyA/union", "chr1\t0\t10\nchr2\t0\t10", "chr1\t5\t15", OpUnion, []string{"chr1:0-15", "chr2:0-10"}},
		// multi-chrom: chrom only in B
		{"chromOnlyB/inter", "chr1\t0\t10", "chr1\t5\t15\nchr2\t0\t10", OpInter, []string{"chr1:5-10"}},
		{"chromOnlyB/union", "chr1\t0\t10", "chr1\t5\t15\nchr2\t0\t10", OpUnion, []string{"chr1:0-15", "chr2:0-10"}},
		{"chromOnlyB/xor", "chr1\t0\t10", "chr1\t5\t15\nchr2\t0\t10", OpXor, []string{"chr1:0-5", "chr1:10-15", "chr2:0-10"}},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			segs := runSetOp(t, c.a, c.b, SetOpts{Op: c.op, IgnoreStrand: true})
			if got := coords(segs); !reflect.DeepEqual(got, c.want) {
				t.Errorf("got %v, want %v", got, c.want)
			}
		})
	}
}

func TestSetOpFlanking(t *testing.T) {
	// gap of exactly 5 between [0,10) and [15,25).
	if got := coords(runSetOp(t, "chr1\t0\t10", "chr1\t15\t25", SetOpts{Op: OpUnion, IgnoreStrand: true, FlankX: 5})); !reflect.DeepEqual(got, []string{"chr1:0-25"}) {
		t.Errorf("flank=5: got %v, want [chr1:0-25]", got)
	}
	if got := coords(runSetOp(t, "chr1\t0\t10", "chr1\t15\t25", SetOpts{Op: OpUnion, IgnoreStrand: true, FlankX: 4})); !reflect.DeepEqual(got, []string{"chr1:0-10", "chr1:15-25"}) {
		t.Errorf("flank=4: got %v, want two segments", got)
	}
	// output coordinates are the true (unpadded) extents.
	if got := coords(runSetOp(t, "chr1\t0\t10\nchr1\t15\t25", "", SetOpts{Op: OpUnion, IgnoreStrand: true, FlankX: 5})); !reflect.DeepEqual(got, []string{"chr1:0-25"}) {
		t.Errorf("flank true-coords: got %v, want [chr1:0-25]", got)
	}
}

func TestSetOpStrandAware(t *testing.T) {
	aPlus := "chr1\t0\t20\tx\t0\t+"
	bMinus := "chr1\t10\t30\ty\t0\t-"
	bPlus := "chr1\t10\t30\ty\t0\t+"
	aNone := "chr1\t0\t20\tx\t0\t."

	// different strands never intersect
	if got := coords(runSetOp(t, aPlus, bMinus, SetOpts{Op: OpInter})); len(got) != 0 {
		t.Errorf("cross-strand inter: got %v, want none", got)
	}
	// same strand intersects, strand preserved
	if got := coordsStrand(runSetOp(t, aPlus, bPlus, SetOpts{Op: OpInter})); !reflect.DeepEqual(got, []string{"chr1:10-20:+"}) {
		t.Errorf("same-strand inter: got %v", got)
	}
	// sub: minus B does not subtract from plus A
	if got := coordsStrand(runSetOp(t, aPlus, bMinus, SetOpts{Op: OpSub})); !reflect.DeepEqual(got, []string{"chr1:0-20:+"}) {
		t.Errorf("cross-strand sub: got %v", got)
	}
	// union keeps both strands as separate records, sorted by start
	if got := coordsStrand(runSetOp(t, aPlus, bMinus, SetOpts{Op: OpUnion})); !reflect.DeepEqual(got, []string{"chr1:0-20:+", "chr1:10-30:-"}) {
		t.Errorf("cross-strand union: got %v", got)
	}
	// "." only matches "."
	if got := coords(runSetOp(t, aNone, bPlus, SetOpts{Op: OpInter})); len(got) != 0 {
		t.Errorf("dot vs plus inter: got %v, want none", got)
	}
	// --ignore-strand collapses → they intersect, strand dropped (StrandNone)
	if got := coordsStrand(runSetOp(t, aPlus, bMinus, SetOpts{Op: OpInter, IgnoreStrand: true})); !reflect.DeepEqual(got, []string{"chr1:10-20:."}) {
		t.Errorf("ignore-strand inter: got %v", got)
	}
}

func TestSetOpContribs(t *testing.T) {
	// inter: segment covered by both → contribs from A and B.
	segs := runSetOp(t, "chr1\t0\t20\tA1\t2\t+", "chr1\t10\t30\tB1\t3\t+", SetOpts{Op: OpInter})
	if len(segs) != 1 {
		t.Fatalf("want 1 segment, got %d", len(segs))
	}
	names, sum := contribSummary(segs[0])
	if !reflect.DeepEqual(names, []string{"A1", "B1"}) {
		t.Errorf("inter names = %v, want [A1 B1]", names)
	}
	if sum != 5 || len(segs[0].Contribs) != 2 {
		t.Errorf("inter sum=%v count=%d, want 5 and 2", sum, len(segs[0].Contribs))
	}

	// union coalesced run accumulates contributors from the whole run.
	segs = runSetOp(t, "chr1\t0\t20\tA1\t1\t+", "chr1\t15\t30\tB1\t1\t+", SetOpts{Op: OpUnion})
	if len(segs) != 1 || segs[0].Start != 0 || segs[0].End != 30 {
		t.Fatalf("union coalesce: %v", coords(segs))
	}
	names, sum = contribSummary(segs[0])
	if !reflect.DeepEqual(names, []string{"A1", "B1"}) || sum != 2 {
		t.Errorf("union contribs names=%v sum=%v, want [A1 B1] 2", names, sum)
	}

	// sub: A−B output bases are A's only → contrib only from A.
	segs = runSetOp(t, "chr1\t0\t20\tA1\t7\t+", "chr1\t10\t30\tB1\t3\t+", SetOpts{Op: OpSub})
	if len(segs) != 1 {
		t.Fatalf("sub want 1 segment, got %v", coords(segs))
	}
	names, _ = contribSummary(segs[0])
	if !reflect.DeepEqual(names, []string{"A1"}) || segs[0].Contribs[0].Source != 0 {
		t.Errorf("sub contribs = %v (sources), want only A1 from source 0", segs[0].Contribs)
	}
}

func contribSummary(s *OutSegment) (names []string, sum float64) {
	for _, c := range s.Contribs {
		names = append(names, c.Name)
		sum += c.Score
	}
	sort.Strings(names)
	return names, sum
}
