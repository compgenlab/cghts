package varstore

import (
	"path/filepath"
	"reflect"
	"testing"
)

// The fixture behind every test in this file, which is also the worked example
// the reference-call gate was argued from:
//
//	sites chr1:1000, 2000, 3000 -- every one in the catalog, S1 carries all three
//	S2  0/0 at DP 25, 40, 45    -> one run 1000-3000, min_dp 25 (band 20..49)
//	S3  0/0, run written before min_dp existed (0), conversion --min-dp 10
//
// The source VCF would report S2 as 0/0 at chr1:2000 under any gate up to 40.
// The store keeps only the run's floor, so it can vouch for S2 up to a gate of
// 25 and no further.
func buildRunGateStore(t *testing.T) string {
	t.Helper()
	base := filepath.Join(t.TempDir(), "cohort")
	w, err := NewWriter(base, WriterOpts{
		Samples:      []string{"S1", "S2", "S3"},
		MinDP:        10,
		DepthBands:   []int32{10, 20, 50},
		RowGroupSize: 4,
	})
	if err != nil {
		t.Fatal(err)
	}
	for _, pos := range []int32{1000, 2000, 3000} {
		if err := w.WriteSite(Site{
			Chrom: "chr1", Pos: pos, Ref: "G", Alt: "A", AC: 1, AN: 6, NCarriers: 1, NCalled: 3,
		}); err != nil {
			t.Fatal(err)
		}
		if err := w.WriteCall(Call{
			SampleID: "S1", Chrom: "chr1", Pos: pos, Ref: "G", Alt: "A",
			GT: "0/1", DP: 60, ADRef: 30, ADAlt: 30, GQ: 99,
		}); err != nil {
			t.Fatal(err)
		}
	}
	for _, r := range []CalledSiteRun{
		{SampleID: "S1", Chrom: "chr1", Start: 1000, End: 3000, NSites: 3, MinDP: 60},
		{SampleID: "S2", Chrom: "chr1", Start: 1000, End: 3000, NSites: 3, MinDP: 25},
		{SampleID: "S3", Chrom: "chr1", Start: 1000, End: 3000, NSites: 3, MinDP: 0},
	} {
		if err := w.WriteRegion(r); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
	return base
}

var runGateLocus = Locus{Chrom: "chr1", Pos: 2000, Ref: "G", Alt: "A"}

func statesByName(states []SampleState) map[string]State {
	out := map[string]State{}
	for _, st := range states {
		out[st.SampleID] = st.State
	}
	return out
}

// THE BUG. A run's min_dp was written but never compared with the query gate,
// so a sample whose run bottomed out at 25x was reported reference under
// --min-dp 30 -- while a 25x ALT call was correctly held below the gate. The
// VCF backend gates the same reference call, so the two disagreed.
func TestAReferenceCallIsGatedOnItsRunsDepth(t *testing.T) {
	s, err := OpenParquet(buildRunGateStore(t))
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	for _, c := range []struct {
		gate int32
		want State
	}{
		{0, StateNonCarrier},  // no gate
		{10, StateNonCarrier}, // the conversion gate
		{20, StateNonCarrier}, // a band boundary below the floor
		{25, StateNonCarrier}, // exactly the floor
		{26, StateNotAssayed}, // inside the band, above the floor
		{30, StateNotAssayed}, // the worked example
		{50, StateNotAssayed}, // the next boundary
	} {
		g := Gate{MinDP: c.gate}
		one, err := s.Classify(runGateLocus, g)
		if err != nil {
			t.Fatal(err)
		}
		many, err := s.ClassifyMany([]Locus{runGateLocus}, g)
		if err != nil {
			t.Fatal(err)
		}
		for path, states := range map[string][]SampleState{"Classify": one, "ClassifyMany": many[runGateLocus]} {
			got := statesByName(states)
			if got["S2"] != c.want {
				t.Errorf("%s gate %d: S2 = %s, want %s", path, c.gate, got["S2"], c.want)
			}
			// A failed run is NOT uncertain: that state means a below-gate ALT
			// call, and S2 has no ALT call here at all.
			if got["S2"] == StateUncertain {
				t.Errorf("%s gate %d: a refused run read as uncertain", path, c.gate)
			}
			if got["S1"] != StateCarrier {
				t.Errorf("%s gate %d: S1 = %s, want carrier", path, c.gate, got["S1"])
			}
		}
	}
}

// IncludeRef goes through a different path -- runsFor and covers rather than
// Classify's scan -- and has to reach the same answer, with the synthesized row
// carrying the floor it was judged on.
func TestIncludeRefGatesOnRunDepthAndCarriesTheFloor(t *testing.T) {
	s, err := OpenParquet(buildRunGateStore(t))
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	refsAt := func(gate int32) map[string]Call {
		calls, err := CollectCalls(s, Query{
			Loci: []Locus{runGateLocus}, IncludeRef: true, Gate: Gate{MinDP: gate},
		})
		if err != nil {
			t.Fatal(err)
		}
		out := map[string]Call{}
		for _, c := range calls {
			if c.GT == HomRefGT {
				out[c.SampleID] = c
			}
		}
		return out
	}

	got := refsAt(20)
	if c, ok := got["S2"]; !ok {
		t.Fatal("gate 20: S2's run floor is 25, so it should be reported 0/0")
	} else if c.MinDP != 25 {
		t.Errorf("gate 20: S2's reference call carries min_dp %d, want the run floor 25", c.MinDP)
	}

	if _, ok := refsAt(30)["S2"]; ok {
		t.Error("gate 30: S2 reported 0/0 from a run whose floor is 25")
	}
}

// A run written before min_dp existed reads back 0. That must fall back to the
// conversion gate -- every site in any run cleared it -- and never to "unknown",
// which would let the run pass any gate at all.
func TestARunWithNoRecordedDepthFallsBackToTheConversionGate(t *testing.T) {
	s, err := OpenParquet(buildRunGateStore(t))
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	at := func(gate int32) State {
		states, err := s.Classify(runGateLocus, Gate{MinDP: gate})
		if err != nil {
			t.Fatal(err)
		}
		return statesByName(states)["S3"]
	}
	if got := at(10); got != StateNonCarrier {
		t.Errorf("gate 10 (the conversion gate): S3 = %s, want non_carrier", got)
	}
	if got := at(11); got != StateNotAssayed {
		t.Errorf("gate 11: S3 = %s, want not_assayed -- a run with no recorded depth "+
			"vouches for the conversion gate and nothing more", got)
	}

	calls, err := CollectCalls(s, Query{
		Loci: []Locus{runGateLocus}, Samples: []string{"S3"}, IncludeRef: true,
	})
	if err != nil {
		t.Fatal(err)
	}
	if len(calls) != 1 || calls[0].MinDP != 10 {
		t.Errorf("S3's reference call = %+v, want one row carrying min_dp 10", calls)
	}
}

func TestGeometricDepthBands(t *testing.T) {
	got, err := GeometricDepthBands(10, DefaultDepthRatio, MaxDepthBand)
	if err != nil {
		t.Fatal(err)
	}
	want := []int32{10, 15, 22, 33, 50, 75, 113, 170, 256, 384, 576, 864}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("GeometricDepthBands(10, 1.5) = %v, want %v", got, want)
	}

	// A ratio near 1 at a low start floors to the same integer repeatedly; the
	// list must still be strictly increasing.
	got, err = GeometricDepthBands(2, 1.1, 20)
	if err != nil {
		t.Fatal(err)
	}
	if err := ValidateDepthBands(got); err != nil {
		t.Errorf("GeometricDepthBands(2, 1.1) = %v: %v", got, err)
	}

	for _, bad := range []struct {
		start int32
		ratio float64
	}{{0, 1.5}, {10, 1}, {10, 0.5}} {
		if _, err := GeometricDepthBands(bad.start, bad.ratio, MaxDepthBand); err == nil {
			t.Errorf("GeometricDepthBands(%d, %v) accepted", bad.start, bad.ratio)
		}
	}
}

func TestExplicitDepthBandsMustIncrease(t *testing.T) {
	for _, ok := range [][]int32{nil, {10}, {10, 20, 30}} {
		if err := ValidateDepthBands(ok); err != nil {
			t.Errorf("%v refused: %v", ok, err)
		}
	}
	for _, bad := range [][]int32{{10, 10}, {20, 10}, {0, 10}, {-5}} {
		if err := ValidateDepthBands(bad); err == nil {
			t.Errorf("%v accepted", bad)
		}
	}
	if _, err := NewWriter(filepath.Join(t.TempDir(), "s"), WriterOpts{
		Samples: []string{"S1"}, DepthBands: []int32{20, 10},
	}); err == nil {
		t.Error("NewWriter accepted a decreasing band list")
	}
}

func TestBandContaining(t *testing.T) {
	bands := []int32{10, 20, 50}
	for _, c := range []struct {
		gate   int32
		lo, hi int32
		ok     bool
	}{
		{5, 0, 0, false},  // below the conversion gate: exact
		{10, 0, 0, false}, // on a boundary: exact
		{20, 0, 0, false},
		{30, 20, 50, true}, // inside a band
		{60, 50, 0, true},  // above the last boundary: no upper edge
	} {
		lo, hi, ok := BandContaining(bands, c.gate)
		if lo != c.lo || hi != c.hi || ok != c.ok {
			t.Errorf("BandContaining(%d) = %d, %d, %v; want %d, %d, %v",
				c.gate, lo, hi, ok, c.lo, c.hi, c.ok)
		}
	}
	if _, _, ok := BandContaining(nil, 30); ok {
		t.Error("unbanded store reported a band")
	}
}
