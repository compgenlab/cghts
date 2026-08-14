package varstore

import (
	"path/filepath"
	"testing"
)

// ClassifyMany must agree with Classify, locus by locus and sample by sample.
//
// It is an optimisation, so the only thing that matters is that it changed
// nothing. A differential test against the function it replaces is the whole
// point -- and the case worth writing it for is the GATED-OUT CALL, where the
// obvious implementation (build it on Calls, which already takes a locus set)
// silently reports StateNotAssayed where Classify reports StateUncertain.
// "Nobody looked" and "we looked and cannot vouch for what we saw" are
// different claims about a person.

func manyFixture(t *testing.T) string {
	t.Helper()
	dir := filepath.Join(t.TempDir(), "store")
	w, err := NewWriter(dir, WriterOpts{Samples: []string{"S1", "S2", "S3"}, MinDP: 10})
	if err != nil {
		t.Fatal(err)
	}

	loci := []Locus{
		{Chrom: "chr1", Pos: 100, Ref: "G", Alt: "A"},
		{Chrom: "chr1", Pos: 200, Ref: "C", Alt: "T"},
		{Chrom: "chr1", Pos: 300, Ref: "A", Alt: "G"},
	}
	for _, l := range loci {
		if err := w.WriteSite(Site{
			Chrom: l.Chrom, Pos: l.Pos, Ref: l.Ref, Alt: l.Alt, AN: 6, NCalled: 3,
		}); err != nil {
			t.Fatal(err)
		}
	}

	// S1 carries at 100 with ample depth, and at 200 BELOW the gate -- the case
	// that separates uncertain from not-assayed.
	if err := w.WriteCall(Call{
		SampleID: "S1", Chrom: "chr1", Pos: 100, Ref: "G", Alt: "A",
		GT: "0/1", DP: 40, ADRef: 20, ADAlt: 20, GQ: Missing,
	}); err != nil {
		t.Fatal(err)
	}
	if err := w.WriteCall(Call{
		SampleID: "S1", Chrom: "chr1", Pos: 200, Ref: "C", Alt: "T",
		GT: "1/1", DP: 4, ADRef: 0, ADAlt: 4, GQ: Missing,
	}); err != nil {
		t.Fatal(err)
	}
	// S2 carries at 300.
	if err := w.WriteCall(Call{
		SampleID: "S2", Chrom: "chr1", Pos: 300, Ref: "A", Alt: "G",
		GT: "0/1", DP: 30, ADRef: 15, ADAlt: 15, GQ: Missing,
	}); err != nil {
		t.Fatal(err)
	}

	// S1 and S2 are covered across the whole span; S3 only reaches 150, so it is
	// reference at 100 and NOT ASSAYED at 200 and 300.
	for _, r := range []CalledSiteRun{
		{SampleID: "S1", Chrom: "chr1", Start: 50, End: 400, NSites: 3, MinDP: 30},
		{SampleID: "S2", Chrom: "chr1", Start: 50, End: 400, NSites: 3, MinDP: 30},
		{SampleID: "S3", Chrom: "chr1", Start: 50, End: 150, NSites: 1, MinDP: 30},
	} {
		if err := w.WriteRegion(r); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
	return dir
}

func TestClassifyManyAgreesWithClassify(t *testing.T) {
	s, err := OpenParquet(manyFixture(t))
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	loci := []Locus{
		{Chrom: "chr1", Pos: 100, Ref: "G", Alt: "A"},
		{Chrom: "chr1", Pos: 200, Ref: "C", Alt: "T"},
		{Chrom: "chr1", Pos: 300, Ref: "A", Alt: "G"},
		// One the store never interrogated: every sample not assayed, and the
		// run intervals must not be consulted even though they bracket it.
		{Chrom: "chr1", Pos: 250, Ref: "T", Alt: "C"},
	}
	gate := Gate{MinDP: 10}

	batched, err := s.ClassifyMany(loci, gate)
	if err != nil {
		t.Fatal(err)
	}
	for _, l := range loci {
		one, err := s.Classify(l, gate)
		if err != nil {
			t.Fatal(err)
		}
		many := batched[l]
		if len(many) != len(one) {
			t.Fatalf("%s: batched returned %d states, single %d", l, len(many), len(one))
		}
		byName := map[string]State{}
		for _, st := range many {
			byName[st.SampleID] = st.State
		}
		for _, st := range one {
			if got := byName[st.SampleID]; got != st.State {
				t.Errorf("%s %s: batched %q, single %q", l, st.SampleID, got, st.State)
			}
		}
	}
}

// The case that would have been silently wrong. S1's call at 200 is real and
// below the gate: we looked, and cannot vouch for what we saw.
func TestGatedCallIsUncertainNotAbsent(t *testing.T) {
	s, err := OpenParquet(manyFixture(t))
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	l := Locus{Chrom: "chr1", Pos: 200, Ref: "C", Alt: "T"}
	got, err := s.ClassifyMany([]Locus{l}, Gate{MinDP: 10})
	if err != nil {
		t.Fatal(err)
	}
	for _, st := range got[l] {
		if st.SampleID != "S1" {
			continue
		}
		switch st.State {
		case StateUncertain:
			return // correct
		case StateNotAssayed:
			t.Fatal("a gated-out ALT call reported as NEVER ASSAYED -- " +
				"'we cannot vouch for this' became 'nobody looked', which are different claims")
		case StateNonCarrier:
			t.Fatal("a gated-out ALT call reported as REFERENCE")
		default:
			t.Fatalf("S1 at the gated locus = %q, want uncertain", st.State)
		}
	}
	t.Fatal("S1 missing from the result entirely")
}

// A locus outside the catalog is not assayed for everyone, even where the
// callable runs bracket the position. A run says a sample was covered across a
// span; it cannot say a position the source never reported was examined.
func TestUninterrogatedLocusIgnoresRuns(t *testing.T) {
	s, err := OpenParquet(manyFixture(t))
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	l := Locus{Chrom: "chr1", Pos: 250, Ref: "T", Alt: "C"}
	got, err := s.ClassifyMany([]Locus{l}, Gate{MinDP: 10})
	if err != nil {
		t.Fatal(err)
	}
	for _, st := range got[l] {
		if st.State != StateNotAssayed {
			t.Errorf("%s at an off-catalog locus = %q, want not assayed -- "+
				"the runs bracket it, and that is not the same as having looked",
				st.SampleID, st.State)
		}
	}
}
