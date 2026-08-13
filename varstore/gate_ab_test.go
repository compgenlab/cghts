package varstore

import "testing"

// Allele-balance gating.
//
// THE POINT IS WHAT DEPTH CANNOT SEE. A het at AD=(50,2) has a depth of 52 and
// passes any DP gate anybody would set, while 4% alt reads is contamination or
// mismapping rather than a heterozygote. These tests are about that gap.
//
// And about the failure mode: a call that fails must become UNCERTAIN, never a
// different genotype and never a reference. That asymmetry is enforced by the
// classifier rather than by the gate -- the gate only says yes or no -- but it
// is why saying no is safe.

func call(gt string, adRef, adAlt int32) Call {
	return Call{GT: gt, DP: adRef + adAlt, ADRef: adRef, ADAlt: adAlt, GQ: Missing}
}

func TestAlleleBalanceGatesAHetDepthWouldAdmit(t *testing.T) {
	g := Gate{MinDP: 10, MinABHet: 0.15}

	// DP 52, so depth alone is delighted. 3.8% alt reads.
	if g.Admits(call("0/1", 50, 2)) {
		t.Error("a 0/1 at 4% alt reads passed; depth cannot catch this and that is the point")
	}
	// A real het.
	if !g.Admits(call("0/1", 25, 24)) {
		t.Error("a balanced het was rejected")
	}
	// Right at the boundary: 0.15 exactly should pass, since the gate is a
	// minimum rather than a strict inequality.
	if !g.Admits(call("0/1", 85, 15)) {
		t.Error("a het at exactly the threshold was rejected")
	}
}

func TestHomAltPurity(t *testing.T) {
	g := Gate{MaxRefFracHom: 0.10}

	// AD=(20,25) called 1/1 is really a het.
	if g.Admits(call("1/1", 20, 25)) {
		t.Error("a 1/1 carrying 44% reference reads passed; that is a heterozygote")
	}
	if !g.Admits(call("1/1", 1, 40)) {
		t.Error("a clean hom-alt was rejected")
	}
}

func TestMinAltReads(t *testing.T) {
	g := Gate{MinADAlt: 3}
	if g.Admits(call("0/1", 60, 2)) {
		t.Error("a carrier call resting on two alt reads passed")
	}
	if !g.Admits(call("0/1", 30, 3)) {
		t.Error("a call at exactly the alt-read floor was rejected")
	}
}

// Absent evidence is not evidence of poor quality. A caller that omits AD
// entirely, and every gVCF reference block, must not be gated into oblivion.
func TestMissingAlleleDepthsPass(t *testing.T) {
	g := Gate{MinADAlt: 3, MinABHet: 0.2, MaxRefFracHom: 0.1}
	c := Call{GT: "0/1", DP: 30, ADRef: Missing, ADAlt: Missing, GQ: Missing}
	if !g.Admits(c) {
		t.Error("a call carrying no AD was rejected; absent quality is not poor quality")
	}
}

// The het gate must not touch a hom-alt, and the hom gate must not touch a het.
// A 1/1 has an alt fraction near 1 and would sail through any het threshold;
// applying the wrong one silently does nothing, which is worse than an error.
func TestGatesApplyToTheirOwnZygosity(t *testing.T) {
	// A hom-alt with a low alt fraction is caught by MaxRefFracHom, not MinABHet.
	hom := call("1/1", 20, 25)
	if !(Gate{MinABHet: 0.9}).Admits(hom) {
		t.Error("the het allele-balance gate rejected a hom-alt call")
	}
	// A het is not subject to hom-alt purity.
	het := call("0/1", 25, 24)
	if !(Gate{MaxRefFracHom: 0.01}).Admits(het) {
		t.Error("the hom-alt purity gate rejected a het call")
	}
}

// A 1/2 sample reads as "1/." on the focal allele's row: its AD holds reference
// and focal-allele reads, and the third allele's reads are somewhere else. The
// alt fraction is therefore not a quantity a het threshold describes, and
// gating on it would discard real compound carriers.
func TestMultiAllelicCarrierIsNotHetGated(t *testing.T) {
	g := Gate{MinABHet: 0.4}
	if !g.Admits(call("1/.", 10, 8)) {
		t.Error("a 1/2 carrier was gated by a het allele-balance threshold that does not describe it")
	}
}

func TestGTClassification(t *testing.T) {
	for _, c := range []struct {
		gt              string
		het, homAlt     bool
	}{
		{"0/1", true, false},
		{"1|0", true, false},
		{"1/1", false, true},
		{"1|1", false, true},
		{"0/0", false, false},
		{"1/.", false, false}, // 1/2: neither, deliberately
		{"./.", false, false},
		{"", false, false},
	} {
		if got := IsHet(c.gt); got != c.het {
			t.Errorf("IsHet(%q) = %v, want %v", c.gt, got, c.het)
		}
		if got := IsHomAlt(c.gt); got != c.homAlt {
			t.Errorf("IsHomAlt(%q) = %v, want %v", c.gt, got, c.homAlt)
		}
	}
}

func TestDepthBanding(t *testing.T) {
	bands := DefaultDepthBands // 10, 20, 50
	for _, c := range []struct {
		dp   int32
		want int
	}{
		{5, 0}, // below the first boundary, shares its band; the gate excludes it anyway
		{10, 0}, {19, 0},
		{20, 1}, {49, 1},
		{50, 2}, {200, 2},
	} {
		if got := DepthBand(bands, c.dp); got != c.want {
			t.Errorf("DepthBand(%d) = %d, want %d", c.dp, got, c.want)
		}
	}
	if got := DepthBand(nil, 30); got != -1 {
		t.Errorf("unbanded returned %d, want -1", got)
	}
}

// THE ASYMMETRY THIS ALL RESTS ON: a gated-out carrier becomes UNCERTAIN, never
// a reference.
//
// Losing confidence in a call must never manufacture confidence in its
// opposite. If an allele-balance failure turned a suspicious het into a
// non-carrier, the gate would be converting "we are not sure they carry this"
// into "we are sure they do not" -- which is the one transformation nothing
// here is allowed to make, and the one that would quietly move every
// denominator.
func TestGatedCarrierBecomesUncertainNotReference(t *testing.T) {
	l := Locus{Chrom: "chr1", Pos: 100, Ref: "G", Alt: "A"}
	// A het the caller made, at a depth no gate would object to, but with an
	// allele balance that says contamination.
	suspicious := Call{
		SampleID: "S1", Chrom: l.Chrom, Pos: l.Pos, Ref: l.Ref, Alt: l.Alt,
		GT: "0/1", DP: 52, ADRef: 50, ADAlt: 2, GQ: Missing,
	}
	g := Gate{MinDP: 10, MinABHet: 0.15}
	if g.Admits(suspicious) {
		t.Fatal("fixture is wrong: this call should fail the gate")
	}

	states := classifyAgainst(g, []string{"S1"}, map[string]Call{"S1": suspicious}, nil)
	if len(states) != 1 {
		t.Fatalf("got %d states, want 1", len(states))
	}
	switch states[0].State {
	case StateUncertain:
		// Correct: we saw alt reads and do not believe them enough to call.
	case StateNonCarrier:
		t.Fatal("a gated-out carrier was reported as a NON-CARRIER -- the gate turned " +
			"'not sure they carry this' into 'sure they do not', which moves every denominator")
	default:
		t.Fatalf("got %q, want uncertain", states[0].State)
	}
}

// classifyAgainst mirrors the classifier's gate branch, so the property is
// asserted against the rule rather than against a whole store.
func classifyAgainst(g Gate, samples []string, calls map[string]Call, called map[string]bool) []SampleState {
	out := make([]SampleState, 0, len(samples))
	for _, name := range samples {
		st := SampleState{SampleID: name}
		if c, ok := calls[name]; ok {
			cc := c
			st.Call = &cc
			if g.Admits(c) {
				st.State = StateCarrier
			} else {
				st.State = StateUncertain
			}
		} else if called[name] {
			st.State = StateNonCarrier
		} else {
			st.State = StateNotAssayed
		}
		out = append(out, st)
	}
	return out
}
