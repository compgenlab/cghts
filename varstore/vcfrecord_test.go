package varstore

import "testing"

// SplitGT is the whole multiallelic recoding rule, and had no direct test. The
// load-bearing decision is that a non-focal alternate becomes "." and not "0":
// a 1/2 sample carries neither allele as reference, and writing 0 would invent a
// reference observation the data does not support.
func TestSplitGT(t *testing.T) {
	cases := []struct {
		gt      string
		altIdx  int
		want    string
		carrier bool
	}{
		{"0/1", 1, "0/1", true},
		{"1/1", 1, "1/1", true},
		{"0/0", 1, "0/0", false},
		{"./.", 1, "./.", false},

		// A 1/2 sample is a carrier of BOTH split rows, and in each the other
		// alternate is masked to "." rather than called reference.
		{"1/2", 1, "1/.", true},
		{"1/2", 2, "./1", true},

		// A sample carrying neither of the focal alleles.
		{"2/3", 1, "./.", false},
		{"0/2", 1, "0/.", false},
		{"0/2", 2, "0/1", true},

		// Phasing survives.
		{"0|1", 1, "0|1", true},
		{"1|0", 1, "1|0", true},
		{"1|2", 2, ".|1", true},

		// Haploid stays haploid.
		{"1", 1, "1", true},
		{"0", 1, "0", false},
		{".", 1, ".", false},

		// Missing alleles normalize; an empty GT is not a call.
		{"", 1, ".", false},
		{"./1", 1, "./1", true},

		// An unparseable allele is masked rather than guessed at.
		{"0/x", 1, "0/.", false},

		// Polyploid.
		{"0/1/1", 1, "0/1/1", true},
	}
	for _, c := range cases {
		got, carrier := SplitGT(c.gt, c.altIdx)
		if got != c.want || carrier != c.carrier {
			t.Errorf("SplitGT(%q, %d) = (%q, %v), want (%q, %v)",
				c.gt, c.altIdx, got, carrier, c.want, c.carrier)
		}
	}
}

// AD is taken per allele and never summed: at a multiallelic site the depth
// supporting allele 1 says nothing about allele 2.
func TestSplitAD(t *testing.T) {
	cases := []struct {
		ad             string
		altIdx         int
		wantRef, wantA int32
	}{
		{"28,2", 1, 28, 2},
		{"10,3,7", 1, 10, 3},
		{"10,3,7", 2, 10, 7},

		// Beyond the recorded alleles: unknown, not zero.
		{"10,3", 2, 10, Missing},
		// A truncated AD is not padded with zeroes.
		{"10", 1, 10, Missing},

		{"", 1, Missing, Missing},
		{".", 1, Missing, Missing},
		{".,.", 1, Missing, Missing},
		{"x,y", 1, Missing, Missing},
	}
	for _, c := range cases {
		gotRef, gotAlt := SplitAD(c.ad, c.altIdx)
		if gotRef != c.wantRef || gotAlt != c.wantA {
			t.Errorf("SplitAD(%q, %d) = (%d, %d), want (%d, %d)",
				c.ad, c.altIdx, gotRef, gotAlt, c.wantRef, c.wantA)
		}
	}
}

// The three GT predicates draw different lines, and the differences are the
// point: HasCall and IsHomRef separate "declined to call" from "called
// reference", which is what separates not-assayed from non-carrier.
func TestGTPredicates(t *testing.T) {
	cases := []struct {
		gt                       string
		altCarrier, hasCall, ref bool
	}{
		{"0/0", false, true, true},
		{"0|0", false, true, true},
		{"0", false, true, true},
		{"0/1", true, true, false},
		{"1/1", true, true, false},
		{"1/2", true, true, false},
		{"./.", false, false, false},
		{"", false, false, false},
		// A half-call: one allele known and reference. Not an all-reference
		// observation, so IsHomRef declines it.
		{"./0", false, true, false},
		{"./1", true, true, false},
	}
	for _, c := range cases {
		if got := IsAltCarrier(c.gt); got != c.altCarrier {
			t.Errorf("IsAltCarrier(%q) = %v, want %v", c.gt, got, c.altCarrier)
		}
		if got := HasCall(c.gt); got != c.hasCall {
			t.Errorf("HasCall(%q) = %v, want %v", c.gt, got, c.hasCall)
		}
		if got := IsHomRef(c.gt); got != c.ref {
			t.Errorf("IsHomRef(%q) = %v, want %v", c.gt, got, c.ref)
		}
	}
}
