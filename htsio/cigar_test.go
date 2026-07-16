package htsio

import (
	"strings"
	"testing"
)

// CigarRefLen is covered by TestCigarRefLen in sam_test.go.

func TestParseCigar(t *testing.T) {
	ok := []struct {
		cigar string
		want  []CigarOp
	}{
		{"*", nil},
		{"", nil},
		{"10M", []CigarOp{{10, 'M'}}},
		{"3S10M2S", []CigarOp{{3, 'S'}, {10, 'M'}, {2, 'S'}}},
		{"5H3S10M2S5H", []CigarOp{{5, 'H'}, {3, 'S'}, {10, 'M'}, {2, 'S'}, {5, 'H'}}},
		{"0M", []CigarOp{{0, 'M'}}}, // zero-length op is well-formed
		{"10=5X", []CigarOp{{10, '='}, {5, 'X'}}},
		{"4M2P4M", []CigarOp{{4, 'M'}, {2, 'P'}, {4, 'M'}}},
		{"100M", []CigarOp{{100, 'M'}}}, // multi-digit length
	}
	for _, c := range ok {
		got, err := ParseCigar(c.cigar)
		if err != nil {
			t.Errorf("ParseCigar(%q) unexpected error: %v", c.cigar, err)
			continue
		}
		if len(got) != len(c.want) {
			t.Errorf("ParseCigar(%q) = %v, want %v", c.cigar, got, c.want)
			continue
		}
		for i := range got {
			if got[i] != c.want[i] {
				t.Errorf("ParseCigar(%q)[%d] = %v, want %v", c.cigar, i, got[i], c.want[i])
			}
		}
	}

	bad := []string{
		"10",   // trailing length with no operation
		"M10",  // operation with no length
		"10Z",  // operation outside MIDNSHP=X
		"10B",  // B is not a SAM operation, though malformed BAM can decode to it
		"-5M",  // '-' reads as an operation with no length
		"10M5", // trailing length after a valid op
		"MM",   // no lengths at all
	}
	for _, cigar := range bad {
		if _, err := ParseCigar(cigar); err == nil {
			t.Errorf("ParseCigar(%q) expected an error, got nil", cigar)
		}
	}
}

func TestCigarOpConsumes(t *testing.T) {
	cases := []struct {
		op               byte
		wantRef, wantQry bool
	}{
		{'M', true, true},
		{'I', false, true},
		{'D', true, false},
		{'N', true, false},
		{'S', false, true},
		{'H', false, false},
		{'P', false, false},
		{'=', true, true},
		{'X', true, true},
	}
	for _, c := range cases {
		op := CigarOp{Len: 1, Op: c.op}
		if got := op.ConsumesRef(); got != c.wantRef {
			t.Errorf("CigarOp{%q}.ConsumesRef() = %v, want %v", string(c.op), got, c.wantRef)
		}
		if got := op.ConsumesQuery(); got != c.wantQry {
			t.Errorf("CigarOp{%q}.ConsumesQuery() = %v, want %v", string(c.op), got, c.wantQry)
		}
	}
}

func TestParseCigarRoundTrip(t *testing.T) {
	for _, cigar := range []string{"10M", "3S10M2S", "5H3S10M2S5H", "10M5D10M", "10=5X", "4M2P4M", "100M"} {
		ops, err := ParseCigar(cigar)
		if err != nil {
			t.Errorf("ParseCigar(%q): %v", cigar, err)
			continue
		}
		var sb strings.Builder
		for _, op := range ops {
			sb.WriteString(op.String())
		}
		if got := sb.String(); got != cigar {
			t.Errorf("round trip of %q = %q", cigar, got)
		}
	}
}

func TestCigarAlignEnd(t *testing.T) {
	cases := []struct {
		pos   int
		cigar string
		want  int
	}{
		{100, "50M", 149},      // 0-based half-open: 99+50
		{100, "50M20S", 149},   // soft clip consumes no reference
		{100, "10M5D40M", 154}, // deletion consumes reference
		{1, "10M", 10},         // contig start
		{100, "*", 99},         // unspecified CIGAR consumes no reference
		{100, "5H50M5H", 149},  // hard clips consume no reference
	}
	for _, c := range cases {
		if got := CigarAlignEnd(c.pos, c.cigar); got != c.want {
			t.Errorf("CigarAlignEnd(%d, %q) = %d, want %d", c.pos, c.cigar, got, c.want)
		}
	}
}

func TestCigarSoftClips(t *testing.T) {
	cases := []struct {
		cigar               string
		wantLead, wantTrail int
	}{
		{"50M20S", 0, 20},
		{"20S50M", 20, 0},
		{"20S50M10S", 20, 10},
		{"5H20S50M10S5H", 20, 10}, // hard clips ignored
		{"50M", 0, 0},
		{"*", 0, 0},
		{"", 0, 0},
		{"10S", 10, 0},   // a lone S is leading, never both
		{"5H10S", 10, 0}, // same, behind a hard clip
		{"10Z", 0, 0},    // malformed
	}
	for _, c := range cases {
		lead, trail := CigarSoftClips(c.cigar)
		if lead != c.wantLead || trail != c.wantTrail {
			t.Errorf("CigarSoftClips(%q) = (%d, %d), want (%d, %d)",
				c.cigar, lead, trail, c.wantLead, c.wantTrail)
		}
	}
}

func TestCigarQueryToRef(t *testing.T) {
	cases := []struct {
		cigar   string
		pos     int
		qIdx    int
		wantPos int
		wantOk  bool
	}{
		{"50M20S", 100, 45, 144, true},      // aligned
		{"50M20S", 100, 50, 0, false},       // first base of the trailing clip
		{"50M20S", 100, 60, 0, false},       // deeper in the trailing clip
		{"20S50M", 100, 19, 0, false},       // last base of the leading clip
		{"20S50M", 100, 25, 104, true},      // aligned, 5 bases in
		{"50M5I15S", 100, 52, 0, false},     // inside the insertion
		{"10M5D40M20S", 100, 30, 134, true}, // past a deletion
		{"10M5N40M20S", 100, 30, 134, true}, // past a skip
		{"5H50M20S", 100, 50, 0, false},     // H absent from SEQ; clip
		{"5H50M20S", 100, 40, 139, true},    // H absent from SEQ; aligned
		{"70M", 100, 60, 159, true},         // no clips at all
		{"*", 0, 5, 0, false},               // unspecified CIGAR
		{"50M", 100, 50, 0, false},          // qIdx past the query
		{"50M", 100, -1, 0, false},          // negative qIdx
		{"10Z", 100, 5, 0, false},           // malformed CIGAR
	}
	for _, c := range cases {
		got, ok := CigarQueryToRef(c.pos, c.cigar, c.qIdx)
		if ok != c.wantOk || (ok && got != c.wantPos) {
			t.Errorf("CigarQueryToRef(%d, %q, %d) = (%d, %v), want (%d, %v)",
				c.pos, c.cigar, c.qIdx, got, ok, c.wantPos, c.wantOk)
		}
	}

	// For an all-M CIGAR every query index maps to pos-1+qIdx.
	for q := 0; q < 50; q++ {
		got, ok := CigarQueryToRef(100, "50M", q)
		if !ok || got != 99+q {
			t.Errorf("CigarQueryToRef(100, \"50M\", %d) = (%d, %v), want (%d, true)", q, got, ok, 99+q)
		}
	}

	// Reference positions increase strictly across successive mappable indices.
	for _, cigar := range []string{"20S10M5D40M20S", "10M5I40M", "5H20S50M10S5H", "10M3N40M"} {
		prev := -1
		for q := 0; q < CigarQueryLen(cigar); q++ {
			got, ok := CigarQueryToRef(100, cigar, q)
			if !ok {
				continue
			}
			if got <= prev {
				t.Errorf("CigarQueryToRef(100, %q, %d) = %d, not strictly increasing (prev %d)",
					cigar, q, got, prev)
			}
			prev = got
		}
	}
}

// TestParseCigarConsistency pins ParseCigar against the two hand-rolled scanners
// it deliberately does not share code with.
func TestParseCigarConsistency(t *testing.T) {
	corpus := []string{
		"10M", "5M1I4M", "3S10M2S", "10M5D10M", "10M3N10M", "5H10M5H",
		"4M2P4M", "10=5X", "100M", "20S10M5D40M20S", "5H20S50M10S5H", "*", "",
	}
	for _, cigar := range corpus {
		ops, err := ParseCigar(cigar)
		if err != nil {
			t.Errorf("ParseCigar(%q): %v", cigar, err)
			continue
		}
		refLen, qryLen := 0, 0
		for _, op := range ops {
			if op.ConsumesRef() {
				refLen += op.Len
			}
			if op.ConsumesQuery() {
				qryLen += op.Len
			}
		}
		if want := CigarRefLen(cigar); refLen != want {
			t.Errorf("%q: ParseCigar ref len %d != CigarRefLen %d", cigar, refLen, want)
		}
		if want := CigarQueryLen(cigar); qryLen != want {
			t.Errorf("%q: ParseCigar query len %d != CigarQueryLen %d", cigar, qryLen, want)
		}
	}
}

func TestCigarQueryLen(t *testing.T) {
	cases := map[string]int{
		"*":        0,
		"":         0,
		"10M":      10,
		"5M1I4M":   10, // M+I+M consume query
		"3S10M2S":  15, // soft clips consume query
		"10M5D10M": 20, // D does not consume query
		"10M3N10M": 20, // N does not consume query
		"5H10M5H":  10, // hard clips do not consume query
		"4M2P4M":   8,  // padding does not consume query
		"100M":     100,
	}
	for cigar, want := range cases {
		if got := CigarQueryLen(cigar); got != want {
			t.Errorf("CigarQueryLen(%q) = %d, want %d", cigar, got, want)
		}
	}
}

func TestValidateCigarSeq(t *testing.T) {
	cases := []struct {
		cigar, seq string
		wantErr    bool
	}{
		{"10M", "ACGTACGTAC", false},    // 10 == 10
		{"5M1I4M", "ACGTAACGTA", false}, // query len 10 == 10
		{"3S7M", "ACGTACGTAC", false},   // 10 == 10
		{"*", "ACGT", false},            // unspecified CIGAR — skip
		{"10M", "*", false},             // no sequence — skip
		{"10M", "", false},              // empty sequence — skip
		{"10M", "ACGT", true},           // 10 != 4
		{"5M", "ACGTACGTAC", true},      // 5 != 10
		{"5M5D", "ACGTA", false},        // query len 5 (D excluded) == 5
		{"5M5D", "ACGTACGTAC", true},    // query len 5 != 10
	}
	for _, c := range cases {
		err := ValidateCigarSeq(c.cigar, c.seq)
		if (err != nil) != c.wantErr {
			t.Errorf("ValidateCigarSeq(%q, %q) err = %v, wantErr = %v", c.cigar, c.seq, err, c.wantErr)
		}
	}
}
