package htsio

import "testing"

func TestCigarQueryLen(t *testing.T) {
	cases := map[string]int{
		"*":         0,
		"":          0,
		"10M":       10,
		"5M1I4M":    10, // M+I+M consume query
		"3S10M2S":   15, // soft clips consume query
		"10M5D10M":  20, // D does not consume query
		"10M3N10M":  20, // N does not consume query
		"5H10M5H":   10, // hard clips do not consume query
		"4M2P4M":    8,  // padding does not consume query
		"100M":      100,
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
		{"10M", "ACGTACGTAC", false},   // 10 == 10
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
