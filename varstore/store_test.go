package varstore

import "testing"

// Gate.Admits fails open on purpose: an absent value is not evidence of poor
// quality, and a VCF with no DP column should still return its ALT calls under
// --min-dp. That is right for an ALT call, which is evidence of itself, and
// wrong for a reference call, which claims something was observed -- so the
// reference paths do not use this. See blockVouchesReference.
func TestGateAdmits(t *testing.T) {
	cases := []struct {
		name string
		gate Gate
		call Call
		want bool
	}{
		{"no gate admits anything", Gate{}, Call{DP: 1, GQ: 1}, true},
		{"no gate admits missing", Gate{}, Call{DP: Missing, GQ: Missing}, true},

		{"dp above", Gate{MinDP: 10}, Call{DP: 30, GQ: Missing}, true},
		{"dp equal", Gate{MinDP: 10}, Call{DP: 10, GQ: Missing}, true},
		{"dp below", Gate{MinDP: 10}, Call{DP: 3, GQ: Missing}, false},
		{"dp zero", Gate{MinDP: 10}, Call{DP: 0, GQ: Missing}, false},

		// Fail-open: absent depth passes. Deliberate, and the reason a gate can
		// silently do nothing against a source that carries no DP.
		{"dp missing fails open", Gate{MinDP: 10}, Call{DP: Missing, GQ: Missing}, true},

		// MinDP wins over DP when the call carries one: a gVCF block records a
		// floor across its whole span and no depth at any single base, so gating
		// on the recorded DP would assert the depth at POS across thousands of
		// bases that were never that well covered.
		{"min_dp preferred over dp", Gate{MinDP: 10}, Call{MinDP: 30, DP: 3}, true},
		{"min_dp below rejects despite dp", Gate{MinDP: 10}, Call{MinDP: 3, DP: 99}, false},

		// The documented hole this convention leaves. MinDP == 0 is
		// indistinguishable from unknown here, so Admits falls back to DP; a
		// block-derived call pins DP to Missing, and Missing passes. This is why
		// reference calls are gated by blockVouchesReference instead -- see
		// TestGvcfZeroMinDPNeverVouchesForReference.
		{"min_dp zero falls back to dp", Gate{MinDP: 10}, Call{MinDP: 0, DP: Missing}, true},

		{"gq above", Gate{MinGQ: 20}, Call{DP: Missing, GQ: 50}, true},
		{"gq below", Gate{MinGQ: 20}, Call{DP: Missing, GQ: 5}, false},
		{"gq missing fails open", Gate{MinGQ: 20}, Call{DP: Missing, GQ: Missing}, true},

		{"both must pass", Gate{MinDP: 10, MinGQ: 20}, Call{DP: 30, GQ: 5}, false},
		{"both pass", Gate{MinDP: 10, MinGQ: 20}, Call{DP: 30, GQ: 50}, true},
	}
	for _, c := range cases {
		if got := c.gate.Admits(c.call); got != c.want {
			t.Errorf("%s: Admits(%+v) under %+v = %v, want %v", c.name, c.call, c.gate, got, c.want)
		}
	}
}

// blockVouchesReference is the opposite rule, and the asymmetry is the point:
// here an unknown or zero depth CANNOT establish a reference call.
func TestBlockVouchesReference(t *testing.T) {
	cases := []struct {
		name string
		gate Gate
		call Call
		want bool
	}{
		{"no gate vouches", Gate{}, Call{MinDP: 0, GQ: Missing}, true},
		{"above", Gate{MinDP: 10}, Call{MinDP: 30, GQ: Missing}, true},
		{"equal", Gate{MinDP: 10}, Call{MinDP: 10, GQ: Missing}, true},
		{"below", Gate{MinDP: 10}, Call{MinDP: 3, GQ: Missing}, false},

		// The bug. Zero coverage and unknown coverage both fail, and neither
		// needs telling apart: neither vouches for anything.
		{"zero", Gate{MinDP: 10}, Call{MinDP: 0, DP: Missing, GQ: Missing}, false},
		{"unknown", Gate{MinDP: 10}, Call{MinDP: 0, DP: 99, GQ: Missing}, false},

		{"gq below", Gate{MinGQ: 20}, Call{MinDP: 99, GQ: 5}, false},
		{"gq missing cannot vouch", Gate{MinGQ: 20}, Call{MinDP: 99, GQ: Missing}, false},
		{"gq above", Gate{MinGQ: 20}, Call{MinDP: 99, GQ: 50}, true},
	}
	for _, c := range cases {
		if got := blockVouchesReference(c.call, c.gate); got != c.want {
			t.Errorf("%s: blockVouchesReference(%+v, %+v) = %v, want %v",
				c.name, c.call, c.gate, got, c.want)
		}
	}
}

// Span is 0-based half-open; a record's pos is 1-based and refEnd is a 0-based
// exclusive end. Both backends share this predicate, and it had no direct test.
func TestSpanOverlaps(t *testing.T) {
	s := Span{Chrom: "chr1", Start: 1000, End: 2000}
	cases := []struct {
		name   string
		chrom  string
		pos    int32
		refEnd int32
		want   bool
	}{
		{"inside", "chr1", 1500, 0, true},
		{"at first base", "chr1", 1001, 0, true},
		{"at last base", "chr1", 2000, 0, true},
		{"one before", "chr1", 1000, 0, false},
		{"one after", "chr1", 2001, 0, false},

		// refEnd 0 means unknown and must fall back to pos, never be read as a
		// zero-length record.
		{"unknown refEnd falls back to pos", "chr1", 1500, 0, true},

		// A deletion beginning before the span and reaching into it.
		{"reaches in", "chr1", 500, 1500, true},
		{"ends exactly at span start", "chr1", 500, 1000, false},
		{"ends one into span", "chr1", 500, 1001, true},
		{"ends well before", "chr1", 100, 200, false},

		// Chromosome identity is canonical, so spellings agree.
		{"ensembl spelling", "1", 1500, 0, true},
		{"refseq spelling", "NC_000001.11", 1500, 0, true},
		{"wrong chromosome", "chr2", 1500, 0, false},
	}
	for _, c := range cases {
		if got := s.Overlaps(c.chrom, c.pos, c.refEnd); got != c.want {
			t.Errorf("%s: Overlaps(%q, %d, %d) = %v, want %v",
				c.name, c.chrom, c.pos, c.refEnd, got, c.want)
		}
	}
}

func TestParseLocus(t *testing.T) {
	ok := []struct {
		in   string
		want Locus
	}{
		{"chr1:100:A:T", Locus{"chr1", 100, "A", "T"}},
		{"1:100:A:T", Locus{"1", 100, "A", "T"}},
		{"chr1:100:AT:A", Locus{"chr1", 100, "AT", "A"}},
		{"chr1:100:A:<DEL>", Locus{"chr1", 100, "A", "<DEL>"}},
		{" chr1 :100:A:T", Locus{"chr1", 100, "A", "T"}}, // chrom is trimmed
	}
	for _, c := range ok {
		got, err := ParseLocus(c.in)
		if err != nil {
			t.Errorf("ParseLocus(%q): %v", c.in, err)
			continue
		}
		if got != c.want {
			t.Errorf("ParseLocus(%q) = %+v, want %+v", c.in, got, c.want)
		}
	}

	for _, in := range []string{
		"",
		"chr1",
		"chr1:100",
		"chr1:100:A",
		"chr1:100:A:T:extra",
		"chr1:x:A:T",
		"chr1:0:A:T", // 1-based, so 0 is not a position
		"chr1:-1:A:T",
		"chr1:100::T", // ref required
		"chr1:100:A:", // alt required
	} {
		if got, err := ParseLocus(in); err == nil {
			t.Errorf("ParseLocus(%q) = %+v, want an error", in, got)
		}
	}
}

func TestParseSpan(t *testing.T) {
	// An empty region is "no restriction", not an error.
	got, err := ParseSpan("")
	if err != nil || got != nil {
		t.Errorf(`ParseSpan("") = (%+v, %v), want (nil, nil)`, got, err)
	}

	sp, err := ParseSpan("chr1:1000-2000")
	if err != nil {
		t.Fatal(err)
	}
	// htsio parses 1-based inclusive into 0-based half-open.
	if sp.Chrom != "chr1" || sp.Start != 999 || sp.End != 2000 {
		t.Errorf("ParseSpan(chr1:1000-2000) = %+v", sp)
	}

	// A bare contig is the whole thing, with the open end clamped so it stays
	// representable in the int32 the Span carries.
	whole, err := ParseSpan("chr1")
	if err != nil {
		t.Fatal(err)
	}
	if whole.Chrom != "chr1" || whole.Start != 0 || whole.End <= 0 {
		t.Errorf("ParseSpan(chr1) = %+v, want the whole contig", whole)
	}
}
