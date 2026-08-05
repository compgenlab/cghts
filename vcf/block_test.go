package vcf

import (
	"reflect"
	"testing"
)

// BlockDepth and BlockRGQ are the accessors a reference assertion rests on: a
// gVCF block claims coverage across a span, and these say how much. Getting
// "unknown" and "zero" the wrong way round here is what admits an uncovered
// region as a confident 0/0.

func TestBlockDepth(t *testing.T) {
	for _, tc := range []struct {
		name    string
		format  string
		sample  string
		want    int32
		wantOK  bool
		comment string
	}{
		{name: "MIN_DP is the block-wide floor", format: "GT:DP:MIN_DP", sample: "0/0:45:12", want: 12, wantOK: true},
		// DP alone is the depth at POS, a weaker claim than MIN_DP but the only
		// one available.
		{name: "DP is the fallback", format: "GT:DP", sample: "0/0:30", want: 30, wantOK: true},
		{name: "MIN_DP wins over DP", format: "GT:MIN_DP:DP", sample: "0/0:5:900", want: 5, wantOK: true},
		// Zero is a real claim -- the caller looked and found no reads -- and
		// must be reported as known, so a min-dp gate can reject it. Reporting
		// it as unknown is how an uncovered block passes every gate.
		{name: "zero depth is known, not unknown", format: "GT:MIN_DP", sample: "0/0:0", want: 0, wantOK: true},
		{name: "absent is unknown", format: "GT", sample: "0/0", want: 0, wantOK: false},
		// "." and "" are unknown, and fall through to DP rather than reporting
		// a spurious zero.
		{name: "missing MIN_DP falls through to DP", format: "GT:MIN_DP:DP", sample: "0/0:.:22", want: 22, wantOK: true},
		{name: "missing both is unknown", format: "GT:MIN_DP:DP", sample: "0/0:.:.", want: 0, wantOK: false},
		{name: "non-numeric is unknown", format: "GT:MIN_DP", sample: "0/0:abc", want: 0, wantOK: false},
		// A short sample column pads with "." rather than erroring, so this
		// reads as unknown too.
		{name: "truncated sample column", format: "GT:DP:MIN_DP", sample: "0/0", want: 0, wantOK: false},
	} {
		t.Run(tc.name, func(t *testing.T) {
			line := "chr1\t1000\t.\tA\t<NON_REF>\t.\t.\tEND=2000\t" + tc.format + "\t" + tc.sample
			rec := editRec(t, line)
			got, ok := rec.BlockDepth(0)
			if got != tc.want || ok != tc.wantOK {
				t.Errorf("BlockDepth(0) = (%d, %v), want (%d, %v)", got, ok, tc.want, tc.wantOK)
			}
		})
	}
}

func TestBlockRGQ(t *testing.T) {
	for _, tc := range []struct {
		name   string
		format string
		sample string
		want   int32
		wantOK bool
	}{
		// GATK writes RGQ on reference blocks in place of GQ, so a caller
		// reading only GQ finds nothing at exactly the records it cares about.
		{"RGQ preferred", "GT:RGQ", "0/0:60", 60, true},
		{"GQ is the fallback", "GT:GQ", "0/0:99", 99, true},
		{"RGQ wins over GQ", "GT:RGQ:GQ", "0/0:1:99", 1, true},
		{"zero is known", "GT:RGQ", "0/0:0", 0, true},
		{"absent is unknown", "GT", "0/0", 0, false},
		{"missing RGQ falls through to GQ", "GT:RGQ:GQ", "0/0:.:44", 44, true},
		{"non-numeric is unknown", "GT:RGQ", "0/0:x", 0, false},
	} {
		t.Run(tc.name, func(t *testing.T) {
			line := "chr1\t1000\t.\tA\t<NON_REF>\t.\t.\tEND=2000\t" + tc.format + "\t" + tc.sample
			rec := editRec(t, line)
			got, ok := rec.BlockRGQ(0)
			if got != tc.want || ok != tc.wantOK {
				t.Errorf("BlockRGQ(0) = (%d, %v), want (%d, %v)", got, ok, tc.want, tc.wantOK)
			}
		})
	}
}

// An out-of-range sample index is unknown rather than a panic or a zero claim --
// the accessors are called per sample over a roster that may not match the record.
func TestBlockAccessorsOutOfRange(t *testing.T) {
	rec := editRec(t, "chr1\t1000\t.\tA\t<NON_REF>\t.\t.\tEND=2000\tGT:MIN_DP:RGQ\t0/0:12:60")
	for _, i := range []int{-1, 1, 99} {
		if v, ok := rec.BlockDepth(i); ok {
			t.Errorf("BlockDepth(%d) = (%d, true), want unknown", i, v)
		}
		if v, ok := rec.BlockRGQ(i); ok {
			t.Errorf("BlockRGQ(%d) = (%d, true), want unknown", i, v)
		}
	}
	// A record with no sample columns at all.
	bare := editRec(t, "chr1\t1000\t.\tA\t<NON_REF>\t.\t.\tEND=2000")
	if v, ok := bare.BlockDepth(0); ok {
		t.Errorf("BlockDepth(0) on a sample-less record = (%d, true), want unknown", v)
	}
}

func TestIsRefBlockAlt(t *testing.T) {
	for _, tc := range []struct {
		alt  string
		want bool
	}{
		{"<NON_REF>", true},
		{"<*>", true},
		{"T", false},
		{"<DEL>", false},
		{".", false}, // a bare "." ALT is handled by IsRefBlock, not per-allele
		{"", false},
		{"<non_ref>", false}, // the spelling is fixed; no case folding
		{"NON_REF", false},
		{"<NON_REF", false},
	} {
		if got := IsRefBlockAlt(tc.alt); got != tc.want {
			t.Errorf("IsRefBlockAlt(%q) = %v, want %v", tc.alt, got, tc.want)
		}
	}
}

// RefSpanEnd is RefSpan's second return, and must stay that way -- the tabix
// index writer and an interval query both go through one of the two.
func TestRefSpanEndMatchesRefSpan(t *testing.T) {
	for _, tc := range spanCases {
		t.Run(tc.name, func(t *testing.T) {
			rec := editRec(t, tc.line)
			_, end := rec.RefSpan()
			if got := rec.RefSpanEnd(); got != end {
				t.Errorf("RefSpanEnd() = %d but RefSpan() end = %d", got, end)
			}
			if got := rec.RefSpanEnd(); got != tc.want {
				t.Errorf("RefSpanEnd() = %d, want %d", got, tc.want)
			}
		})
	}
}

// recordFields is the adapter that lets the shared span rules read a parsed
// record. Its two list accessors are the ones with an interesting contract.
func TestRecordFieldsAccessors(t *testing.T) {
	rec := editRec(t, "chr1\t1000\t.\tA\tG,<NON_REF>\t.\t.\tEND=2000\tGT:LEN:DP\t0/1:100:30\t0/0:.:20\t./.::5")
	f := recordFields{rec}

	if got := f.Ref(); got != "A" {
		t.Errorf("Ref() = %q, want A", got)
	}
	if got := f.Alts(); !reflect.DeepEqual(got, []string{"G", "<NON_REF>"}) {
		t.Errorf("Alts() = %v, want [G <NON_REF>]", got)
	}
	if v, ok := f.Info("END"); !ok || v != "2000" {
		t.Errorf("Info(END) = (%q, %v), want (2000, true)", v, ok)
	}
	if _, ok := f.Info("NOPE"); ok {
		t.Error("Info reported an absent key as present")
	}

	// SampleValues gathers one FORMAT key across all samples, dropping "." and
	// "" -- an unknown length must not read as a length of nothing.
	if got := f.SampleValues("LEN"); !reflect.DeepEqual(got, []string{"100"}) {
		t.Errorf("SampleValues(LEN) = %v, want [100]", got)
	}
	if got := f.SampleValues("DP"); !reflect.DeepEqual(got, []string{"30", "20", "5"}) {
		t.Errorf("SampleValues(DP) = %v, want [30 20 5]", got)
	}
	// A key no sample carries yields an empty, non-nil slice when there are
	// samples -- distinct from the nil returned when there are none.
	empty := f.SampleValues("NOPE")
	if len(empty) != 0 {
		t.Errorf("SampleValues(NOPE) = %v, want empty", empty)
	}
	if empty == nil {
		t.Error("SampleValues returned nil for a record that has samples; want an empty slice")
	}

	bare := recordFields{editRec(t, "chr1\t1000\t.\tA\tG\t.\t.\t.")}
	if got := bare.SampleValues("DP"); got != nil {
		t.Errorf("SampleValues on a sample-less record = %v, want nil", got)
	}
	// ALT "." parses to no alternates at all, which is what makes IsRefBlock
	// treat such a record as pure reference.
	dotAlt := recordFields{editRec(t, "chr1\t1000\t.\tA\t.\t.\t.\t.")}
	if got := dotAlt.Alts(); got != nil {
		t.Errorf("Alts() for ALT \".\" = %v, want nil", got)
	}
}
