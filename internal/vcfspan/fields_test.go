package vcfspan

import (
	"reflect"
	"strings"
	"testing"
)

// LineFields is the adapter the tabix index writer reads a record through,
// working on raw split columns rather than a parsed record. It had no direct
// test, so the accessors below were covered only incidentally through End.

func fields(line string) LineFields {
	return LineFields(strings.Split(line, "\t"))
}

const spanLine = "chr1\t1000\t.\tACGT\tA,<DEL>\t50\tPASS\tEND=2000;SVLEN=-1000;DB\tGT:DP:LEN\t0/1:30:900\t1/1:.:800"

func TestLineFieldsAccessors(t *testing.T) {
	f := fields(spanLine)

	if got := f.Ref(); got != "ACGT" {
		t.Errorf("Ref() = %q, want %q", got, "ACGT")
	}
	if got := f.Alts(); !reflect.DeepEqual(got, []string{"A", "<DEL>"}) {
		t.Errorf("Alts() = %v, want [A <DEL>]", got)
	}

	if v, ok := f.Info("END"); !ok || v != "2000" {
		t.Errorf("Info(END) = (%q, %v), want (2000, true)", v, ok)
	}
	// A key that is a prefix of another must not match it.
	if v, ok := f.Info("EN"); ok {
		t.Errorf("Info(EN) matched %q; keys must match whole, not by prefix", v)
	}
	// SVLEN sits after END, so this also covers walking past the first pair.
	if v, ok := f.Info("SVLEN"); !ok || v != "-1000" {
		t.Errorf("Info(SVLEN) = (%q, %v), want (-1000, true)", v, ok)
	}
	if _, ok := f.Info("NOPE"); ok {
		t.Error("Info(NOPE) reported present")
	}
	// A flag has no "=", so this accessor cannot see it. That is correct for
	// every key it is actually used for -- END, SVLEN and LEN are all typed --
	// but it is the reason a flag key would silently read as absent.
	if _, ok := f.Info("DB"); ok {
		t.Error("Info(DB) matched a valueless flag; the accessor only reads KEY=VALUE")
	}
}

func TestLineFieldsSampleValues(t *testing.T) {
	f := fields(spanLine)

	if got := f.SampleValues("LEN"); !reflect.DeepEqual(got, []string{"900", "800"}) {
		t.Errorf("SampleValues(LEN) = %v, want [900 800]", got)
	}
	// "." is a recorded value here and reaches the caller; End then fails to
	// parse it and moves on. The record-backed adapter in vcf/span.go drops it
	// instead, which is a difference between the two Fields implementations
	// rather than a difference in outcome.
	if got := f.SampleValues("DP"); !reflect.DeepEqual(got, []string{"30", "."}) {
		t.Errorf("SampleValues(DP) = %v, want [30 .]", got)
	}
	if got := f.SampleValues("NOPE"); got != nil {
		t.Errorf("SampleValues(NOPE) = %v, want nil", got)
	}
}

// A truncated sample column is legal where the trailing subfields are absent,
// and must not read as a value or panic.
func TestLineFieldsTruncatedSampleColumn(t *testing.T) {
	f := fields("chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT:DP:LEN\t0/1:30\t0/1")
	if got := f.SampleValues("LEN"); len(got) != 0 {
		t.Errorf("SampleValues(LEN) = %v, want no values for columns with no LEN subfield", got)
	}
	if got := f.SampleValues("DP"); !reflect.DeepEqual(got, []string{"30"}) {
		t.Errorf("SampleValues(DP) = %v, want [30]", got)
	}
}

// Short lines must be answered, not panicked on: the tabix writer sees whatever
// the input file contains, including malformed rows.
func TestLineFieldsShortLines(t *testing.T) {
	for _, line := range []string{"", "chr1", "chr1\t100", "chr1\t100\t.", "chr1\t100\t.\tA"} {
		f := fields(line)
		f.Ref()
		f.Alts()
		f.Info("END")
		f.SampleValues("LEN")
		if got := FieldsEnd(strings.Split(line, "\t"), 100); got < 100 {
			t.Errorf("FieldsEnd(%q) = %d, want an end at or past the start", line, got)
		}
	}
}

// FieldsEnd is the entry point the index writer uses; a record with no REF
// column at all still has to yield a one-base span rather than a zero or a panic.
func TestFieldsEnd(t *testing.T) {
	const beg = 999 // 0-based for POS=1000
	// END=2000 and SVLEN=-1000 describe the same 1001-base extent, so they
	// reconcile to one answer rather than one overriding the other.
	if got := FieldsEnd(strings.Split(spanLine, "\t"), beg); got != 2000 {
		t.Errorf("FieldsEnd = %d, want 2000 (1001 bases from %d)", got, beg)
	}
	if got := FieldsEnd([]string{"chr1", "100"}, 99); got != 100 {
		t.Errorf("FieldsEnd on a REF-less line = %d, want a one-base span", got)
	}
}
