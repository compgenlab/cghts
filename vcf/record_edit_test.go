package vcf

import (
	"reflect"
	"strings"
	"testing"
)

// The record mutators all share one hazard: they set the dirty flag, which
// switches the writer from emitting the source line verbatim to reconstructing
// it. So each of these asserts on the *serialized* result, not just on the
// accessor -- an edit the accessor reports but serialize drops is exactly the
// failure mode that ships wrong output.

const editLine = "chr1\t100\trs1\tA\tG\t50.0\tPASS\tDP=30;AF=0.5\tGT:AD\t0/0:28,2\t0/1:15,15"

func editRec(t *testing.T, line string) *VcfRecord {
	t.Helper()
	rec, err := newRecord(line, nil)
	if err != nil {
		t.Fatalf("newRecord(%q): %v", line, err)
	}
	return rec
}

// filterCol returns the FILTER column of the serialized record.
func filterCol(t *testing.T, rec *VcfRecord) string {
	t.Helper()
	cols := strings.Split(strings.TrimRight(writeRec(t, rec), "\n"), "\t")
	if len(cols) < 7 {
		t.Fatalf("serialized record has %d columns: %v", len(cols), cols)
	}
	return cols[6]
}

// FILTER has three states, not two, and the mutators have to preserve the
// difference: no codes and PASS mean "assessed, passed", while "." means "not
// assessed at all". Collapsing them turns an unfiltered call into a passing one.
func TestRecordFilterMutators(t *testing.T) {
	for _, tc := range []struct {
		name string
		src  string // FILTER column of the source line
		edit func(*VcfRecord)
		want string
	}{
		{"AddFilter to a PASS record", "PASS", func(r *VcfRecord) { r.AddFilter("lowqual") }, "lowqual"},
		{"AddFilter accumulates", "lowqual", func(r *VcfRecord) { r.AddFilter("lowdp") }, "lowqual;lowdp"},
		{"AddFilter to an unassessed record", ".", func(r *VcfRecord) { r.AddFilter("lowqual") }, "lowqual"},
		{"ClearFilters returns to PASS", "lowqual;lowdp", func(r *VcfRecord) { r.ClearFilters() }, "PASS"},
		{"ClearFilters on PASS is a no-op", "PASS", func(r *VcfRecord) { r.ClearFilters() }, "PASS"},
		// ClearFilters is stronger than RetainFilters-to-empty: it discards the
		// "was assessed" fact along with the codes.
		{"ClearFilters erases the unassessed marker", ".", func(r *VcfRecord) { r.ClearFilters() }, "PASS"},
		{"SetFilters replaces", "lowqual", func(r *VcfRecord) { r.SetFilters([]string{"a", "b"}) }, "a;b"},
		{"SetFilters on PASS", "PASS", func(r *VcfRecord) { r.SetFilters([]string{"a"}) }, "a"},
		{"SetFilters with an empty slice clears to PASS", "lowqual",
			func(r *VcfRecord) { r.SetFilters([]string{}) }, "PASS"},
		{"SetFilters with nil clears to PASS", "lowqual",
			func(r *VcfRecord) { r.SetFilters(nil) }, "PASS"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			line := "chr1\t100\trs1\tA\tG\t50.0\t" + tc.src + "\tDP=30"
			rec := editRec(t, line)
			tc.edit(rec)
			if !rec.Dirty() {
				t.Error("the edit did not mark the record dirty, so a writer would emit the source line unchanged")
			}
			if got := filterCol(t, rec); got != tc.want {
				t.Errorf("FILTER = %q, want %q", got, tc.want)
			}
		})
	}
}

// SetFilters copies rather than aliasing: a caller reusing its slice must not be
// able to rewrite a record it already handed off.
func TestSetFiltersCopies(t *testing.T) {
	rec := editRec(t, editLine)
	codes := []string{"lowqual"}
	rec.SetFilters(codes)
	codes[0] = "MUTATED"
	if got := filterCol(t, rec); got != "lowqual" {
		t.Errorf("FILTER = %q, want lowqual -- SetFilters aliased the caller's slice", got)
	}
}

// Filters distinguishes the three FILTER states by nil-ness, which is subtle
// enough that every caller gets it wrong at least once.
func TestFiltersNilness(t *testing.T) {
	for _, tc := range []struct {
		src     string
		wantNil bool
		want    []string
	}{
		{"PASS", true, nil},
		{".", false, []string{}},
		{"lowqual", false, []string{"lowqual"}},
		{"lowqual;lowdp", false, []string{"lowqual", "lowdp"}},
		// A "." mixed in with real codes is dropped, not kept as a code.
		{"lowqual;.", false, []string{"lowqual"}},
	} {
		rec := editRec(t, "chr1\t100\t.\tA\tG\t50.0\t"+tc.src+"\t.")
		got := rec.Filters()
		if (got == nil) != tc.wantNil {
			t.Errorf("FILTER %q: Filters() nil = %v, want %v", tc.src, got == nil, tc.wantNil)
		}
		if !reflect.DeepEqual(got, tc.want) {
			t.Errorf("FILTER %q: Filters() = %#v, want %#v", tc.src, got, tc.want)
		}
		if want := len(tc.want) > 0; rec.IsFiltered() != want {
			t.Errorf("FILTER %q: IsFiltered() = %v, want %v", tc.src, rec.IsFiltered(), want)
		}
	}
}

func TestClearID(t *testing.T) {
	rec := editRec(t, editLine)
	if rec.ID() != "rs1" {
		t.Fatalf("ID() = %q, want rs1", rec.ID())
	}
	rec.ClearID()
	if rec.ID() != "" {
		t.Errorf("after ClearID, ID() = %q, want empty", rec.ID())
	}
	// The ID column renders as ".", not as an empty column -- an empty column
	// would shift every field after it.
	cols := strings.Split(strings.TrimRight(writeRec(t, rec), "\n"), "\t")
	if cols[2] != "." {
		t.Errorf("serialized ID column = %q, want %q", cols[2], ".")
	}
}

// AltOrig is the raw column; Alt() is the parsed alleles with "." dropped. The
// pair exists because a caller writing a record out needs the source spelling
// while a caller counting alleles needs the parsed list.
func TestAltOrigVsAlt(t *testing.T) {
	for _, tc := range []struct {
		alt      string
		wantOrig string
		wantAlt  []string
	}{
		{"G", "G", []string{"G"}},
		{"G,T", "G,T", []string{"G", "T"}},
		{".", ".", nil},
		// A "." among real alleles is dropped from Alt() but kept verbatim in
		// AltOrig, so the two lists differ in length.
		{"G,.", "G,.", []string{"G"}},
		{"<NON_REF>", "<NON_REF>", []string{"<NON_REF>"}},
		{"T[chr5:2000[", "T[chr5:2000[", []string{"T[chr5:2000["}},
	} {
		rec := editRec(t, "chr1\t100\t.\tA\t"+tc.alt+"\t50.0\tPASS\t.")
		if got := rec.AltOrig(); got != tc.wantOrig {
			t.Errorf("AltOrig() for %q = %q, want %q", tc.alt, got, tc.wantOrig)
		}
		if got := rec.Alt(); !reflect.DeepEqual(got, tc.wantAlt) {
			t.Errorf("Alt() for %q = %v, want %v", tc.alt, got, tc.wantAlt)
		}
	}
}

func TestInfoValue(t *testing.T) {
	rec := editRec(t, "chr1\t100\t.\tA\tG\t50.0\tPASS\tDP=30;DB;MISS=.")

	if v, ok := rec.InfoValue("DP"); !ok || v.String() != "30" {
		t.Errorf("InfoValue(DP) = (%q, %v), want (30, true)", v.String(), ok)
	}
	// A bare flag is present with an empty value -- absence and a valueless
	// presence are different answers.
	if v, ok := rec.InfoValue("DB"); !ok || !v.IsEmpty() {
		t.Errorf("InfoValue(DB) = (%q, %v), want an empty value, true", v.String(), ok)
	}
	if v, ok := rec.InfoValue("MISS"); !ok || !v.IsMissing() {
		t.Errorf("InfoValue(MISS) = (%q, %v), want the missing marker, true", v.String(), ok)
	}
	if _, ok := rec.InfoValue("NOPE"); ok {
		t.Error("InfoValue reported an absent key as present")
	}
}

func TestFormatKeysAndZeroBasedStart(t *testing.T) {
	rec := editRec(t, editLine)
	if got := rec.FormatKeys(); !reflect.DeepEqual(got, []string{"GT", "AD"}) {
		t.Errorf("FormatKeys() = %v, want [GT AD]", got)
	}
	// POS is 1-based; BED output is 0-based, and off-by-one here shifts every
	// exported interval.
	if got := rec.ZeroBasedStart(); got != 99 {
		t.Errorf("ZeroBasedStart() = %d, want 99 (POS 100)", got)
	}

	// A record with no sample columns has no FORMAT keys at all, rather than one
	// empty key from splitting "".
	bare := editRec(t, "chr1\t100\t.\tA\tG\t50.0\tPASS\tDP=30")
	if got := bare.FormatKeys(); len(got) != 0 {
		t.Errorf("FormatKeys() on a sample-less record = %v, want empty", got)
	}
	if got := bare.NumSamples(); got != 0 {
		t.Errorf("NumSamples() on a sample-less record = %d, want 0", got)
	}
}

func TestIsIndel(t *testing.T) {
	for _, tc := range []struct {
		ref, alt string
		want     bool
	}{
		{"A", "G", false},
		{"A", "G,T", false},
		{"G", "GA", true},  // insertion
		{"GA", "G", true},  // deletion
		{"AC", "GT", true}, // MNV counts, since neither allele is one base
		// One long allele among short ones is enough.
		{"A", "G,GTT", true},
		// A symbolic alternate is longer than one base, so it reads as an indel.
		{"N", "<DEL>", true},
		{"A", ".", false},
	} {
		rec := editRec(t, "chr1\t100\t.\t"+tc.ref+"\t"+tc.alt+"\t50.0\tPASS\t.")
		if got := rec.IsIndel(); got != tc.want {
			t.Errorf("IsIndel() for %s>%s = %v, want %v", tc.ref, tc.alt, got, tc.want)
		}
	}
}

// Line is the source line, and stays the source line after an edit. That is what
// makes it usable as provenance; a caller wanting the edited form asks the writer.
func TestLineIsTheSourceLine(t *testing.T) {
	rec := editRec(t, editLine)
	if got := rec.Line(); got != editLine {
		t.Errorf("Line() = %q, want the source line", got)
	}
	rec.AddInfo("NEW", "1")
	if got := rec.Line(); got != editLine {
		t.Errorf("Line() changed after an edit: %q -- it must stay the source line", got)
	}
	if !strings.Contains(writeRec(t, rec), "NEW=1") {
		t.Error("the edit did not reach the serialized output")
	}
}

// A trailing newline in the source is stripped, so Line() concatenates cleanly.
func TestLineStripsLineEndings(t *testing.T) {
	for _, suffix := range []string{"\n", "\r\n", ""} {
		rec := editRec(t, editLine+suffix)
		if got := rec.Line(); got != editLine {
			t.Errorf("Line() with suffix %q = %q, want the bare line", suffix, got)
		}
	}
}

func TestReorderSamplesLine(t *testing.T) {
	rec := editRec(t, editLine)
	for _, tc := range []struct {
		name  string
		order []int
		want  string
	}{
		{"swap", []int{1, 0}, "chr1\t100\trs1\tA\tG\t50.0\tPASS\tDP=30;AF=0.5\tGT:AD\t0/1:15,15\t0/0:28,2"},
		{"identity", []int{0, 1}, editLine},
		{"subset", []int{1}, "chr1\t100\trs1\tA\tG\t50.0\tPASS\tDP=30;AF=0.5\tGT:AD\t0/1:15,15"},
		{"duplicate", []int{0, 0}, "chr1\t100\trs1\tA\tG\t50.0\tPASS\tDP=30;AF=0.5\tGT:AD\t0/0:28,2\t0/0:28,2"},
		// An out-of-range index emits "." rather than failing, which is what lets
		// a caller widen a cohort with placeholder columns.
		{"out of range", []int{0, 7}, "chr1\t100\trs1\tA\tG\t50.0\tPASS\tDP=30;AF=0.5\tGT:AD\t0/0:28,2\t."},
		{"negative", []int{-1, 0}, "chr1\t100\trs1\tA\tG\t50.0\tPASS\tDP=30;AF=0.5\tGT:AD\t.\t0/0:28,2"},
		// An empty order drops every sample column but keeps FORMAT, which is a
		// malformed record -- worth pinning so a caller knows to pass nil samples
		// through a different route.
		{"empty order keeps FORMAT", []int{}, "chr1\t100\trs1\tA\tG\t50.0\tPASS\tDP=30;AF=0.5\tGT:AD"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			if got := rec.ReorderSamplesLine(tc.order); got != tc.want {
				t.Errorf("ReorderSamplesLine(%v):\n got: %q\nwant: %q", tc.order, got, tc.want)
			}
		})
	}
}

// A record with no sample columns has nothing to permute, so the line comes back
// whole -- rather than gaining a stray FORMAT column.
func TestReorderSamplesLineWithoutSamples(t *testing.T) {
	const line = "chr1\t100\trs1\tA\tG\t50.0\tPASS\tDP=30"
	rec := editRec(t, line)
	if got := rec.ReorderSamplesLine([]int{1, 0}); got != line {
		t.Errorf("ReorderSamplesLine on a sample-less record = %q, want the line unchanged", got)
	}
}

// ReorderSamplesLine works from the raw line and does not consult the parsed
// model, so edits made to the record are absent from its result. That is
// deliberate -- it is what makes permuting a 3,000-sample cohort cheap -- but it
// means a caller must not combine it with record mutation. cgkit's vcf-reorder
// passes a no-op record hook for exactly this reason.
func TestReorderSamplesLineIgnoresEdits(t *testing.T) {
	rec := editRec(t, editLine)
	rec.AddInfo("NEW", "1")
	rec.ClearID()

	got := rec.ReorderSamplesLine([]int{0, 1})
	if strings.Contains(got, "NEW=1") {
		t.Error("ReorderSamplesLine picked up an INFO edit; the raw-line contract has changed")
	}
	if got != editLine {
		t.Errorf("ReorderSamplesLine = %q, want the unedited source line", got)
	}
}
