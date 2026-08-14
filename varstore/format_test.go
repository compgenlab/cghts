package varstore

import (
	"path/filepath"
	"testing"
)

// Captured FORMAT fields.
//
// The same shape as captured INFO, with the costs the other way round: a call
// row is one sample at one ALT, so the cardinality fits better, but calls is
// the large member so each column is expensive. The tests that matter are the
// compatibility ones -- an existing reader must not notice the extra columns --
// and the refusals, since a field that cannot be stored must say so rather than
// arriving as a column of zeros.

func formatWriter(t *testing.T, dir string, fields []FormatField) *Writer {
	t.Helper()
	w, err := NewWriter(dir, WriterOpts{Samples: []string{"S1", "S2"}, MinDP: 10, Format: fields})
	if err != nil {
		t.Fatalf("NewWriter: %v", err)
	}
	return w
}

var pidField = FormatField{Name: "PID", Type: InfoString, Number: "1"}
var vafField = FormatField{Name: "VAF", Type: InfoFloat, Number: "A"}

func TestFormatRoundTrips(t *testing.T) {
	dir := filepath.Join(t.TempDir(), "store")
	w := formatWriter(t, dir, []FormatField{pidField, vafField})

	if err := w.WriteSite(Site{Chrom: "chr1", Pos: 100, Ref: "G", Alt: "A", AN: 4, NCalled: 2}); err != nil {
		t.Fatal(err)
	}
	if err := w.WriteCallFormat(Call{
		SampleID: "S1", Chrom: "chr1", Pos: 100, Ref: "G", Alt: "A",
		GT: "0/1", DP: 40, ADRef: 20, ADAlt: 20, GQ: 99,
	}, map[string]any{"PID": "100_G_A", "VAF": 0.51}); err != nil {
		t.Fatal(err)
	}
	// S2 carries the variant but the source published no VAF for it.
	if err := w.WriteCallFormat(Call{
		SampleID: "S2", Chrom: "chr1", Pos: 100, Ref: "G", Alt: "A",
		GT: "1/1", DP: 30, ADRef: 0, ADAlt: 30, GQ: 99,
	}, map[string]any{"PID": "100_G_A"}); err != nil {
		t.Fatal(err)
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}

	s, err := OpenParquet(dir)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	if got := s.FormatFields(); len(got) != 2 {
		t.Fatalf("manifest records %d format fields, want 2: %+v", len(got), got)
	}

	seen := map[string]struct {
		call Call
		vaf  *float64
		pid  string
	}{}
	if err := s.CallsFormat(func(c Call, f FormatRow) bool {
		e := seen[c.SampleID]
		e.call = c
		e.pid, _ = f.String("fmt_pid")
		if v, ok := f.Float("fmt_vaf"); ok {
			e.vaf = &v
		}
		seen[c.SampleID] = e
		return true
	}); err != nil {
		t.Fatal(err)
	}
	if len(seen) != 2 {
		t.Fatalf("read %d calls, want 2", len(seen))
	}
	// The Call itself must survive the dynamic read intact.
	if c := seen["S1"].call; c.GT != "0/1" || c.DP != 40 || c.ADAlt != 20 || c.GQ != 99 {
		t.Errorf("S1 rebuilt wrong: %+v", c)
	}
	if seen["S1"].pid != "100_G_A" {
		t.Errorf("S1 PID = %q", seen["S1"].pid)
	}
	if seen["S1"].vaf == nil || *seen["S1"].vaf != 0.51 {
		t.Errorf("S1 VAF = %v, want 0.51", seen["S1"].vaf)
	}
	// ABSENT IS NOT ZERO. S2 had no VAF, and a 0.0 there would be a claim the
	// source never made.
	if seen["S2"].vaf != nil {
		t.Errorf("S2 published no VAF but it reads as %v", *seen["S2"].vaf)
	}
}

// An existing reader must not notice that calls grew columns. This is the whole
// premise of putting them in calls.parquet rather than a fourth member.
func TestTypedCallReaderIgnoresCapturedColumns(t *testing.T) {
	dir := filepath.Join(t.TempDir(), "store")
	w := formatWriter(t, dir, []FormatField{pidField})
	if err := w.WriteSite(Site{Chrom: "chr1", Pos: 100, Ref: "G", Alt: "A"}); err != nil {
		t.Fatal(err)
	}
	if err := w.WriteCallFormat(Call{
		SampleID: "S1", Chrom: "chr1", Pos: 100, Ref: "G", Alt: "A",
		GT: "0/1", DP: 40, ADRef: 20, ADAlt: 20, GQ: 50,
	}, map[string]any{"PID": "x"}); err != nil {
		t.Fatal(err)
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}

	s, err := OpenParquet(dir)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	seq, err := s.Calls(Query{})
	if err != nil {
		t.Fatal(err)
	}
	var n int
	for c, err := range seq {
		if err != nil {
			t.Fatal(err)
		}
		n++
		if c.SampleID != "S1" || c.GT != "0/1" || c.DP != 40 {
			t.Errorf("typed read got %+v", c)
		}
	}
	if n != 1 {
		t.Fatalf("typed read saw %d calls, want 1", n)
	}
}

func TestValidateFormatRefusals(t *testing.T) {
	for _, c := range []struct {
		name  string
		field FormatField
		why   string
	}{
		{"per-allele-with-ref", FormatField{Name: "AD", Type: InfoInteger, Number: "R"},
			"Number=R has no room in a call row, and AD is already stored"},
		{"per-genotype", FormatField{Name: "PL", Type: InfoInteger, Number: "G"},
			"Number=G is three values for a biallelic site"},
		{"variable", FormatField{Name: "SAC", Type: InfoInteger, Number: "."},
			"Number=. varies per record"},
		{"flag", FormatField{Name: "X", Type: InfoFlag, Number: "0"},
			"a per-sample field with no value cannot be told from an absent one"},
		{"reserved-column", FormatField{Name: "DP", Column: "dp", Type: InfoInteger, Number: "1"},
			"dp is the call's own"},
	} {
		t.Run(c.name, func(t *testing.T) {
			if err := ValidateFormat([]FormatField{c.field}); err == nil {
				t.Fatalf("accepted %+v, but %s", c.field, c.why)
			}
		})
	}
	// The prefix is what makes capturing DP safe rather than forbidden.
	if err := ValidateFormat([]FormatField{{Name: "DP", Type: InfoInteger, Number: "1"}}); err != nil {
		t.Errorf("DP should be capturable as fmt_dp: %v", err)
	}
	if got := FormatColumn("DP"); got != "fmt_dp" {
		t.Errorf("FormatColumn(DP) = %q, want fmt_dp", got)
	}
}

func TestWritingFormatToAStoreThatCapturesNoneFails(t *testing.T) {
	dir := filepath.Join(t.TempDir(), "store")
	w := formatWriter(t, dir, nil)
	err := w.WriteCallFormat(Call{
		SampleID: "S1", Chrom: "chr1", Pos: 1, Ref: "A", Alt: "T", GT: "0/1",
	}, map[string]any{"PID": "x"})
	if err == nil {
		t.Fatal("accepted FORMAT values into a store with no captured columns")
	}
	_ = w.Close()
}
