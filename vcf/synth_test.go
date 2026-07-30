package vcf

import (
	"bytes"
	"strings"
	"testing"
)

// Tests for synthesizing a VCF from something that is not a VCF: a header built
// from nothing, and records built from nothing that carry genotypes. Callers that
// need this are writing a VCF out of a genotype store, a simulator or a merger.

// TestNewRecordHasNoSamples pins the gap these primitives close. NewRecord's
// per-sample paths are all bounded by a sample count of zero, so without
// EnsureSamples a synthesized record cannot express a genotype at all.
func TestNewRecordHasNoSamples(t *testing.T) {
	r := NewRecord("chr1", 100, "A", []string{"G"})
	if n := r.NumSamples(); n != 0 {
		t.Errorf("NumSamples on a bare NewRecord = %d, want 0", n)
	}
	if _, err := r.Sample(0); err == nil {
		t.Error("Sample(0) on a bare NewRecord should fail, not return empty attributes")
	}
	if err := r.AddFormat(0, "GT", "0/1"); err == nil {
		t.Error("AddFormat on a bare NewRecord should fail")
	}
}

func TestEnsureSamplesAllowsGenotypes(t *testing.T) {
	r := NewRecordWithSamples("chr1", 100, "A", []string{"G"}, 3)
	if n := r.NumSamples(); n != 3 {
		t.Fatalf("NumSamples = %d, want 3", n)
	}
	for i, gt := range []string{"0/1", "0/0", "./."} {
		if err := r.AddFormat(i, "GT", gt); err != nil {
			t.Fatalf("AddFormat(%d): %v", i, err)
		}
	}
	got := r.serialize()
	// FILTER is PASS, not ".": NewRecord documents nil filters as PASS.
	want := "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\t./."
	if got != want {
		t.Errorf("serialize:\n got %q\nwant %q", got, want)
	}
}

// TestEnsureSamplesIsIdempotent: growing to a size already reached must not add
// columns or clear the genotypes already set.
func TestEnsureSamplesIsIdempotent(t *testing.T) {
	r := NewRecordWithSamples("chr1", 100, "A", []string{"G"}, 2)
	_ = r.AddFormat(0, "GT", "1/1")
	r.EnsureSamples(1)
	r.EnsureSamples(2)
	if n := r.NumSamples(); n != 2 {
		t.Errorf("NumSamples after shrinking calls = %d, want 2", n)
	}
	s0, err := r.Sample(0)
	if err != nil {
		t.Fatal(err)
	}
	if v, _ := s0.Get("GT"); v.String() != "1/1" {
		t.Errorf("sample 0 GT = %q, want 1/1 (EnsureSamples must not reset existing slots)", v.String())
	}
}

// TestEnsureSamplesOnParsedRecord: growing a record that came from a real line
// must keep the parsed samples intact and render the new ones as missing.
func TestEnsureSamplesOnParsedRecord(t *testing.T) {
	h := headerWithSamples(t, "S1", "S2")
	r := parseInto(t, h, "chr1\t100\t.\tA\tG\t.\t.\t.\tGT:DP\t0/1:30\t1/1:25")
	r.EnsureSamples(4)
	if n := r.NumSamples(); n != 4 {
		t.Fatalf("NumSamples = %d, want 4", n)
	}
	s1, err := r.Sample(1)
	if err != nil {
		t.Fatal(err)
	}
	if v, _ := s1.Get("DP"); v.String() != "25" {
		t.Errorf("parsed sample 1 DP = %q, want 25", v.String())
	}
	line := r.serialize()
	if !strings.HasSuffix(line, "\t0/1:30\t1/1:25\t.\t.") {
		t.Errorf("new slots should render as missing, got %q", line)
	}
}

// TestNoCallFirstSampleKeepsOthersFields is the genotype-matrix shape of the
// same bug: the first sample is the one with no call, so it has no DP -- and
// deriving the column from it used to drop every other sample's depth.
func TestNoCallFirstSampleKeepsOthersFields(t *testing.T) {
	r := NewRecordWithSamples("chr1", 100, "A", []string{"G"}, 2)
	_ = r.AddFormat(0, "GT", "./.") // no call: nothing else to report
	_ = r.AddFormat(1, "GT", "0/1")
	_ = r.AddFormat(1, "DP", "30")

	got := r.serialize()
	// Sample 0 renders as a bare "./." -- trailing missing is trimmed, though GT
	// is always kept -- while sample 1 keeps its depth. FILTER is PASS because
	// NewRecord documents nil filters as PASS.
	want := "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT:DP\t./.\t0/1:30"
	if got != want {
		t.Errorf("serialize:\n got %q\nwant %q", got, want)
	}
}

// TestSetFormatKeysForcesStableColumn covers what SetFormatKeys is still for now
// that the union handles truncation: guaranteeing a column exists even when no
// sample carries it, so a matrix writer emits the same FORMAT on every record
// rather than one that shrinks whenever a row happens to lack a field.
func TestSetFormatKeysForcesStableColumn(t *testing.T) {
	r := NewRecordWithSamples("chr1", 100, "A", []string{"G"}, 2)
	r.SetFormatKeys([]string{"GT", "DP"})
	_ = r.AddFormat(0, "GT", "./.")
	_ = r.AddFormat(1, "GT", "0/0") // nobody carries DP at this site

	got := r.serialize()
	if !strings.Contains(got, "\tGT:DP\t") {
		t.Errorf("SetFormatKeys should force DP into the column even unused, got %q", got)
	}

	// And a key some sample does carry is never silently dropped, even if the
	// forced list omits it.
	r2 := NewRecordWithSamples("chr1", 100, "A", []string{"G"}, 2)
	r2.SetFormatKeys([]string{"GT"})
	_ = r2.AddFormat(0, "GT", "0/1")
	_ = r2.AddFormat(1, "GT", "0/0")
	_ = r2.AddFormat(1, "GQ", "42")
	if got := r2.serialize(); !strings.Contains(got, "\tGT:GQ\t") {
		t.Errorf("a carried key must be appended rather than dropped, got %q", got)
	}
}

// TestDerivedFormatKeysStillDefault guards the annotate package: adding a
// per-sample key to sample 0 must still reach the FORMAT column, which only
// works while deriving is the default.
func TestDerivedFormatKeysStillDefault(t *testing.T) {
	h := headerWithSamples(t, "S1")
	r := parseInto(t, h, "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1")
	if err := r.AddFormat(0, "XX", "7"); err != nil {
		t.Fatal(err)
	}
	if got := r.serialize(); !strings.HasSuffix(got, "\tGT:XX\t0/1:7") {
		t.Errorf("an added FORMAT key must appear without SetFormatKeys, got %q", got)
	}
}

// TestFormatKeyOnLaterSampleSurvives is the regression test for the bug that
// deriving the FORMAT column from sample 0 caused.
//
// bigwig, bigbed and tabix annotators all write their field to a single
// --sample. Whenever that was not the first sample, the field never reached the
// output at all: sample 0 did not have the key, so the key was not in FORMAT, so
// every sample rendered without it.
func TestFormatKeyOnLaterSampleSurvives(t *testing.T) {
	h := headerWithSamples(t, "S1", "S2")
	r := parseInto(t, h, "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/0")
	if err := r.AddFormat(1, "CG_BW", "0.75"); err != nil {
		t.Fatal(err)
	}
	got := r.serialize()
	want := "chr1\t100\t.\tA\tG\t.\t.\t.\tGT:CG_BW\t0/1\t0/0:0.75"
	if got != want {
		t.Errorf("a FORMAT field written to sample 1 must reach the output:\n got %q\nwant %q", got, want)
	}
}

// TestFormatKeyRemovedFromEverySampleIsDropped pins the other direction, which
// cgkit's vcf-strip depends on: dropping a field means removing it from every
// sample, and the column must then lose it.
func TestFormatKeyRemovedFromEverySampleIsDropped(t *testing.T) {
	h := headerWithSamples(t, "S1", "S2")
	r := parseInto(t, h, "chr1\t100\t.\tA\tG\t.\t.\t.\tGT:DP:GQ\t0/1:30:99\t0/0:25:80")
	for i := 0; i < r.NumSamples(); i++ {
		s, err := r.Sample(i)
		if err != nil {
			t.Fatal(err)
		}
		s.Remove("DP")
	}
	r.MarkDirty()
	got := r.serialize()
	want := "chr1\t100\t.\tA\tG\t.\t.\t.\tGT:GQ\t0/1:99\t0/0:80"
	if got != want {
		t.Errorf("removing a key from every sample must drop it from FORMAT:\n got %q\nwant %q", got, want)
	}
}

// TestFormatColumnOrderIsStable: the record's own FORMAT order is preserved and
// added keys append, so a fix to the key set does not reshuffle columns.
func TestFormatColumnOrderIsStable(t *testing.T) {
	h := headerWithSamples(t, "S1", "S2")
	r := parseInto(t, h, "chr1\t100\t.\tA\tG\t.\t.\t.\tGT:DP:GQ\t0/1:30:99\t0/0:25:80")
	if err := r.AddFormat(1, "ZZ", "1"); err != nil {
		t.Fatal(err)
	}
	if got := r.serialize(); !strings.Contains(got, "\tGT:DP:GQ:ZZ\t") {
		t.Errorf("FORMAT order should stay GT:DP:GQ then append ZZ, got %q", got)
	}
}

// TestNewVcfHeaderWritesValidHeader: a header built from nothing must produce a
// usable ##fileformat line and #CHROM line.
func TestNewVcfHeaderWritesValidHeader(t *testing.T) {
	h := NewVcfHeader()
	h.SetSamples([]string{"S1", "S2"})
	h.AddFormat(&AnnotationDef{ID: "GT", Number: "1", Type: "String", Description: "Genotype"})

	var buf bytes.Buffer
	if _, err := h.WriteTo(&buf); err != nil {
		t.Fatalf("WriteTo: %v", err)
	}
	out := buf.String()
	for _, want := range []string{
		"##fileformat=VCFv4.2",
		`##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`,
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
	} {
		if !strings.Contains(out, want) {
			t.Errorf("header missing %q; got:\n%s", want, out)
		}
	}
	if got := h.SampleIndex("S2"); got != 1 {
		t.Errorf("SampleIndex(S2) = %d, want 1", got)
	}
}

// TestSynthesizedVcfRoundTrips writes a synthesized VCF through the real writer
// and reads it back, which is the end-to-end guard on all of the above.
func TestSynthesizedVcfRoundTrips(t *testing.T) {
	samples := []string{"S1", "S2", "S3"}
	h := NewVcfHeader()
	h.SetSamples(samples)
	h.AddFormat(&AnnotationDef{ID: "GT", Number: "1", Type: "String", Description: "Genotype"})
	h.AddFormat(&AnnotationDef{ID: "DP", Number: "1", Type: "Integer", Description: "Read depth"})
	h.AddContig(&ContigDef{ID: "chr1", Length: 1000})

	var buf bytes.Buffer
	w := NewVcfWriter(&buf)
	if err := w.WriteHeader(h); err != nil {
		t.Fatal(err)
	}
	gts := [][]string{{"0/1", "0/0", "./."}, {"1/1", "./.", "0/0"}}
	for i, pos := range []int{100, 200} {
		r := NewRecordWithSamples("chr1", pos, "A", []string{"G"}, len(samples))
		r.SetFormatKeys([]string{"GT", "DP"})
		for j, gt := range gts[i] {
			if err := r.AddFormat(j, "GT", gt); err != nil {
				t.Fatal(err)
			}
			if gt != "./." {
				if err := r.AddFormat(j, "DP", "30"); err != nil {
					t.Fatal(err)
				}
			}
		}
		if err := w.WriteRecord(r); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	rd, err := NewVcfReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("reading back: %v", err)
	}
	rh, err := rd.Header()
	if err != nil {
		t.Fatal(err)
	}
	if got := rh.Samples(); len(got) != 3 || got[2] != "S3" {
		t.Fatalf("round-tripped samples = %v, want %v", got, samples)
	}
	for i := 0; i < 2; i++ {
		rec, err := rd.NextRecord()
		if err != nil {
			t.Fatalf("record %d: %v", i, err)
		}
		for j, want := range gts[i] {
			attrs, err := rec.Sample(j)
			if err != nil {
				t.Fatalf("record %d sample %d: %v", i, j, err)
			}
			v, ok := attrs.Get("GT")
			if !ok || v.String() != want {
				t.Errorf("record %d sample %d GT = %q, want %q", i, j, v.String(), want)
			}
		}
	}
}

// headerWithSamples builds a parsed header carrying the given samples.
func headerWithSamples(t *testing.T, samples ...string) *VcfHeader {
	t.Helper()
	lines := "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
		strings.Join(samples, "\t") + "\n"
	r, err := NewVcfReader(strings.NewReader(lines))
	if err != nil {
		t.Fatalf("building header: %v", err)
	}
	h, err := r.Header()
	if err != nil {
		t.Fatalf("building header: %v", err)
	}
	return h
}

// parseInto parses one data line against a header, going through the reader so
// the record is built exactly as a real one would be.
func parseInto(t *testing.T, h *VcfHeader, line string) *VcfRecord {
	t.Helper()
	var b strings.Builder
	for _, l := range h.Lines() {
		b.WriteString(l + "\n")
	}
	b.WriteString(h.ChromLine() + "\n")
	b.WriteString(line + "\n")
	r, err := NewVcfReader(strings.NewReader(b.String()))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := r.Header(); err != nil {
		t.Fatal(err)
	}
	rec, err := r.NextRecord()
	if err != nil {
		t.Fatalf("parsing %q: %v", line, err)
	}
	return rec
}
