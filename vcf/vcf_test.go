package vcf

import (
	"io"
	"os"
	"testing"
)

func openTestFile(t *testing.T) *VcfReader {
	t.Helper()
	r, err := NewVcfFile("testdata/sample.vcf")
	if err != nil {
		t.Fatalf("NewVcfFile: %v", err)
	}
	return r
}

func readAll(t *testing.T, r *VcfReader) []*VcfRecord {
	t.Helper()
	var recs []*VcfRecord
	for {
		rec, err := r.NextRecord()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("NextRecord: %v", err)
		}
		recs = append(recs, rec)
	}
	return recs
}

func TestHeaderParse(t *testing.T) {
	r := openTestFile(t)
	defer r.Close()
	h, err := r.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if got := h.Samples(); len(got) != 2 || got[0] != "NORMAL" || got[1] != "TUMOR" {
		t.Errorf("Samples = %v, want [NORMAL TUMOR]", got)
	}
	if h.SampleIndex("TUMOR") != 1 {
		t.Errorf("SampleIndex(TUMOR) = %d, want 1", h.SampleIndex("TUMOR"))
	}
	if h.SampleIndex("2") != 1 {
		t.Errorf("SampleIndex(\"2\") = %d, want 1", h.SampleIndex("2"))
	}
	if d, ok := h.InfoDef("AF"); !ok || d.Number != "A" || d.Type != "Float" {
		t.Errorf("InfoDef(AF) = %+v, ok=%v", d, ok)
	}
	if names := h.ContigNames(); len(names) != 2 || names[0] != "chr1" {
		t.Errorf("ContigNames = %v", names)
	}
}

func TestRecordEagerFields(t *testing.T) {
	r := openTestFile(t)
	defer r.Close()
	recs := readAll(t, r)
	if len(recs) != 5 {
		t.Fatalf("got %d records, want 5", len(recs))
	}
	rec := recs[0]
	if rec.Chrom != "chr1" || rec.Pos != 100 || rec.Ref != "A" {
		t.Errorf("eager fields = %s:%d %s", rec.Chrom, rec.Pos, rec.Ref)
	}
	if rec.ID() != "rs1" {
		t.Errorf("ID = %q, want rs1", rec.ID())
	}
	if alt := rec.Alt(); len(alt) != 1 || alt[0] != "G" {
		t.Errorf("Alt = %v", alt)
	}
	if rec.Qual() != 50.0 {
		t.Errorf("Qual = %v", rec.Qual())
	}
	if rec.IsFiltered() {
		t.Errorf("rec[0] should be PASS")
	}
	if recs[1].ID() != "" {
		t.Errorf("rec[1] ID should be empty (.)")
	}
	if !recs[1].IsFiltered() {
		t.Errorf("rec[1] should be filtered (lowqual)")
	}
}

func TestInfoLazy(t *testing.T) {
	r := openTestFile(t)
	defer r.Close()
	rec, _ := r.NextRecord()
	info := rec.Info()
	if v, ok := info.Get("DP"); !ok {
		t.Fatal("DP missing")
	} else if n, _ := v.Int(); n != 30 {
		t.Errorf("DP = %d, want 30", n)
	}
	if v, ok := info.Get("AF"); !ok {
		t.Fatal("AF missing")
	} else if f, _ := v.Float(); f != 0.5 {
		t.Errorf("AF = %v, want 0.5", f)
	}
	if v, ok := info.Get("DB"); !ok || !v.IsEmpty() {
		t.Errorf("DB flag = %+v, ok=%v", v, ok)
	}
	if _, ok := info.Get("NOPE"); ok {
		t.Errorf("absent key should report ok=false")
	}
}

func TestSampleLazy(t *testing.T) {
	r := openTestFile(t)
	defer r.Close()
	rec, _ := r.NextRecord()
	if rec.NumSamples() != 2 {
		t.Fatalf("NumSamples = %d", rec.NumSamples())
	}
	tumor, err := rec.SampleByName("TUMOR")
	if err != nil {
		t.Fatalf("SampleByName: %v", err)
	}
	ad, ok := tumor.Get("AD")
	if !ok {
		t.Fatal("AD missing")
	}
	if s, _ := ad.StringFor("alt1"); s != "15" {
		t.Errorf("AD alt1 = %q, want 15", s)
	}
	if sum, _ := ad.FloatFor("sum"); sum != 30 {
		t.Errorf("AD sum = %v, want 30", sum)
	}
}

func TestCalcTsTv(t *testing.T) {
	r := openTestFile(t)
	defer r.Close()
	recs := readAll(t, r)
	if recs[0].CalcTsTv() != -1 { // A>G transition
		t.Errorf("rec0 A>G TsTv = %d, want -1", recs[0].CalcTsTv())
	}
	if recs[1].CalcTsTv() != -1 { // C>T transition
		t.Errorf("rec1 C>T TsTv = %d, want -1", recs[1].CalcTsTv())
	}
	if recs[2].CalcTsTv() != 0 { // indel G>GA
		t.Errorf("rec2 indel TsTv = %d, want 0", recs[2].CalcTsTv())
	}
}

func TestAltPositions(t *testing.T) {
	r := openTestFile(t)
	defer r.Close()
	recs := readAll(t, r)

	// SNV
	ap := recs[0].AltPositions("", "", "", "")
	if len(ap) != 1 || ap[0].Type != VarSNV || ap[0].Pos != 100 {
		t.Errorf("SNV AltPositions = %+v", ap)
	}
	// <DEL> with END=900
	ap = recs[3].AltPositions("", "", "", "")
	if len(ap) != 1 || ap[0].Type != VarDEL || ap[0].Pos != 900 {
		t.Errorf("DEL AltPositions = %+v", ap)
	}
	// BND T[chr5:2000[
	ap = recs[4].AltPositions("", "", "", "")
	if len(ap) != 1 || ap[0].Type != VarBND || ap[0].Chrom != "chr5" || ap[0].Pos != 2000 {
		t.Errorf("BND AltPositions = %+v", ap)
	}
}

func TestLazyDoesNotParseBadInfo(t *testing.T) {
	// A line with malformed INFO should still yield CHROM/POS/REF/ALT because
	// INFO is only parsed on demand.
	const line = "chr1\t100\t.\tA\tG\t.\t.\tmalformed=info=field=stuff\tGT\t0/0"
	rec, err := newRecord(line, nil)
	if err != nil {
		t.Fatalf("newRecord: %v", err)
	}
	if rec.Chrom != "chr1" || rec.Pos != 100 || rec.Ref != "A" {
		t.Errorf("eager fields wrong: %s:%d %s", rec.Chrom, rec.Pos, rec.Ref)
	}
}

func TestGzipFileRead(t *testing.T) {
	// NewVcfFile transparently reads a BGZF-compressed file via gzip.
	r, err := NewVcfFile("testdata/sample.vcf.gz")
	if err != nil {
		t.Fatalf("NewVcfFile(.gz): %v", err)
	}
	defer r.Close()
	recs := readAll(t, r)
	if len(recs) != 5 {
		t.Fatalf("gz: got %d records, want 5", len(recs))
	}
}

func TestIndexedQuery(t *testing.T) {
	r, err := NewIndexedVcfReader("testdata/sample.vcf.gz")
	if err != nil {
		t.Fatalf("NewIndexedVcfReader: %v", err)
	}
	defer r.Close()

	h, err := r.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if len(h.Samples()) != 2 {
		t.Errorf("indexed header samples = %v", h.Samples())
	}

	// 0-based half-open [199, 1000) -> chr1:200 and chr1:300.
	seq, err := r.Query("chr1", 199, 1000)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}
	var got []int
	for rec, err := range seq {
		if err != nil {
			t.Fatalf("query iter: %v", err)
		}
		got = append(got, rec.Pos)
	}
	if len(got) != 2 || got[0] != 200 || got[1] != 300 {
		t.Errorf("query positions = %v, want [200 300]", got)
	}
}

func TestMain(m *testing.M) {
	os.Exit(m.Run())
}
