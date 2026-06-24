package annotate

import (
	"bytes"
	"io"
	"strings"
	"testing"

	"github.com/compgenlab/hts/support/stats"
	"github.com/compgenlab/hts/vcf"
)

const hdr = "##fileformat=VCFv4.2\n" +
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
	"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n" +
	"##FORMAT=<ID=SAC,Number=.,Type=Integer,Description=\"Strand allele counts\">\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n"

// readRecs parses an in-memory VCF (header text + data lines) and returns the
// header and records.
func readRecs(t *testing.T, dataLines ...string) (*vcf.VcfHeader, []*vcf.VcfRecord) {
	t.Helper()
	r, err := vcf.NewVcfReader(strings.NewReader(hdr + strings.Join(dataLines, "\n") + "\n"))
	if err != nil {
		t.Fatal(err)
	}
	h, err := r.Header()
	if err != nil {
		t.Fatal(err)
	}
	var recs []*vcf.VcfRecord
	for {
		rec, err := r.NextRecord()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatal(err)
		}
		recs = append(recs, rec)
	}
	return h, recs
}

func infoVal(t *testing.T, rec *vcf.VcfRecord, key string) string {
	t.Helper()
	v, ok := rec.Info().Get(key)
	if !ok {
		t.Fatalf("INFO %s absent", key)
	}
	return v.String()
}

func sampleVal(t *testing.T, rec *vcf.VcfRecord, idx int, key string) string {
	t.Helper()
	s, err := rec.Sample(idx)
	if err != nil {
		t.Fatal(err)
	}
	v, ok := s.Get(key)
	if !ok {
		t.Fatalf("sample %d FORMAT %s absent", idx, key)
	}
	return v.String()
}

func setupAnnotate(t *testing.T, h *vcf.VcfHeader, recs []*vcf.VcfRecord, a Annotator) {
	t.Helper()
	if err := a.SetupHeader(h); err != nil {
		t.Fatalf("SetupHeader: %v", err)
	}
	for _, rec := range recs {
		if err := a.Annotate(rec); err != nil {
			t.Fatalf("Annotate: %v", err)
		}
	}
}

func TestIndel(t *testing.T) {
	h, recs := readRecs(t,
		"chr1\t300\t.\tG\tGA\t.\tPASS\t.\tGT\t0/0\t0/1",
		"chr2\t500\t.\tCAT\tC\t.\tPASS\t.\tGT\t0/0\t0/1")
	setupAnnotate(t, h, recs, NewIndel())

	ins := recs[0]
	if _, ok := ins.Info().Get("CG_INSERT"); !ok {
		t.Error("insertion missing CG_INSERT")
	}
	if got := infoVal(t, ins, "CG_INSLEN"); got != "1" {
		t.Errorf("CG_INSLEN = %q, want 1", got)
	}
	if got := infoVal(t, ins, "CG_INDELLEN"); got != "1" {
		t.Errorf("CG_INDELLEN = %q, want 1", got)
	}
	del := recs[1]
	if _, ok := del.Info().Get("CG_DELETE"); !ok {
		t.Error("deletion missing CG_DELETE")
	}
	if got := infoVal(t, del, "CG_DELLEN"); got != "2" {
		t.Errorf("CG_DELLEN = %q, want 2", got)
	}
	if got := infoVal(t, del, "CG_INDELLEN"); got != "-2" {
		t.Errorf("CG_INDELLEN = %q, want -2", got)
	}
}

func TestTsTv(t *testing.T) {
	h, recs := readRecs(t,
		"chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1",  // transition
		"chr1\t160\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/1",  // transversion
		"chr1\t300\t.\tG\tGA\t.\tPASS\t.\tGT\t0/0\t0/1") // indel -> skipped
	setupAnnotate(t, h, recs, NewTsTv())
	if got := infoVal(t, recs[0], "CG_TSTV"); got != "TS" {
		t.Errorf("A>G = %q, want TS", got)
	}
	if got := infoVal(t, recs[1], "CG_TSTV"); got != "TV" {
		t.Errorf("A>C = %q, want TV", got)
	}
	if _, ok := recs[2].Info().Get("CG_TSTV"); ok {
		t.Error("indel should not get CG_TSTV")
	}
}

func TestAutoID(t *testing.T) {
	h, recs := readRecs(t, "chr1\t100\trs9\tA\tG,T\t.\tPASS\t.\tGT\t0/0\t0/1")
	setupAnnotate(t, h, recs, NewAutoID())
	if got := recs[0].ID(); got != "chr1_100_A_G;chr1_100_A_T" {
		t.Errorf("auto-id = %q", got)
	}
}

func TestConstantTag(t *testing.T) {
	h, recs := readRecs(t, "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1")
	setupAnnotate(t, h, recs, NewConstantFlag("MYFLAG"))
	setupAnnotate(t, h, recs, NewConstantTag("SRC", "panelA"))
	if v, ok := recs[0].Info().Get("MYFLAG"); !ok || !v.IsEmpty() {
		t.Errorf("MYFLAG flag = %+v ok=%v", v, ok)
	}
	if got := infoVal(t, recs[0], "SRC"); got != "panelA" {
		t.Errorf("SRC = %q", got)
	}
}

func TestDosage(t *testing.T) {
	h, recs := readRecs(t, "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1")
	setupAnnotate(t, h, recs, NewDosage())
	if got := sampleVal(t, recs[0], 0, "CG_DS"); got != "0" {
		t.Errorf("NORMAL 0/0 dosage = %q, want 0", got)
	}
	if got := sampleVal(t, recs[0], 1, "CG_DS"); got != "2" {
		t.Errorf("TUMOR 1/1 dosage = %q, want 2", got)
	}
}

func TestVAFAndMinorStrand(t *testing.T) {
	h, recs := readRecs(t,
		"chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT:AD:SAC\t0/0:14,1:15,13,1,1\t0/1:15,15:8,7,8,7")
	setupAnnotate(t, h, recs, NewVAF())
	setupAnnotate(t, h, recs, NewMinorStrand())

	// NORMAL: total=30, alt pair=2 -> 0.067 ; minor: plus==minus -> 1/2 = 0.500
	if got := sampleVal(t, recs[0], 0, "CG_VAF"); got != "0.067" {
		t.Errorf("NORMAL VAF = %q, want 0.067", got)
	}
	if got := sampleVal(t, recs[0], 0, "CG_SBPCT"); got != "0.500" {
		t.Errorf("NORMAL SBPCT = %q, want 0.500", got)
	}
	// TUMOR: total=30, alt pair=15 -> 0.500 ; minor: 7/15 = 0.467
	if got := sampleVal(t, recs[0], 1, "CG_VAF"); got != "0.500" {
		t.Errorf("TUMOR VAF = %q, want 0.500", got)
	}
	if got := sampleVal(t, recs[0], 1, "CG_SBPCT"); got != "0.467" {
		t.Errorf("TUMOR SBPCT = %q, want 0.467", got)
	}
}

func TestFisherStrandBias(t *testing.T) {
	h, recs := readRecs(t,
		"chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT:AD:SAC\t0/0:14,1:15,13,1,1\t0/1:15,15:8,7,8,7")
	setupAnnotate(t, h, recs, NewFisherSB())
	// TUMOR alt pair plus=8,minus=7 -> half=7 -> phred(fisher(7,7,8,7)).
	f := stats.NewFisherExact()
	want := round(stats.Phred(f.TwoTailedPvalue(7, 7, 8, 7)), 3)
	if got := sampleVal(t, recs[0], 1, "CG_FSB"); got != want {
		t.Errorf("TUMOR FSB = %q, want %q", got, want)
	}
}

func TestCopyLogRatio(t *testing.T) {
	h, recs := readRecs(t,
		"chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/0:14,1\t0/1:15,15")
	// somatic=TUMOR (sum 30), germline=NORMAL (sum 15) -> log2(30)-log2(15)=1.
	setupAnnotate(t, h, recs, NewCopyLogRatio("TUMOR", "NORMAL", -1, -1))
	if got := infoVal(t, recs[0], "CG_CNLR"); got != "1.000000" {
		t.Errorf("CNLR = %q, want 1.000000", got)
	}
}

func TestCopyLogRatioMissingAD(t *testing.T) {
	// Header with GT only (no AD def) -> SetupHeader must fail.
	noAD := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n"
	r, _ := vcf.NewVcfReader(strings.NewReader(noAD))
	h, err := r.Header()
	if err != nil {
		t.Fatal(err)
	}
	if err := NewCopyLogRatio("TUMOR", "NORMAL", -1, -1).SetupHeader(h); err == nil {
		t.Error("expected error when AD FORMAT is missing")
	}
}

func TestVariantDistance(t *testing.T) {
	h, recs := readRecs(t,
		"chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1",
		"chr1\t140\t.\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/1",
		"chr1\t150\t.\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/1",
		"chr2\t500\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1")

	vd := NewVariantDistance()
	if err := vd.SetupHeader(h); err != nil {
		t.Fatal(err)
	}
	i := 0
	src := vd.Wrap(func() (*vcf.VcfRecord, error) {
		if i >= len(recs) {
			return nil, io.EOF
		}
		rec := recs[i]
		i++
		return rec, nil
	})
	var got []string
	for {
		rec, err := src()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatal(err)
		}
		v, _ := rec.Info().Get("CG_VARDIST")
		got = append(got, rec.Chrom+":"+v.String())
	}
	// 100 -> nearest is 140 (40); 140 -> min(40,10)=10; 150 -> 10; chr2/500 -> -1 (alone).
	want := []string{"chr1:40", "chr1:10", "chr1:10", "chr2:-1"}
	if strings.Join(got, ",") != strings.Join(want, ",") {
		t.Errorf("vardist = %v, want %v", got, want)
	}
}

func TestPipelineSerialization(t *testing.T) {
	// Run a multi-annotator pipeline and check the reconstructed VCF output.
	r, _ := vcf.NewVcfReader(strings.NewReader(hdr +
		"chr1\t100\t.\tA\tG\t50.0\tPASS\t.\tGT:AD\t0/0:14,1\t1/1:0,30\n"))
	h, _ := r.Header()

	p := NewPipeline()
	p.Add(NewAutoID())
	p.Add(NewTsTv())
	p.Add(NewDosage())
	if err := p.SetupHeaders(h); err != nil {
		t.Fatal(err)
	}
	src := p.Build(r.NextRecord)

	var buf bytes.Buffer
	w := vcf.NewVcfWriter(&buf)
	for {
		rec, err := src()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatal(err)
		}
		if err := w.WriteRecord(rec); err != nil {
			t.Fatal(err)
		}
	}
	w.Close()

	// ID set, QUAL reformatted (50.0->50), INFO has CG_TSTV, FORMAT gains CG_DS.
	want := "chr1\t100\tchr1_100_A_G\tA\tG\t50\tPASS\tCG_TSTV=TS\tGT:AD:CG_DS\t0/0:14,1:0\t1/1:0,30:2\n"
	if buf.String() != want {
		t.Errorf("pipeline output mismatch.\n got: %q\nwant: %q", buf.String(), want)
	}
}

func TestTupleBuilderAnnotation(t *testing.T) {
	// Locus annotators work on a bare tuple; sample-based ones no-op.
	rec := vcf.NewRecord("chr1", 100, "A", []string{"G"})
	for _, a := range []Annotator{NewIndel(), NewTsTv(), NewAutoID(), NewDosage()} {
		if err := a.Annotate(rec); err != nil {
			t.Fatalf("annotate tuple: %v", err)
		}
	}
	if rec.ID() != "chr1_100_A_G" {
		t.Errorf("tuple auto-id = %q", rec.ID())
	}
	if got := infoVal(t, rec, "CG_TSTV"); got != "TS" {
		t.Errorf("tuple tstv = %q", got)
	}
	if rec.NumSamples() != 0 {
		t.Errorf("tuple should have 0 samples, got %d", rec.NumSamples())
	}
}
