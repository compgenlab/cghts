package annotate

import (
	"strings"
	"testing"

	"github.com/compgenlab/hts/vcf"
)

const bedHdr = "##fileformat=VCFv4.2\n" +
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"

func bedRecs(t *testing.T, lines ...string) (*vcf.VcfHeader, []*vcf.VcfRecord) {
	t.Helper()
	return readRecsHdr(t, bedHdr, lines...)
}

// readRecsHdr is like readRecs but with a caller-supplied header.
func readRecsHdr(t *testing.T, header string, lines ...string) (*vcf.VcfHeader, []*vcf.VcfRecord) {
	t.Helper()
	var text string
	for _, l := range lines {
		text += l + "\n"
	}
	r, err := vcf.NewVcfReader(strings.NewReader(header + text))
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
		if err != nil {
			break
		}
		recs = append(recs, rec)
	}
	return h, recs
}

func runTabix(t *testing.T, opts TabixOptions, h *vcf.VcfHeader, recs []*vcf.VcfRecord) *TabixAnnotator {
	t.Helper()
	a, err := NewTabixAnnotator(opts)
	if err != nil {
		t.Fatalf("NewTabixAnnotator: %v", err)
	}
	if err := a.SetupHeader(h); err != nil {
		t.Fatalf("SetupHeader: %v", err)
	}
	for _, rec := range recs {
		if err := a.Annotate(rec); err != nil {
			t.Fatalf("Annotate: %v", err)
		}
	}
	return a
}

func TestTabixBED(t *testing.T) {
	h, recs := bedRecs(t,
		"chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",
		"chr1\t150\t.\tA\tC\t.\tPASS\t.\tGT\t0/1",
		"chr2\t500\t.\tC\tA\t.\tPASS\t.\tGT\t0/1",
		"chr1\t200\t.\tA\tG\t.\tPASS\t.\tGT\t0/1") // outside any region
	a := runTabix(t, TabixOptions{Name: "REGION", Filename: "testdata/regions.bed.gz", Col: 4}, h, recs)
	defer a.Close()
	want := []string{"geneA", "enhB", "geneC", ""}
	for i, w := range want {
		v, ok := recs[i].Info().Get("REGION")
		if w == "" {
			if ok {
				t.Errorf("rec %d should have no REGION, got %q", i, v.String())
			}
			continue
		}
		if !ok || v.String() != w {
			t.Errorf("rec %d REGION = %q (ok=%v), want %q", i, v.String(), ok, w)
		}
	}
	// Header def is present and tabix-style.
	if d, ok := h.InfoDef("REGION"); !ok || d.Type != "String" || d.Number != "." {
		t.Errorf("REGION def = %+v ok=%v", d, ok)
	}
}

func TestTabixBedFlag(t *testing.T) {
	h, recs := bedRecs(t,
		"chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",
		"chr1\t200\t.\tA\tG\t.\tPASS\t.\tGT\t0/1")
	a := runTabix(t, TabixOptions{Name: "INREGION", Filename: "testdata/regions.bed.gz", Col: 0}, h, recs)
	defer a.Close()
	if v, ok := recs[0].Info().Get("INREGION"); !ok || !v.IsEmpty() {
		t.Errorf("rec0 INREGION flag = %+v ok=%v", v, ok)
	}
	if _, ok := recs[1].Info().Get("INREGION"); ok {
		t.Errorf("rec1 should not have INREGION flag")
	}
	if d, ok := h.InfoDef("INREGION"); !ok || d.Type != "Flag" {
		t.Errorf("INREGION def = %+v ok=%v", d, ok)
	}
}

func TestTabixFormatBed(t *testing.T) {
	h, recs := bedRecs(t, "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1")
	a := runTabix(t, TabixOptions{Name: "REGION", Filename: "testdata/regions.bed.gz", Col: 4, Sample: "S1"}, h, recs)
	defer a.Close()
	s, _ := recs[0].Sample(0)
	if v, ok := s.Get("REGION"); !ok || v.String() != "geneA" {
		t.Errorf("FORMAT REGION = %q ok=%v, want geneA", v.String(), ok)
	}
	if d, ok := h.FormatDef("REGION"); !ok || d.Number != "." {
		t.Errorf("REGION FORMAT def = %+v ok=%v", d, ok)
	}
}

func TestTabixColMultiAndAggregate(t *testing.T) {
	// chr1:100 has two tab rows (A>G score 0.9, A>T score 0.5).
	mk := func() (*vcf.VcfHeader, []*vcf.VcfRecord) {
		return bedRecs(t, "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1")
	}
	type tc struct {
		opts TabixOptions
		want string
	}
	for _, c := range []tc{
		{TabixOptions{Name: "S", Filename: "testdata/scores.tab.gz", Col: 5}, "0.9,0.5"},
		{TabixOptions{Name: "S", Filename: "testdata/scores.tab.gz", Col: 5, First: true}, "0.9"},
		{TabixOptions{Name: "S", Filename: "testdata/scores.tab.gz", Col: 5, Collapse: true}, "0.9,0.5"},
		{TabixOptions{Name: "S", Filename: "testdata/scores.tab.gz", Col: 5, Max: true, IsNumber: true}, "0.9"},
		// Restrict to the matching alt allele (G) -> only the 0.9 row.
		{TabixOptions{Name: "S", Filename: "testdata/scores.tab.gz", Col: 5, AltCol: 4}, "0.9"},
		// Match ref+alt exactly.
		{TabixOptions{Name: "S", Filename: "testdata/scores.tab.gz", Col: 6, AltCol: 4, RefCol: 3}, "hotspot"},
	} {
		h, recs := mk()
		a := runTabix(t, c.opts, h, recs)
		v, ok := recs[0].Info().Get("S")
		if !ok || v.String() != c.want {
			t.Errorf("opts %+v -> %q (ok=%v), want %q", c.opts, v.String(), ok, c.want)
		}
		a.Close()
	}
}

func TestTabixColumnByName(t *testing.T) {
	// scores_hdr.tab.gz has a header (#chrom pos ref alt score label) indexed
	// with Skip(1), so columns can be referenced by name.
	h, recs := bedRecs(t, "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1")
	a := runTabix(t, TabixOptions{
		Name:     "LBL",
		Filename: "testdata/scores_hdr.tab.gz",
		ColName:  "label",
		AltName:  "alt",
		RefName:  "ref",
	}, h, recs)
	defer a.Close()
	// chr1:100 A>G -> the "hot" row (alt=G), label column.
	if v, ok := recs[0].Info().Get("LBL"); !ok || v.String() != "hot" {
		t.Errorf("LBL by name = %q ok=%v, want hot", v.String(), ok)
	}
}

func TestTabixColumnByNameMissing(t *testing.T) {
	_, err := NewTabixAnnotator(TabixOptions{
		Name: "X", Filename: "testdata/scores_hdr.tab.gz", ColName: "nope",
	})
	if err == nil {
		t.Error("expected error for unknown column name")
	}
}

func TestTabixColumnByNameNoHeaderFile(t *testing.T) {
	// scores.tab.gz has no header line (no skipped line); requesting a column by
	// name must fail.
	_, err := NewTabixAnnotator(TabixOptions{
		Name: "X", Filename: "testdata/scores.tab.gz", ColName: "score",
	})
	if err == nil {
		t.Error("expected error: cannot resolve a column name on a headerless file")
	}
}

func TestTabixExtend(t *testing.T) {
	// chr1:115 is just outside the BED region [90,110); --extend reaches it.
	h, recs := bedRecs(t, "chr1\t115\t.\tA\tG\t.\tPASS\t.\tGT\t0/1")
	a := runTabix(t, TabixOptions{Name: "REGION", Filename: "testdata/regions.bed.gz", Col: 4, Extend: 10}, h, recs)
	defer a.Close()
	if v, ok := recs[0].Info().Get("REGION"); !ok || v.String() != "geneA" {
		t.Errorf("extended REGION = %q ok=%v, want geneA", v.String(), ok)
	}

	h2, recs2 := bedRecs(t, "chr1\t115\t.\tA\tG\t.\tPASS\t.\tGT\t0/1")
	a2 := runTabix(t, TabixOptions{Name: "REGION", Filename: "testdata/regions.bed.gz", Col: 4}, h2, recs2)
	defer a2.Close()
	if _, ok := recs2[0].Info().Get("REGION"); ok {
		t.Errorf("without extend, chr1:115 should be outside the region")
	}
}

func TestTabixMissingSample(t *testing.T) {
	h, _ := bedRecs(t, "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1")
	a, err := NewTabixAnnotator(TabixOptions{Name: "X", Filename: "testdata/regions.bed.gz", Col: 4, Sample: "NOPE"})
	if err != nil {
		t.Fatal(err)
	}
	defer a.Close()
	if err := a.SetupHeader(h); err == nil {
		t.Error("expected error for missing sample")
	}
}

func TestBaseCoordResolution(t *testing.T) {
	_, recs := bedRecs(t, "chr1\t100\t.\tCAT\tC\t.\tPASS\t.\tGT\t0/1") // deletion
	var b base
	// Deletion: query position is pos+1, end is pos-1+len(ref) = 100+3 = 103.
	if p, ok := b.Pos(recs[0]); !ok || p != 101 {
		t.Errorf("Pos(deletion) = %d ok=%v, want 101", p, ok)
	}
	if e, ok := b.EndPos(recs[0]); !ok || e != 103 {
		t.Errorf("EndPos(deletion) = %d ok=%v, want 103", e, ok)
	}
	// SNV path.
	_, snv := bedRecs(t, "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1")
	if p, _ := b.Pos(snv[0]); p != 100 {
		t.Errorf("Pos(snv) = %d, want 100", p)
	}
	// alt-chrom/alt-pos from INFO.
	_, sv := readRecsHdr(t,
		"##fileformat=VCFv4.2\n##INFO=<ID=CHR2,Number=1,Type=String,Description=\"x\">\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"x\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
		"chr1\t100\t.\tA\t<DEL>\t.\tPASS\tCHR2=chr5;END=900\n")
	var b2 base
	b2.SetAltChrom("CHR2")
	b2.SetAltPos("END")
	if c, ok := b2.Chrom(sv[0]); !ok || c != "chr5" {
		t.Errorf("alt Chrom = %q ok=%v, want chr5", c, ok)
	}
	if p, ok := b2.Pos(sv[0]); !ok || p != 900 {
		t.Errorf("alt Pos = %d ok=%v, want 900", p, ok)
	}
}
