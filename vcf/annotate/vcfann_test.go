package annotate

import (
	"testing"

	"github.com/compgenlab/hts/vcf"
)

func runVcfAnn(t *testing.T, opts VcfOptions, h *vcf.VcfHeader, recs []*vcf.VcfRecord) *VcfAnnotation {
	t.Helper()
	a, err := NewVcfAnnotation(opts)
	if err != nil {
		t.Fatalf("NewVcfAnnotation: %v", err)
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

// vcfTargets returns target records matching source.vcf.gz positions plus a
// control on a chromosome absent from the source.
func vcfTargets(t *testing.T) (*vcf.VcfHeader, []*vcf.VcfRecord) {
	return bedRecs(t,
		"chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",   // matches rsA (PASS)
		"chr1\t150\t.\tA\tC\t.\tPASS\t.\tGT\t0/1",   // matches rsB (rej)
		"chr2\t500\t.\tCAT\tC\t.\tPASS\t.\tGT\t0/1", // matches rsC (deletion)
		"chr9\t900\t.\tA\tT\t.\tPASS\t.\tGT\t0/1")   // chrom absent from source
}

func TestVcfAnnInfoField(t *testing.T) {
	h, recs := vcfTargets(t)
	a := runVcfAnn(t, VcfOptions{Name: "KAF", Field: "AF", Filename: "testdata/source.vcf.gz"}, h, recs)
	defer a.Close()
	if v, _ := recs[0].Info().Get("KAF"); v.String() != "0.20" {
		t.Errorf("rec0 KAF = %q, want 0.20", v.String())
	}
	if v, _ := recs[1].Info().Get("KAF"); v.String() != "0.10" {
		t.Errorf("rec1 KAF = %q, want 0.10", v.String())
	}
	if _, ok := recs[2].Info().Get("KAF"); ok {
		t.Error("rec2 (no AF in source) should have no KAF")
	}
	if _, ok := recs[3].Info().Get("KAF"); ok {
		t.Error("rec3 (absent chrom) should have no KAF")
	}
	if d, ok := h.InfoDef("KAF"); !ok || d.Type != "String" {
		t.Errorf("KAF def = %+v ok=%v", d, ok)
	}
}

func TestVcfAnnPassingOnly(t *testing.T) {
	h, recs := vcfTargets(t)
	a := runVcfAnn(t, VcfOptions{Name: "KAF", Field: "AF", Filename: "testdata/source.vcf.gz", Passing: true}, h, recs)
	defer a.Close()
	if v, _ := recs[0].Info().Get("KAF"); v.String() != "0.20" {
		t.Errorf("rec0 KAF = %q, want 0.20", v.String())
	}
	if _, ok := recs[1].Info().Get("KAF"); ok {
		t.Error("rec1 (rsB is rej) should be skipped with Passing")
	}
}

func TestVcfAnnFlag(t *testing.T) {
	h, recs := vcfTargets(t)
	a := runVcfAnn(t, VcfOptions{Name: "KNOWN", Filename: "testdata/source.vcf.gz"}, h, recs)
	defer a.Close()
	for i, want := range []bool{true, true, true, false} {
		_, ok := recs[i].Info().Get("KNOWN")
		if ok != want {
			t.Errorf("rec%d KNOWN present=%v, want %v", i, ok, want)
		}
	}
	if d, ok := h.InfoDef("KNOWN"); !ok || d.Type != "Flag" {
		t.Errorf("KNOWN def = %+v ok=%v", d, ok)
	}
}

func TestVcfAnnIDCopy(t *testing.T) {
	h, recs := vcfTargets(t)
	a := runVcfAnn(t, VcfOptions{Name: "@ID", Filename: "testdata/source.vcf.gz"}, h, recs)
	defer a.Close()
	for i, want := range []string{"rsA", "rsB", "rsC", ""} {
		if recs[i].ID() != want {
			t.Errorf("rec%d ID = %q, want %q", i, recs[i].ID(), want)
		}
	}
	// @ID adds no header def.
	if _, ok := h.InfoDef("@ID"); ok {
		t.Error("@ID should not add a header def")
	}
}

func TestVcfAnnExactMatch(t *testing.T) {
	// A target whose ALT differs from the source must not match under exact.
	h, recs := bedRecs(t, "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/1") // source rsA is A>G
	a := runVcfAnn(t, VcfOptions{Name: "KAF", Field: "AF", Filename: "testdata/source.vcf.gz", Exact: true}, h, recs)
	defer a.Close()
	if _, ok := recs[0].Info().Get("KAF"); ok {
		t.Error("exact match should reject A>T against source A>G")
	}
	// Without exact, position match alone is enough.
	h2, recs2 := bedRecs(t, "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/1")
	a2 := runVcfAnn(t, VcfOptions{Name: "KAF", Field: "AF", Filename: "testdata/source.vcf.gz"}, h2, recs2)
	defer a2.Close()
	if v, ok := recs2[0].Info().Get("KAF"); !ok || v.String() != "0.20" {
		t.Errorf("non-exact KAF = %q ok=%v, want 0.20", v.String(), ok)
	}
}
