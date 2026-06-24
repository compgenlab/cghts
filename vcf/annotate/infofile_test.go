package annotate

import (
	"testing"

	"github.com/compgenlab/hts/vcf"
)

func runInFile(t *testing.T, opts InfoFileOptions, h *vcf.VcfHeader, recs []*vcf.VcfRecord) {
	t.Helper()
	a, err := NewInfoInFile(opts)
	if err != nil {
		t.Fatalf("NewInfoInFile: %v", err)
	}
	if err := a.SetupHeader(h); err != nil {
		t.Fatalf("SetupHeader: %v", err)
	}
	for _, rec := range recs {
		if err := a.Annotate(rec); err != nil {
			t.Fatalf("Annotate: %v", err)
		}
	}
}

func TestInfoInFileFlag(t *testing.T) {
	h, recs := bedRecs(t,
		"chr1\t100\t.\tA\tG\t.\tPASS\tDP=30\tGT\t0/1",
		"chr1\t140\t.\tC\tT\t.\tPASS\tDP=10\tGT\t0/1",
		"chr1\t150\t.\tA\tC\t.\tPASS\tDP=25\tGT\t0/1")
	runInFile(t, InfoFileOptions{Filename: "testdata/dp_set.txt", Tag: "DP", FlagName: "HIT"}, h, recs)
	for i, want := range []bool{true, false, true} {
		_, ok := recs[i].Info().Get("HIT")
		if ok != want {
			t.Errorf("rec%d HIT present=%v, want %v", i, ok, want)
		}
	}
	if d, ok := h.InfoDef("HIT"); !ok || d.Type != "Flag" {
		t.Errorf("HIT def = %+v ok=%v", d, ok)
	}
}

func TestInfoInFileTabcol(t *testing.T) {
	h, recs := bedRecs(t,
		"chr1\t100\t.\tA\tG\t.\tPASS\tDP=30\tGT\t0/1",
		"chr1\t150\t.\tA\tC\t.\tPASS\tDP=25\tGT\t0/1",
		"chr1\t300\t.\tG\tGA\t.\tPASS\tDP=99\tGT\t0/1") // not in file
	runInFile(t, InfoFileOptions{Filename: "testdata/dp_val.txt", Tag: "DP", FlagName: "LABEL", Col: 2}, h, recs)
	if v, _ := recs[0].Info().Get("LABEL"); v.String() != "common" {
		t.Errorf("rec0 LABEL = %q, want common", v.String())
	}
	if v, _ := recs[1].Info().Get("LABEL"); v.String() != "rare" {
		t.Errorf("rec1 LABEL = %q, want rare", v.String())
	}
	if _, ok := recs[2].Info().Get("LABEL"); ok {
		t.Error("rec2 (DP=99 not in file) should have no LABEL")
	}
}

func TestInfoInFileCSV(t *testing.T) {
	// The record value is comma-delimited; csv splits it and matches any piece.
	h, recs := bedRecs(t,
		"chr1\t100\t.\tA\tG\t.\tPASS\tDPS=99,30\tGT\t0/1", // 30 is in the file
		"chr1\t140\t.\tC\tT\t.\tPASS\tDPS=11,12\tGT\t0/1") // neither in file
	runInFile(t, InfoFileOptions{Filename: "testdata/dp_set.txt", Tag: "DPS", FlagName: "HIT", Delimiter: ","}, h, recs)
	if _, ok := recs[0].Info().Get("HIT"); !ok {
		t.Error("rec0 should match (30 present in csv value)")
	}
	if _, ok := recs[1].Info().Get("HIT"); ok {
		t.Error("rec1 should not match")
	}
}

func TestInfoInFileMissingTag(t *testing.T) {
	// No DP INFO on the record -> no annotation, no error.
	h, recs := bedRecs(t, "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1")
	runInFile(t, InfoFileOptions{Filename: "testdata/dp_set.txt", Tag: "DP", FlagName: "HIT"}, h, recs)
	if _, ok := recs[0].Info().Get("HIT"); ok {
		t.Error("record without DP should not be flagged")
	}
}
