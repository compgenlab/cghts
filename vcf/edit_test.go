package vcf

import (
	"strings"
	"testing"
)

func TestRenameSample(t *testing.T) {
	r := openTestFile(t)
	defer r.Close()
	h, err := r.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if err := h.RenameSample("NORMAL", "GERMLINE"); err != nil {
		t.Fatalf("RenameSample(name): %v", err)
	}
	if err := h.RenameSample("2", "CANCER"); err != nil { // by 1-based number
		t.Fatalf("RenameSample(number): %v", err)
	}
	if got := h.Samples(); got[0] != "GERMLINE" || got[1] != "CANCER" {
		t.Errorf("Samples = %v, want [GERMLINE CANCER]", got)
	}
	if h.SampleIndex("CANCER") != 1 || h.SampleIndex("NORMAL") != -1 {
		t.Errorf("index not rebuilt: CANCER=%d NORMAL=%d", h.SampleIndex("CANCER"), h.SampleIndex("NORMAL"))
	}
	if err := h.RenameSample("nope", "X"); err == nil {
		t.Errorf("RenameSample(missing) should error")
	}
	if err := h.RenameSample("GERMLINE", "CANCER"); err == nil {
		t.Errorf("RenameSample(collision) should error")
	}
}

func TestRemoveDefs(t *testing.T) {
	r := openTestFile(t)
	defer r.Close()
	h, err := r.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	h.RemoveInfo("DP")
	if _, ok := h.InfoDef("DP"); ok {
		t.Errorf("RemoveInfo(DP) left the def")
	}
	for _, id := range h.InfoIDs() {
		if id == "DP" {
			t.Errorf("InfoIDs still lists DP: %v", h.InfoIDs())
		}
	}
	h.RemoveFormat("AD")
	if _, ok := h.FormatDef("AD"); ok {
		t.Errorf("RemoveFormat(AD) left the def")
	}
	h.RemoveFilter("lowqual")
	if ids := h.FilterIDs(); len(ids) != 1 || ids[0] != "PASS" {
		t.Errorf("FilterIDs after remove = %v, want [PASS]", ids)
	}
	h.RemoveInfo("nope") // no-op, must not panic
}

func TestRetainFilters(t *testing.T) {
	r := openTestFile(t)
	defer r.Close()
	recs := readAll(t, r)

	// chr1:200 carries "lowqual"; dropping it leaves a non-PASS, empty FILTER (".").
	dropAll := func(string) bool { return false }
	recs[1].RetainFilters(dropAll)
	if got := recs[1].serialize(); !strings.Contains(got, "\t.\tDP=10") {
		t.Errorf("emptied filter should render \".\", got %q", got)
	}

	// chr1:100 is PASS; RetainFilters is a no-op and it stays PASS.
	recs[0].RetainFilters(dropAll)
	if got := recs[0].serialize(); !strings.Contains(got, "\tPASS\t") {
		t.Errorf("PASS record should stay PASS, got %q", got)
	}
}

func TestSetChrom(t *testing.T) {
	r := openTestFile(t)
	defer r.Close()
	recs := readAll(t, r)
	rec := recs[0] // chr1:100
	rec.SetChrom("chrX")
	if !rec.Dirty() {
		t.Errorf("SetChrom did not mark dirty")
	}
	if got := rec.serialize(); !strings.HasPrefix(got, "chrX\t100\t") {
		t.Errorf("serialize after SetChrom = %q", got)
	}
}

func TestSubsetSamples(t *testing.T) {
	r := openTestFile(t)
	defer r.Close()
	recs := readAll(t, r)

	keepTumor := recs[0] // chr1 100 ... GT:AD 0/0:28,2 0/1:15,15
	keepTumor.SubsetSamples([]int{1})
	got := keepTumor.serialize()
	if !strings.HasSuffix(got, "\tGT:AD\t0/1:15,15") {
		t.Errorf("SubsetSamples([1]) = %q", got)
	}
	if n := keepTumor.NumSamples(); n != 1 {
		t.Errorf("NumSamples after subset = %d, want 1", n)
	}

	// Dropping all samples leaves an 8-column record (through INFO).
	rec2 := recs[1]
	rec2.SubsetSamples(nil)
	got = rec2.serialize()
	if strings.Contains(got, "GT") || strings.Count(got, "\t") != 7 {
		t.Errorf("SubsetSamples(nil) = %q, want 8 columns, no FORMAT", got)
	}
}

func TestGenotypeBases(t *testing.T) {
	r := openTestFile(t)
	defer r.Close()
	recs := readAll(t, r)

	// chr1:100 A>G, NORMAL 0/0 -> A/A, TUMOR 0/1 -> A/G
	if gt, ok := recs[0].GenotypeBases(0); !ok || gt != "A/A" {
		t.Errorf("GenotypeBases(0) = %q,%v want A/A", gt, ok)
	}
	if gt, ok := recs[0].GenotypeBases(1); !ok || gt != "A/G" {
		t.Errorf("GenotypeBases(1) = %q,%v want A/G", gt, ok)
	}
	// chr1:300 G>GA, TUMOR 0/1 -> G/GA
	if gt, ok := recs[2].GenotypeBases(1); !ok || gt != "G/GA" {
		t.Errorf("GenotypeBases indel = %q,%v want G/GA", gt, ok)
	}
}

func TestResolveGenotypeBases(t *testing.T) {
	alts := []string{"G", "T"}
	cases := []struct {
		gt   string
		want string
		ok   bool
	}{
		{"0/0", "A/A", true},
		{"0/1", "A/G", true},
		{"1/2", "G/T", true},
		{"0|1", "A/G", true}, // phased input, "/"-joined output
		{"./.", "", false},
		{"0", "", false},   // not diploid
		{"0/9", "", false}, // unknown allele index
		{"0/.", "", false}, // missing allele
	}
	for _, c := range cases {
		got, ok := resolveGenotypeBases(c.gt, "A", alts)
		if got != c.want || ok != c.ok {
			t.Errorf("resolveGenotypeBases(%q) = %q,%v want %q,%v", c.gt, got, ok, c.want, c.ok)
		}
	}
}

func TestChromConvert(t *testing.T) {
	cases := []struct{ in, ucsc, ensembl string }{
		{"1", "chr1", "1"},
		{"chr1", "chr1", "1"},
		{"MT", "chrM", "MT"},
		{"chrM", "chrM", "MT"},
		{"X", "chrX", "X"},
	}
	for _, c := range cases {
		if got := ToUCSC(c.in); got != c.ucsc {
			t.Errorf("ToUCSC(%q) = %q want %q", c.in, got, c.ucsc)
		}
		if got := ToEnsembl(c.in); got != c.ensembl {
			t.Errorf("ToEnsembl(%q) = %q want %q", c.in, got, c.ensembl)
		}
	}
	if !IsPrimaryHuman("chr22") || !IsPrimaryHuman("X") || IsPrimaryHuman("chr1_KI270706v1_random") {
		t.Errorf("IsPrimaryHuman classification wrong")
	}
}
