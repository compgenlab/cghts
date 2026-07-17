package annotate

import (
	"testing"

	"github.com/compgenlab/cghts/vcf"
)

// infoKeyEq asserts an INFO key is identical (presence + value) on two records.
func infoKeyEq(t *testing.T, i int, key string, base, grp *vcf.VcfRecord) {
	t.Helper()
	vb, okb := base.Info().Get(key)
	vg, okg := grp.Info().Get(key)
	if okb != okg {
		t.Errorf("rec%d %s: presence baseline=%v group=%v", i, key, okb, okg)
		return
	}
	if okb && vb.String() != vg.String() {
		t.Errorf("rec%d %s: baseline=%q group=%q", i, key, vb.String(), vg.String())
	}
}

// TestVcfGroupMatchesBaseline: NewVcfAnnotationGroup must produce byte-identical
// output to running the equivalent single VcfAnnotation per field — across a value
// field, an exact-match field, a flag, and an @ID copy.
func TestVcfGroupMatchesBaseline(t *testing.T) {
	fields := []VcfFieldOptions{
		{Name: "KAF", Field: "AF"},               // value, non-exact
		{Name: "KG", Field: "GENE", Exact: true}, // value, exact
		{Name: "KNOWN"},                          // flag
		{Name: "@ID"},                            // ID copy (forces exact)
	}

	// Baseline: one single annotator per field, applied in order.
	hB, recsB := vcfTargets(t)
	for _, f := range fields {
		a, err := NewVcfAnnotation(VcfOptions{Name: f.Name, Field: f.Field, Exact: f.Exact, Filename: "testdata/source.vcf.gz"})
		if err != nil {
			t.Fatal(err)
		}
		if err := a.SetupHeader(hB); err != nil {
			t.Fatal(err)
		}
		for _, rec := range recsB {
			if err := a.Annotate(rec); err != nil {
				t.Fatal(err)
			}
		}
		a.Close()
	}

	// Grouped: one reader, one query per record.
	hG, recsG := vcfTargets(t)
	g, err := NewVcfAnnotationGroup(VcfGroupOptions{Filename: "testdata/source.vcf.gz", Fields: fields})
	if err != nil {
		t.Fatal(err)
	}
	defer g.Close()
	if err := g.SetupHeader(hG); err != nil {
		t.Fatal(err)
	}
	for _, rec := range recsG {
		if err := g.Annotate(rec); err != nil {
			t.Fatal(err)
		}
	}

	for i := range recsB {
		if recsB[i].ID() != recsG[i].ID() {
			t.Errorf("rec%d ID: baseline=%q group=%q", i, recsB[i].ID(), recsG[i].ID())
		}
		for _, key := range []string{"KAF", "KG", "KNOWN"} {
			infoKeyEq(t, i, key, recsB[i], recsG[i])
		}
	}
}

// TestTabixGroupMatchesBaseline: NewTabixAnnotationGroup must match two single
// TabixAnnotators sharing the source's ref/alt match columns — one numeric (Max),
// one string — over records incl. a position with multiple alt rows.
func TestTabixGroupMatchesBaseline(t *testing.T) {
	lines := []string{
		"chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1", // chr1:100 has A>G and A>T rows
		"chr1\t150\t.\tA\tC\t.\tPASS\t.\tGT\t0/1",
		"chr2\t500\t.\tC\tA\t.\tPASS\t.\tGT\t0/1",
		"chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/1",
	}

	// Baseline: two single annotators applied to the same records.
	hB, recsB := bedRecs(t, lines...)
	aSC := runTabix(t, TabixOptions{Name: "SC", Filename: "testdata/scores.tab.gz", Col: 5, AltCol: 4, RefCol: 3, IsNumber: true, Max: true}, hB, recsB)
	defer aSC.Close()
	aLB := runTabix(t, TabixOptions{Name: "LB", Filename: "testdata/scores.tab.gz", Col: 6, AltCol: 4, RefCol: 3}, hB, recsB)
	defer aLB.Close()

	// Grouped: one reader, one query per record; shared AltCol/RefCol.
	hG, recsG := bedRecs(t, lines...)
	g, err := NewTabixAnnotationGroup(TabixGroupOptions{
		Filename: "testdata/scores.tab.gz", AltCol: 4, RefCol: 3,
		Fields: []TabixFieldOptions{
			{Name: "SC", Col: 5, IsNumber: true, Max: true},
			{Name: "LB", Col: 6},
		},
	})
	if err != nil {
		t.Fatal(err)
	}
	defer g.Close()
	if err := g.SetupHeader(hG); err != nil {
		t.Fatal(err)
	}
	for _, rec := range recsG {
		if err := g.Annotate(rec); err != nil {
			t.Fatal(err)
		}
	}

	for i := range recsB {
		for _, key := range []string{"SC", "LB"} {
			infoKeyEq(t, i, key, recsB[i], recsG[i])
		}
	}
}
