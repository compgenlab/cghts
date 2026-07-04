package annotate

import "testing"

func TestBigBedAnnotator(t *testing.T) {
	// BED col 4 (name) = geneA over [100,200); col 5 (score) = 42.
	bb := writeTestBigBed(t, "chr1", []tbed{
		{100, 200, "geneA\t42"},
		{300, 400, "geneB\t7"},
	})
	h, recs := bedRecs(t,
		"chr1\t150\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",
		"chr1\t250\t.\tA\tG\t.\tPASS\t.\tGT\t0/1")

	// name column (4)
	a, err := NewBigBedAnnotator(BigBedOptions{Name: "GENE", Filename: bb, Col: 4})
	if err != nil {
		t.Fatalf("NewBigBedAnnotator: %v", err)
	}
	defer a.Close()
	if err := a.SetupHeader(h); err != nil {
		t.Fatal(err)
	}
	for _, rec := range recs {
		if err := a.Annotate(rec); err != nil {
			t.Fatalf("Annotate: %v", err)
		}
	}
	checkInfo(t, recs[0], "GENE", "geneA")
	if _, ok := info(t, recs[1], "GENE"); ok {
		t.Errorf("record at 250 should have no GENE (gap)")
	}
}

func TestBigBedFlag(t *testing.T) {
	bb := writeTestBigBed(t, "chr1", []tbed{{100, 200, "x"}})
	h, recs := bedRecs(t, "chr1\t150\t.\tA\tG\t.\tPASS\t.\tGT\t0/1")
	a, err := NewBigBedAnnotator(BigBedOptions{Name: "INDB", Filename: bb, Col: 0}) // presence flag
	if err != nil {
		t.Fatal(err)
	}
	defer a.Close()
	if err := a.SetupHeader(h); err != nil {
		t.Fatal(err)
	}
	if err := a.Annotate(recs[0]); err != nil {
		t.Fatal(err)
	}
	if !recs[0].Info().Contains("INDB") {
		t.Errorf("expected INDB flag set")
	}
}
