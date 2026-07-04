package annotate

import "testing"

func TestBigWigAnnotator(t *testing.T) {
	// value at chr1:150 (0-based [149,150)) = 0.42; nothing at 500.
	bw := writeTestBigWig(t, "chr1", []twig{
		{149, 150, 0.42},
		{199, 200, 1.5},
	})
	h, recs := bedRecs(t,
		"chr1\t150\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",
		"chr1\t500\t.\tA\tG\t.\tPASS\t.\tGT\t0/1")

	a, err := NewBigWigAnnotator(BigWigOptions{Name: "SCORE", Filename: bw})
	if err != nil {
		t.Fatalf("NewBigWigAnnotator: %v", err)
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
	checkInfo(t, recs[0], "SCORE", "0.42")
	if _, ok := info(t, recs[1], "SCORE"); ok {
		t.Errorf("record at 500 should have no SCORE")
	}
}
