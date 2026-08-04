package varstore

import "testing"

func TestWriterDerivesRefEndFromRef(t *testing.T) {
	base := t.TempDir() + "/d"
	w, _ := NewWriter(base, WriterOpts{Samples: []string{"S1"}, MinDP: 10, RowGroupSize: 1})
	ref := "A"
	for i := 0; i < 199; i++ {
		ref += "C"
	}
	// RefEnd deliberately left zero, as a caller predating the column would.
	w.WriteCall(Call{SampleID: "S1", Chrom: "chr1", Pos: 3000, Ref: ref, Alt: "A", GT: "0/1", DP: 30})
	w.WriteSite(Site{Chrom: "chr1", Pos: 3000, Ref: ref, Alt: "A", AC: 1, AN: 2, NCarriers: 1, NCalled: 1})
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
	s, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	got, err := CollectCalls(s, Query{Spans: []Span{{Chrom: "chr1", Start: 3100, End: 3150}}})
	if err != nil {
		t.Fatal(err)
	}
	if len(got) != 1 {
		t.Fatalf("span inside a 200bp deletion gave %d rows, want 1", len(got))
	}
	if got[0].RefEnd != 2999+200 {
		t.Errorf("RefEnd = %d, want %d", got[0].RefEnd, 2999+200)
	}
}
