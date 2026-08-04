package varstore

import (
	"os"
	"testing"
)

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

// The two backends promise identical answers, and RefEnd is part of an answer:
// a consumer doing overlap arithmetic on returned rows reads 0 as "one base
// wide". The VCF side computed the extent for its own filtering and then
// dropped it when building the row, so every deletion came back one base wide
// from a VCF and its true width from the store converted out of that same VCF.
func TestBackendsAgreeOnRefEnd(t *testing.T) {
	dir := t.TempDir()
	vcfPath := dir + "/in.vcf"
	body := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" +
		"chr1\t1000\t.\tATTTT\tA\t.\t.\t.\tGT\t0/1\n" +
		"chr1\t2000\t.\tC\tG\t.\t.\t.\tGT\t0/1\n"
	if err := os.WriteFile(vcfPath, []byte(body), 0o644); err != nil {
		t.Fatal(err)
	}

	vs, err := OpenVcf(vcfPath)
	if err != nil {
		t.Fatal(err)
	}
	defer vs.Close()
	fromVcf, err := CollectCalls(vs, Query{})
	if err != nil {
		t.Fatal(err)
	}
	if len(fromVcf) != 2 {
		t.Fatalf("got %d calls from the VCF, want 2", len(fromVcf))
	}

	// A 5-base REF starting at 1000 covers 1-based 1000..1004, so the 0-based
	// exclusive end is 1004; a 1-base REF at 2000 ends at 2000.
	for _, want := range []struct {
		pos, refEnd int32
	}{{1000, 1004}, {2000, 2000}} {
		var got int32 = -1
		for _, c := range fromVcf {
			if c.Pos == want.pos {
				got = c.RefEnd
			}
		}
		if got != want.refEnd {
			t.Errorf("VCF call at %d has RefEnd %d, want %d", want.pos, got, want.refEnd)
		}
	}

	// And the store converted from it must say the same.
	base := dir + "/store"
	w, err := NewWriter(base, WriterOpts{Samples: []string{"S1"}, MinDP: 0})
	if err != nil {
		t.Fatal(err)
	}
	for _, c := range fromVcf {
		if err := w.WriteCall(c); err != nil {
			t.Fatal(err)
		}
		if err := w.WriteSite(Site{Chrom: c.Chrom, Pos: c.Pos, Ref: c.Ref, Alt: c.Alt, RefEnd: c.RefEnd}); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
	ps, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer ps.Close()
	fromStore, err := CollectCalls(ps, Query{})
	if err != nil {
		t.Fatal(err)
	}
	if len(fromStore) != len(fromVcf) {
		t.Fatalf("store returned %d calls, VCF %d", len(fromStore), len(fromVcf))
	}
	for i := range fromVcf {
		if fromVcf[i].RefEnd != fromStore[i].RefEnd {
			t.Errorf("row %d: VCF RefEnd %d, store RefEnd %d",
				i, fromVcf[i].RefEnd, fromStore[i].RefEnd)
		}
	}
}
