package annotate

import (
	"testing"

	"github.com/compgenlab/hts/vcf"
)

func runFlanking(t *testing.T, opts FlankingOptions, h *vcf.VcfHeader, recs []*vcf.VcfRecord) {
	t.Helper()
	a, err := NewFlankingBases(opts)
	if err != nil {
		t.Fatalf("NewFlankingBases: %v", err)
	}
	defer a.Close()
	if err := a.SetupHeader(h); err != nil {
		t.Fatalf("SetupHeader: %v", err)
	}
	for _, rec := range recs {
		if err := a.Annotate(rec); err != nil {
			t.Fatalf("Annotate: %v", err)
		}
	}
}

// ref.fa is chr1 = TACGTCAGT (1-based positions T A C G T C A G T).
func TestFlankingSize1(t *testing.T) {
	h, recs := bedRecs(t,
		"chr1\t2\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",  // A>G, REF=A -> revcomp
		"chr1\t3\t.\tC\tA\t.\tPASS\t.\tGT\t0/1",  // C>A, no revcomp
		"chr1\t1\t.\tT\tA\t.\tPASS\t.\tGT\t0/1",  // boundary -> skipped
		"chr1\t5\t.\tCA\tC\t.\tPASS\t.\tGT\t0/1") // indel -> skipped
	runFlanking(t, FlankingOptions{Filename: "testdata/ref.fa", Size: 1}, h, recs)

	check := func(i int, flank, sub string) {
		v, ok := recs[i].Info().Get("CG_FLANKING")
		if (flank == "") != !ok {
			t.Errorf("rec%d FLANKING present=%v, want %q", i, ok, flank)
		} else if ok && v.String() != flank {
			t.Errorf("rec%d FLANKING=%q, want %q", i, v.String(), flank)
		}
		sv, sok := recs[i].Info().Get("CG_FLANKING_SUB")
		if (sub == "") != !sok {
			t.Errorf("rec%d SUB present=%v, want %q", i, sok, sub)
		} else if sok && sv.String() != sub {
			t.Errorf("rec%d SUB=%q, want %q", i, sv.String(), sub)
		}
	}
	// chr1:2 flanking (1..3)="TAC"; REF=A -> revcomp("TAC")="GTA"; alt G -> C.
	check(0, "TAC", "G[T>C]A")
	// chr1:3 flanking (2..4)="ACG"; no revcomp.
	check(1, "ACG", "A[C>A]G")
	// boundary and indel -> nothing.
	check(2, "", "")
	check(3, "", "")
}

func TestFlankingBgzip(t *testing.T) {
	// Same reference, bgzip-compressed + faidx'd; OpenReference handles it.
	h, recs := bedRecs(t, "chr1\t3\t.\tC\tA\t.\tPASS\t.\tGT\t0/1")
	runFlanking(t, FlankingOptions{Filename: "testdata/ref.fa.gz", Size: 1}, h, recs)
	if v, _ := recs[0].Info().Get("CG_FLANKING_SUB"); v.String() != "A[C>A]G" {
		t.Errorf("bgzip SUB=%q, want A[C>A]G", v.String())
	}
}

func TestFlankingSize2(t *testing.T) {
	h, recs := bedRecs(t, "chr1\t5\t.\tT\tA\t.\tPASS\t.\tGT\t0/1")
	runFlanking(t, FlankingOptions{Filename: "testdata/ref.fa", Size: 2}, h, recs)
	// chr1:5 +/-2 flanking (3..7)="CGTCA"; REF=T -> no revcomp.
	if v, _ := recs[0].Info().Get("CG_FLANKING"); v.String() != "CGTCA" {
		t.Errorf("FLANKING=%q, want CGTCA", v.String())
	}
	if v, _ := recs[0].Info().Get("CG_FLANKING_SUB"); v.String() != "CG[T>A]CA" {
		t.Errorf("SUB=%q, want CG[T>A]CA", v.String())
	}
}
