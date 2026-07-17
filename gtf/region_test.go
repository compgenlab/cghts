package gtf

import (
	"testing"

	"github.com/compgenlab/cghts/bed"
)

// TestFindGenicRegionForPos covers every genic-region code on the synthetic GTF
// (see testGTF). G1 (+ strand coding) layout, 0-based:
//
//	exon1 [100,200)  5'UTR (before CDS)
//	intron [200,300) 5'UTR intron (< cdsStart 350)
//	exon2 [300,400)  5'UTR [300,350) + CDS [350,400)
//	intron [400,500) coding intron (between cds bounds)
//	exon3 [500,600)  CDS [500,550) + 3'UTR [550,600)
//	intron [600,700) 3'UTR intron (> cdsEnd 550)
//	exon4 [700,800)  3'UTR
func TestFindGenicRegionForPos(t *testing.T) {
	s := loadTestGTF(t)

	cases := []struct {
		name   string
		ref    string
		pos    int
		strand bed.Strand
		want   GenicRegion
	}{
		{"5utr exon", "chr1", 150, bed.StrandNone, UTR5},
		{"5utr exon2", "chr1", 320, bed.StrandNone, UTR5},
		{"coding exon2", "chr1", 370, bed.StrandNone, Coding},
		{"coding exon3", "chr1", 520, bed.StrandNone, Coding},
		{"3utr exon3", "chr1", 570, bed.StrandNone, UTR3},
		{"3utr exon4", "chr1", 750, bed.StrandNone, UTR3},
		{"5utr intron", "chr1", 250, bed.StrandNone, UTR5Intron},
		{"coding intron", "chr1", 450, bed.StrandNone, CodingIntron},
		{"3utr intron", "chr1", 650, bed.StrandNone, UTR3Intron},

		// G2 non-coding.
		{"nc exon", "chr1", 1050, bed.StrandNone, NCExon},
		{"nc intron", "chr1", 1200, bed.StrandNone, NCIntron},

		// G3 minus strand: UTR sides flip.
		{"minus 3utr (before cds)", "chr1", 2020, bed.StrandNone, UTR3},
		{"minus coding", "chr1", 2070, bed.StrandNone, Coding},
		{"minus 5utr (after cds)", "chr1", 2370, bed.StrandNone, UTR5},

		// Intergenic & mitochondrial.
		{"intergenic", "chr1", 5000, bed.StrandNone, Intergenic},
		{"mitochondrial chrM", "chrM", 10, bed.StrandNone, Mitochondrial},
		{"mitochondrial M", "M", 10, bed.StrandNone, Mitochondrial},

		// Anti-sense: query G1 (a + gene) on the minus strand.
		{"antisense coding", "chr1", 370, bed.StrandMinus, CodingAnti},
		{"antisense 5utr", "chr1", 150, bed.StrandMinus, UTR5Anti},
		{"antisense coding intron", "chr1", 450, bed.StrandMinus, CodingIntronAnti},
		// Same strand as the gene → sense.
		{"sense coding (plus query)", "chr1", 370, bed.StrandPlus, Coding},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got := s.FindGenicRegionForPos(tc.ref, tc.pos, tc.strand, "")
			if got != tc.want {
				t.Errorf("region@%s:%d strand=%v = %q, want %q",
					tc.ref, tc.pos, tc.strand, got.Code, tc.want.Code)
			}
		})
	}
}

func TestFindGenicRegionForRegionJunction(t *testing.T) {
	s := loadTestGTF(t)
	// A region whose start is in a coding exon and end is in a coding intron is
	// crossing a splice junction.
	got := s.FindGenicRegionForRegion("chr1", 370, 460, bed.StrandNone)
	if got != Junction {
		t.Errorf("region [370,460) = %q, want junction", got.Code)
	}
	// A region wholly inside one coding exon agrees on both ends.
	if got := s.FindGenicRegionForRegion("chr1", 360, 380, bed.StrandNone); got != Coding {
		t.Errorf("region [360,380) = %q, want coding_exon", got.Code)
	}
}

func TestJunctionize(t *testing.T) {
	if got := Junctionize(Coding); got != Junction {
		t.Errorf("Junctionize(coding_exon) = %q, want junction", got.Code)
	}
	if got := Junctionize(NCExon); got != NCJunction {
		t.Errorf("Junctionize(nc_exon) = %q, want nc_junction", got.Code)
	}
	if got := Junctionize(CodingAnti); got != JunctionAnti {
		t.Errorf("Junctionize(anti_coding_exon) = %q, want anti_junction", got.Code)
	}
	// Non-gene regions pass through.
	if got := Junctionize(Intergenic); got != Intergenic {
		t.Errorf("Junctionize(intergenic) = %q, want intergenic", got.Code)
	}
}

func TestGenicRegionsOrder(t *testing.T) {
	regions := GenicRegions()
	if len(regions) != 22 {
		t.Fatalf("GenicRegions len = %d, want 22", len(regions))
	}
	if regions[0] != Junction || regions[11] != Mitochondrial {
		t.Errorf("priority order wrong: [0]=%q [11]=%q", regions[0].Code, regions[11].Code)
	}
	if Coding.ord >= NCExon.ord {
		t.Errorf("coding_exon should outrank nc_exon (%d vs %d)", Coding.ord, NCExon.ord)
	}
}

func TestPerGeneRegion(t *testing.T) {
	s := loadTestGTF(t)
	// Restricting to a non-overlapping gene yields intergenic.
	if got := s.FindGenicRegionForPos("chr1", 370, bed.StrandNone, "G2"); got != Intergenic {
		t.Errorf("region@370 restricted to G2 = %q, want intergenic", got.Code)
	}
	if got := s.FindGenicRegionForPos("chr1", 370, bed.StrandNone, "G1"); got != Coding {
		t.Errorf("region@370 restricted to G1 = %q, want coding_exon", got.Code)
	}
}
