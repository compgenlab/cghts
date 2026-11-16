package gtf

import (
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/bed"
)

func parseGTF(t *testing.T, lines ...string) *AnnotationSource {
	t.Helper()
	src, err := parse(strings.NewReader(strings.Join(lines, "\n")+"\n"), nil)
	if err != nil {
		t.Fatal(err)
	}
	return src
}

// twoTranscripts is a gene G on strand with a non-coding transcript nc whose
// exon covers [100,200) and a coding transcript cd whose exons cover
// [300,400) and [500,600) and whose CDS is [320,550). Position 150 is exonic
// in nc and, for cd, upstream of its CDS without being in one of its exons.
func twoTranscripts(strand, nc, cd string) []string {
	g := `gene_id "G"; transcript_id "`
	return []string{
		"chr1\tt\texon\t101\t200\t.\t" + strand + "\t.\t" + g + nc + `";`,
		"chr1\tt\texon\t301\t400\t.\t" + strand + "\t.\t" + g + cd + `";`,
		"chr1\tt\tCDS\t321\t400\t.\t" + strand + "\t0\t" + g + cd + `";`,
		"chr1\tt\texon\t501\t600\t.\t" + strand + "\t.\t" + g + cd + `";`,
		"chr1\tt\tCDS\t501\t550\t.\t" + strand + "\t0\t" + g + cd + `";`,
	}
}

// The region must not depend on what the transcripts are called. It used to:
// visited coding-first (sorted IDs "A_cd" < "Z_nc"), position 150 was marked a
// 5'UTR intron before the exon hit arrived and came out a non-coding exon;
// visited the other way it was a 5'UTR.
func TestRegionDoesNotDependOnTranscriptNames(t *testing.T) {
	for _, tc := range []struct {
		strand string
		want   GenicRegion
	}{
		{"+", UTR5}, // upstream of the CDS on +
		{"-", UTR3}, // the same coordinates are downstream on -
	} {
		for _, names := range [][2]string{{"Z_nc", "A_cd"}, {"A_nc", "Z_cd"}} {
			src := parseGTF(t, twoTranscripts(tc.strand, names[0], names[1])...)
			for _, geneID := range []string{"", "G"} {
				if got := src.FindGenicRegionForPos("chr1", 150, bed.StrandNone, geneID); got != tc.want {
					t.Errorf("strand %s, transcripts %v, geneID %q: got %s, want %s",
						tc.strand, names, geneID, got.Code, tc.want.Code)
				}
			}
		}
	}
}

// Across genes too: a non-coding gene's exon at a position upstream of a
// coding gene's CDS makes it UTR whichever gene is visited first.
func TestRegionDoesNotDependOnGeneNames(t *testing.T) {
	for _, ids := range [][2]string{{"A_nc", "Z_cd"}, {"Z_nc", "A_cd"}} {
		src := parseGTF(t,
			"chr1\tt\texon\t101\t200\t.\t+\t.\tgene_id \""+ids[0]+"\"; transcript_id \"nc1\";",
			"chr1\tt\texon\t51\t80\t.\t+\t.\tgene_id \""+ids[1]+"\"; transcript_id \"cd1\";",
			"chr1\tt\texon\t301\t400\t.\t+\t.\tgene_id \""+ids[1]+"\"; transcript_id \"cd1\";",
			"chr1\tt\tCDS\t321\t400\t.\t+\t0\tgene_id \""+ids[1]+"\"; transcript_id \"cd1\";",
		)
		if got := src.FindGenicRegionForPos("chr1", 150, bed.StrandPlus, ""); got != UTR5 {
			t.Errorf("genes %v: got %s, want 5_utr", ids, got.Code)
		}
	}
}

// Exonic in any transcript is exonic, and coding in any is coding.
func TestExonicInAnyTranscript(t *testing.T) {
	g := `gene_id "G"; transcript_id "`
	src := parseGTF(t,
		// T1 covers [100,400) with one exon; T2 skips [200,300) as an intron.
		"chr1\tt\texon\t101\t400\t.\t+\t.\t"+g+`T1";`,
		"chr1\tt\texon\t101\t200\t.\t+\t.\t"+g+`T2";`,
		"chr1\tt\texon\t301\t400\t.\t+\t.\t"+g+`T2";`,
		// T3 is coding over [120,180), inside T2's first exon.
		"chr1\tt\texon\t101\t200\t.\t+\t.\t"+g+`T3";`,
		"chr1\tt\tCDS\t121\t180\t.\t+\t0\t"+g+`T3";`,
	)
	for _, tc := range []struct {
		pos  int
		want GenicRegion
	}{
		{250, UTR3}, // T2's intron, but T1's exon; past T3's CDS
		{150, Coding},
		{110, UTR5},
	} {
		if got := src.FindGenicRegionForPos("chr1", tc.pos, bed.StrandNone, ""); got != tc.want {
			t.Errorf("pos %d: got %s, want %s", tc.pos, got.Code, tc.want.Code)
		}
	}
}

// With no gene ID the indexed source classifies against every gene at the
// position, as the in-memory one does. It used to look up a gene with an empty
// ID and return Intergenic everywhere.
func TestIndexedWithoutGeneIDMatchesInMemory(t *testing.T) {
	dir := t.TempDir()
	mem, err := NewAnnotationSource(writePlainFixture(t, dir), nil)
	if err != nil {
		t.Fatal(err)
	}
	idx, err := NewIndexedAnnotationSource(writeIndexedFixture(t, filepath.Join(dir)), nil)
	if err != nil {
		t.Fatal(err)
	}
	defer idx.Close()
	genic := 0
	for pos := 900; pos < 6100; pos += 7 {
		for _, st := range []bed.Strand{bed.StrandNone, bed.StrandPlus, bed.StrandMinus} {
			want := mem.FindGenicRegionForPos("chr1", pos, st, "")
			if got := idx.FindGenicRegionForPos("chr1", pos, st, ""); got != want {
				t.Fatalf("chr1:%d strand %v: indexed %s, in-memory %s", pos, st, got.Code, want.Code)
			}
			if want.IsGene {
				genic++
			}
		}
	}
	if genic == 0 {
		t.Fatal("no genic positions compared; the fixture did not load")
	}
}
