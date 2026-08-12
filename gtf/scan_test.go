package gtf

import (
	"errors"
	"strings"
	"testing"
)

// A GENCODE-shaped fragment: gene rows plus the transcript/exon rows that repeat
// the same gene_id, which is what makes deduplication the point of the scan.
const scanFixture = `##description: test
chr1	HAVANA	gene	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; gene_type "lncRNA"; gene_name "DDX11L1";
chr1	HAVANA	transcript	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_name "DDX11L1";
chr1	HAVANA	exon	11869	12227	.	+	.	gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_name "DDX11L1";
chr17	HAVANA	gene	7668402	7687550	.	-	.	gene_id "ENSG00000141510.17"; gene_type "protein_coding"; gene_name "TP53";
chrY	HAVANA	gene	2781479	2781479	.	+	.	gene_id "ENSG00000182378.14_PAR_Y"; gene_name "PLCXD1";
`

func TestScanGenesReportsEachGeneOnce(t *testing.T) {
	var got []GeneRef
	if err := ScanGenes(strings.NewReader(scanFixture), func(g GeneRef) error {
		got = append(got, g)
		return nil
	}); err != nil {
		t.Fatalf("ScanGenes: %v", err)
	}

	want := []GeneRef{
		{GeneID: "ENSG00000223972.5", GeneName: "DDX11L1"},
		{GeneID: "ENSG00000141510.17", GeneName: "TP53"},
		{GeneID: "ENSG00000182378.14_PAR_Y", GeneName: "PLCXD1"},
	}
	if len(got) != len(want) {
		t.Fatalf("got %d genes, want %d: %v", len(got), len(want), got)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("gene %d = %+v, want %+v", i, got[i], want[i])
		}
	}
}

// RefSeq writes the symbol under "gene"; a scan that only read gene_name would
// report every RefSeq gene as nameless, and a gene list built against one would
// match nothing at all.
func TestScanGenesReadsTheRefSeqSymbolKey(t *testing.T) {
	const refseq = `chr1	BestRefSeq	gene	101	800	.	+	.	gene_id "GeneR"; gene "GeneR"; gene_biotype "protein_coding";
`
	var got []GeneRef
	if err := ScanGenes(strings.NewReader(refseq), func(g GeneRef) error {
		got = append(got, g)
		return nil
	}); err != nil {
		t.Fatal(err)
	}
	if len(got) != 1 || got[0].GeneName != "GeneR" {
		t.Errorf("got %+v, want the symbol from the \"gene\" attribute", got)
	}
}

// The caller writing rows into a database is the one that fails, and it has to
// be able to stop a scan of a 1.7M-row file rather than be handed the rest.
func TestScanGenesStopsOnTheVisitorsError(t *testing.T) {
	boom := errors.New("no room")
	n := 0
	err := ScanGenes(strings.NewReader(scanFixture), func(GeneRef) error {
		n++
		return boom
	})
	if !errors.Is(err, boom) {
		t.Errorf("ScanGenes = %v, want the visitor's error", err)
	}
	if n != 1 {
		t.Errorf("visited %d genes after an error, want 1", n)
	}
}

// A row with no gene_id has no identity to report. Skipped rather than reported
// as a gene named "", which would land in a membership set and match nothing.
func TestScanGenesSkipsRowsWithNoGeneID(t *testing.T) {
	const odd = "chr1\tt\texon\t1\t2\t.\t+\t.\ttranscript_id \"T1\";\nshort\tline\n"
	n := 0
	if err := ScanGenes(strings.NewReader(odd), func(GeneRef) error { n++; return nil }); err != nil {
		t.Fatal(err)
	}
	if n != 0 {
		t.Errorf("reported %d genes from rows with no gene_id", n)
	}
}

func TestTrimGeneIDVersion(t *testing.T) {
	for in, want := range map[string]string{
		"ENSG00000141510.17":       "ENSG00000141510",
		"ENSG00000141510":          "ENSG00000141510",
		"ENSG00000223972.5":        "ENSG00000223972",
		"ENSG00000182378.14_PAR_Y": "ENSG00000182378_PAR_Y",
		"TP53":                     "TP53",
		"LOC102723897":             "LOC102723897",
		// No digits after the dot: nothing that looks like a version.
		"ODD.NAME": "ODD.NAME",
		"":         "",
	} {
		if got := TrimGeneIDVersion(in); got != want {
			t.Errorf("TrimGeneIDVersion(%q) = %q, want %q", in, got, want)
		}
	}
}

// The property the rule exists for: whatever it does to an id, it does the same
// to the copy of that id somebody typed into a list.
func TestTrimGeneIDVersionIsIdempotent(t *testing.T) {
	for _, id := range []string{"ENSG00000141510.17", "ENSG00000182378.14_PAR_Y", "TP53", "ODD.NAME"} {
		once := TrimGeneIDVersion(id)
		if twice := TrimGeneIDVersion(once); twice != once {
			t.Errorf("TrimGeneIDVersion(%q) = %q, then %q", id, once, twice)
		}
	}
}
