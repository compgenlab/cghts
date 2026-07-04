package annotate

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/compgenlab/hts/vcf"
)

// annGTF: GeneA (+, coding) over [100,200) and GeneB (-, non-coding lncRNA)
// over [150,250) overlap at [150,200), so a variant there hits both genes.
const annGTF = `chr1	t	exon	101	200	.	+	.	gene_id "GeneA"; gene_name "GeneA"; transcript_id "TA"; gene_type "protein_coding";
chr1	t	CDS	101	200	.	+	0	gene_id "GeneA"; gene_name "GeneA"; transcript_id "TA"; gene_type "protein_coding";
chr1	t	exon	151	250	.	-	.	gene_id "GeneB"; gene_name "GeneB"; transcript_id "TB"; gene_type "lncRNA";
`

func writeAnnGTF(t *testing.T) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), "ann.gtf")
	if err := os.WriteFile(path, []byte(annGTF), 0o644); err != nil {
		t.Fatalf("write gtf: %v", err)
	}
	return path
}

func info(t *testing.T, rec *vcf.VcfRecord, key string) (string, bool) {
	t.Helper()
	v, ok := rec.Info().Get(key)
	if !ok {
		return "", false
	}
	return v.String(), true
}

func TestGtfAnnotator(t *testing.T) {
	h, recs := bedRecs(t,
		"chr1\t120\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",  // GeneA only (coding exon)
		"chr1\t170\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",  // GeneA + GeneB overlap
		"chr1\t9000\t.\tA\tG\t.\tPASS\t.\tGT\t0/1") // intergenic

	a, err := NewGtfAnnotator(GtfOptions{Filename: writeAnnGTF(t)})
	if err != nil {
		t.Fatalf("NewGtfAnnotator: %v", err)
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

	// Record 0: GeneA only.
	checkInfo(t, recs[0], "GTF_GENE", "GeneA")
	checkInfo(t, recs[0], "GTF_GENEID", "GeneA")
	checkInfo(t, recs[0], "GTF_STRAND", "+")
	checkInfo(t, recs[0], "GTF_REGION", "coding_exon")
	checkInfo(t, recs[0], "GTF_BIOTYPE", "protein_coding")
	checkInfo(t, recs[0], "GTF_CODING", "GeneA")
	if v, ok := info(t, recs[0], "GTF_NONCODING"); ok {
		t.Errorf("rec0 CG_NONCODING = %q, want absent", v)
	}

	// Record 1: GeneA (coding exon) + GeneB (nc exon), comma-joined in start order.
	checkInfo(t, recs[1], "GTF_GENE", "GeneA,GeneB")
	checkInfo(t, recs[1], "GTF_STRAND", "+,-")
	checkInfo(t, recs[1], "GTF_REGION", "coding_exon,nc_exon")
	checkInfo(t, recs[1], "GTF_BIOTYPE", "protein_coding,lncRNA")
	checkInfo(t, recs[1], "GTF_CODING", "GeneA")
	checkInfo(t, recs[1], "GTF_NONCODING", "GeneB")

	// Record 2: intergenic — no GTF INFO at all.
	if v, ok := info(t, recs[2], "GTF_GENE"); ok {
		t.Errorf("rec2 CG_GENE = %q, want absent (intergenic)", v)
	}
}

func TestGtfAnnotatorHeader(t *testing.T) {
	h, _ := bedRecs(t)
	a, err := NewGtfAnnotator(GtfOptions{Filename: writeAnnGTF(t), Prefix: "GT_"})
	if err != nil {
		t.Fatalf("NewGtfAnnotator: %v", err)
	}
	defer a.Close()
	if err := a.SetupHeader(h); err != nil {
		t.Fatalf("SetupHeader: %v", err)
	}
	for _, id := range []string{"GT_GENE", "GT_GENEID", "GT_STRAND", "GT_BIOTYPE", "GT_REGION", "GT_CODING", "GT_NONCODING"} {
		if _, ok := h.InfoDef(id); !ok {
			t.Errorf("missing ##INFO def for %s", id)
		}
	}
}

func checkInfo(t *testing.T, rec *vcf.VcfRecord, key, want string) {
	t.Helper()
	got, ok := info(t, rec, key)
	if !ok {
		t.Errorf("%s missing, want %q", key, want)
		return
	}
	if got != want {
		t.Errorf("%s = %q, want %q", key, got, want)
	}
}
