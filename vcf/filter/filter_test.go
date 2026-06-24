package filter

import (
	"io"
	"strings"
	"testing"

	"github.com/compgenlab/hts/vcf"
)

const hdr = "##fileformat=VCFv4.2\n" +
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"

func recs(t *testing.T, lines ...string) (*vcf.VcfHeader, []*vcf.VcfRecord) {
	t.Helper()
	r, err := vcf.NewVcfReader(strings.NewReader(hdr + strings.Join(lines, "\n") + "\n"))
	if err != nil {
		t.Fatal(err)
	}
	h, err := r.Header()
	if err != nil {
		t.Fatal(err)
	}
	var out []*vcf.VcfRecord
	for {
		rec, err := r.NextRecord()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatal(err)
		}
		out = append(out, rec)
	}
	return h, out
}

// filters returns the FILTER codes a record carries after applying f.
func apply(t *testing.T, f Filter, h *vcf.VcfHeader, rec *vcf.VcfRecord) []string {
	t.Helper()
	if err := f.SetupHeader(h); err != nil {
		t.Fatalf("SetupHeader: %v", err)
	}
	if err := f.Filter(rec); err != nil {
		t.Fatalf("Filter: %v", err)
	}
	return rec.Filters()
}

func has(codes []string, want string) bool {
	for _, c := range codes {
		if c == want {
			return true
		}
	}
	return false
}

func TestChromFilters(t *testing.T) {
	h, rs := recs(t,
		"chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",
		"chr2\t200\t.\tA\tG\t.\tPASS\t.\tGT\t0/1")
	f := NewChromPass([]string{"chr1"})
	if codes := apply(t, f, h, rs[0]); len(codes) != 0 {
		t.Errorf("chr1 should pass --chrom-pass chr1, got %v", codes)
	}
	if codes := apply(t, f, h, rs[1]); !has(codes, "CHROM_PASS_chr1") {
		t.Errorf("chr2 should be flagged, got %v", codes)
	}
}

func TestSNVIndelFilters(t *testing.T) {
	h, rs := recs(t,
		"chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",  // SNV
		"chr1\t200\t.\tG\tGA\t.\tPASS\t.\tGT\t0/1") // indel
	if codes := apply(t, NewSNV(), h, rs[0]); !has(codes, "SNV") {
		t.Errorf("SNV not flagged: %v", codes)
	}
	if codes := apply(t, NewSNV(), h, rs[1]); has(codes, "SNV") {
		t.Errorf("indel should not be SNV-flagged: %v", codes)
	}
	h2, rs2 := recs(t, "chr1\t200\t.\tG\tGA\t.\tPASS\t.\tGT\t0/1")
	if codes := apply(t, NewIndel(), h2, rs2[0]); !has(codes, "INDEL") {
		t.Errorf("indel not flagged: %v", codes)
	}
}

func TestQualFilter(t *testing.T) {
	h, rs := recs(t,
		"chr1\t100\t.\tA\tG\t20.0\tPASS\t.\tGT\t0/1",
		"chr1\t200\t.\tA\tG\t50.0\tPASS\t.\tGT\t0/1")
	if codes := apply(t, NewQual(30), h, rs[0]); !has(codes, "QUAL_lt_30.0") {
		t.Errorf("low qual not flagged: %v", codes)
	}
	if codes := apply(t, NewQual(30), h, rs[1]); has(codes, "QUAL_lt_30.0") {
		t.Errorf("high qual should not be flagged: %v", codes)
	}
}

func TestHetHomFilters(t *testing.T) {
	h, rs := recs(t,
		"chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1", // het
		"chr1\t200\t.\tA\tG\t.\tPASS\t.\tGT\t1/1") // hom
	if codes := apply(t, NewHeterozygous(), h, rs[0]); !has(codes, "heterozygous") {
		t.Errorf("0/1 not het-flagged: %v", codes)
	}
	if codes := apply(t, NewHomozygous(), h, rs[1]); !has(codes, "homozygous") {
		t.Errorf("1/1 not hom-flagged: %v", codes)
	}
	if codes := apply(t, NewHomozygous(), h, rs[0]); has(codes, "homozygous") {
		t.Errorf("0/1 should not be hom-flagged: %v", codes)
	}
}

func TestMaxInsDel(t *testing.T) {
	h, rs := recs(t,
		"chr1\t100\t.\tA\tACCC\t.\tPASS\t.\tGT\t0/1", // 3bp insertion
		"chr1\t200\t.\tACGT\tA\t.\tPASS\t.\tGT\t0/1") // 3bp deletion
	if codes := apply(t, NewMaxIns(2), h, rs[0]); !has(codes, "INS_max_2") {
		t.Errorf("3bp insertion not flagged at max 2: %v", codes)
	}
	if codes := apply(t, NewMaxDel(2), h, rs[1]); !has(codes, "DEL_max_2") {
		t.Errorf("3bp deletion not flagged at max 2: %v", codes)
	}
}
