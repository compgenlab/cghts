package vcf

import (
	"strings"
	"testing"

	"github.com/compgenlab/cghts/internal/vcfspan"
)

// spanCases are lines whose reference span is interesting. Shared by the two tests
// below so the parsed and line-oriented paths are held to the same corpus.
var spanCases = []struct {
	name string
	line string
	want int // 0-based exclusive end, for a record at POS 1000 (beg 999)
}{
	{"SNV", "chr1\t1000\t.\tA\tT\t.\tPASS\t.", 1000},
	{"plain deletion spans len(REF)", "chr1\t1000\t.\t" + strings.Repeat("A", 200) + "\tA\t.\tPASS\t.", 1199},
	{"insertion spans only its REF", "chr1\t1000\t.\tA\tACCCCC\t.\tPASS\t.", 1000},
	{"symbolic DEL via INFO/END", "chr1\t1000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=2000", 2000},
	{"INFO/END not past POS ignored", "chr1\t1000\t.\tACGT\tA\t.\tPASS\tEND=500", 1003},
	{"SVLEN negative for a deletion", "chr1\t1000\t.\tN\t<DEL>\t.\tPASS\tSVLEN=-750", 1750},
	{"SVLEN on an insertion does not widen", "chr1\t1000\t.\tN\t<INS>\t.\tPASS\tSVLEN=5000", 1000},
	{"subtyped deletion", "chr1\t1000\t.\tN\t<DEL:ME:ALU>\t.\tPASS\tSVLEN=-300", 1300},
	{"FORMAT/LEN max across samples", "chr1\t1000\t.\tA\t<NON_REF>\t.\t.\t.\tGT:LEN\t0/0:100\t0/0:450\t0/0:20", 1449},
	{"FORMAT/LEN with <*>", "chr1\t1000\t.\tA\t<*>\t.\t.\t.\tGT:DP:LEN\t0/0:30:75", 1074},
	{"FORMAT/LEN ignored without a ref-block ALT", "chr1\t1000\t.\tA\tT\t.\tPASS\t.\tGT:LEN\t0/1:9000", 1000},
	{"multiallelic takes the widest", "chr1\t1000\t.\tACGTACGT\tA,AC\t.\tPASS\t.", 1007},
}

// TestRefSpanMatchesLineFields is the drift guard. The tabix index writer computes a
// record's span from a split line; everything above it computes the same span from a
// parsed record. If those two ever disagree, an index says a record reaches 2 kb
// while a reader says one base, queries silently return the wrong rows, and no test
// anywhere would notice -- so hold both to one corpus.
func TestRefSpanMatchesLineFields(t *testing.T) {
	for _, tc := range spanCases {
		t.Run(tc.name, func(t *testing.T) {
			rec, err := newRecord(tc.line, nil)
			if err != nil {
				t.Fatalf("newRecord: %v", err)
			}
			start, end := rec.RefSpan()
			if start != 999 {
				t.Fatalf("RefSpan start = %d, want 999", start)
			}
			if end != tc.want {
				t.Errorf("(*VcfRecord).RefSpan() end = %d, want %d", end, tc.want)
			}
			if fe := vcfspan.FieldsEnd(strings.Split(tc.line, "\t"), 999); fe != end {
				t.Errorf("parsed path gave %d but line path gave %d -- the two have drifted", end, fe)
			}
		})
	}
}

func TestIsRefBlock(t *testing.T) {
	for _, tc := range []struct {
		line string
		want bool
	}{
		{"chr1\t1000\t.\tA\t<NON_REF>\t.\t.\t.", true},
		{"chr1\t1000\t.\tA\t<*>\t.\t.\t.", true},
		{"chr1\t1000\t.\tA\t.\t.\t.\t.", true},
		{"chr1\t1000\t.\tA\tT\t.\t.\t.", false},
		{"chr1\t1000\t.\tN\t<DEL>\t.\t.\t.", false},
		// GATK emits the real ALT alongside <NON_REF> at a variant site. That record
		// describes a variant, so it must NOT read as a reference block -- otherwise
		// a caller listing variants silently drops every call in a gVCF.
		{"chr1\t1000\t.\tA\tT,<NON_REF>\t.\t.\t.", false},
		{"chr1\t1000\t.\tA\t<NON_REF>,<*>\t.\t.\t.", true},
	} {
		rec, err := newRecord(tc.line, nil)
		if err != nil {
			t.Fatalf("newRecord(%q): %v", tc.line, err)
		}
		if got := rec.IsRefBlock(); got != tc.want {
			t.Errorf("IsRefBlock(%q) = %v, want %v", tc.line, got, tc.want)
		}
	}
}
