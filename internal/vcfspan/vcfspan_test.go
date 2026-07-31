package vcfspan

import (
	"strings"
	"testing"
)

// TestVcfSpanEnd covers the four sources a VCF record's end can come from, and
// the cases where a source must NOT widen the interval.
func TestEnd(t *testing.T) {
	// beg is 0-based; a record at POS 1000 has beg 999.
	const beg = 999

	tests := []struct {
		name string
		line string
		want int
	}{{
		name: "SNV spans one base",
		line: "chr1\t1000\t.\tA\tT\t.\tPASS\t.",
		want: 1000,
	}, {
		name: "len(REF) is the span of a plain deletion",
		line: "chr1\t1000\t.\t" + strings.Repeat("A", 200) + "\tA\t.\tPASS\t.",
		want: 999 + 200,
	}, {
		name: "an insertion still spans only its REF",
		line: "chr1\t1000\t.\tA\tACCCCCCCCC\t.\tPASS\t.",
		want: 1000,
	}, {
		name: "INFO/END on a symbolic ALT",
		line: "chr1\t1000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=2000",
		want: 2000,
	}, {
		name: "INFO/END as the first key",
		line: "chr1\t1000\t.\tN\t<DEL>\t.\tPASS\tEND=1500;SVTYPE=DEL",
		want: 1500,
	}, {
		// A malformed END must not shrink the record, or a query at the true
		// span would miss it -- htslib warns and ignores.
		name: "INFO/END not past POS is ignored",
		line: "chr1\t1000\t.\tACGT\tA\t.\tPASS\tEND=500",
		want: 999 + 4,
	}, {
		name: "SVLEN, reported negative for a deletion",
		line: "chr1\t1000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-750",
		want: 999 + 750,
	}, {
		// SVLEN measures inserted sequence, not reference span, so it must not
		// widen the interval.
		name: "SVLEN on an insertion does not widen",
		line: "chr1\t1000\t.\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=5000",
		want: 1000,
	}, {
		name: "SVLEN on a subtyped deletion still counts",
		line: "chr1\t1000\t.\tN\t<DEL:ME:ALU>\t.\tPASS\tSVLEN=-300",
		want: 999 + 300,
	}, {
		// The per-sample case: one interval per record, so the widest sample
		// wins. Narrowing would lose records; widening only costs a filtered
		// candidate.
		name: "FORMAT/LEN takes the max across samples",
		line: "chr1\t1000\t.\tA\t<NON_REF>\t.\tPASS\t.\tGT:LEN\t0/0:100\t0/0:450\t0/0:20",
		want: 999 + 450,
	}, {
		name: "FORMAT/LEN with the VCF 4.5 <*> allele",
		line: "chr1\t1000\t.\tA\t<*>\t.\tPASS\t.\tGT:DP:LEN\t0/0:30:75",
		want: 999 + 75,
	}, {
		// LEN is only meaningful for a reference block; a normal record that
		// happens to carry the key must not be widened by it.
		name: "FORMAT/LEN ignored without a reference-block ALT",
		line: "chr1\t1000\t.\tA\tT\t.\tPASS\t.\tGT:LEN\t0/1:9000",
		want: 1000,
	}, {
		name: "INFO/END wins over a shorter len(REF)",
		line: "chr1\t1000\t.\tACGT\t<DEL>\t.\tPASS\tEND=3000",
		want: 3000,
	}, {
		name: "len(REF) wins over a shorter INFO/END",
		line: "chr1\t1000\t.\t" + strings.Repeat("A", 50) + "\tA\t.\tPASS\tEND=1010",
		want: 999 + 50,
	}, {
		name: "truncated record does not panic",
		line: "chr1\t1000",
		want: 1000,
	}, {
		name: "a sample column missing its LEN subfield is skipped",
		line: "chr1\t1000\t.\tA\t<NON_REF>\t.\tPASS\t.\tGT:LEN\t0/0\t0/0:60",
		want: 999 + 60,
	}}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got := FieldsEnd(strings.Split(tc.line, "\t"), beg)
			if got != tc.want {
				t.Errorf("FieldsEnd() = %d, want %d", got, tc.want)
			}
		})
	}
}
