package tabix

import (
	"path/filepath"
	"strings"
	"testing"
)

// TestVcfSpanEnd covers the four sources a VCF record's end can come from, and
// the cases where a source must NOT widen the interval.
func TestVcfSpanEnd(t *testing.T) {
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
			got := vcfSpanEnd(strings.Split(tc.line, "\t"), beg)
			if got != tc.want {
				t.Errorf("vcfSpanEnd() = %d, want %d", got, tc.want)
			}
		})
	}
}

// TestVcfPresetRecordedInIndex pins the header field. Without it a reader cannot
// tell that VCF rules apply, so the end derivation silently would not run --
// which is how the whole class of bug went unnoticed.
func TestVcfPresetRecordedInIndex(t *testing.T) {
	o := NewWriterOpts().VCF()
	if got := o.preset & presetMask; got != PresetVCF {
		t.Errorf("VCF() preset = %d, want %d", got, PresetVCF)
	}
	if o.colEnd != 0 {
		t.Errorf("VCF() colEnd = %d, want 0 (the end is derived, not read)", o.colEnd)
	}
	if got := NewWriterOpts().BED().preset & presetMask; got != PresetGeneric {
		t.Errorf("BED() preset = %d, want generic", got)
	}
}

// TestQueryFindsSpanningRecords is the end-to-end case: a region query strictly
// inside a long record must return it. Before the span derivation, both the
// binning and the overlap test treated every VCF record as one base wide, so
// these queries came back empty.
func TestQueryFindsSpanningRecords(t *testing.T) {
	path := filepath.Join(t.TempDir(), "sv.vcf.gz")

	w := NewWriter(path, NewWriterOpts().VCF().AutoIndex())
	w.WriteHeader("##fileformat=VCFv4.2")
	w.WriteHeader("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1")
	lines := []string{
		// a symbolic deletion spanning 1000..2000
		"chr1\t1000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=2000\tGT\t0/1",
		// an explicit-REF deletion spanning 3000..3200
		"chr1\t3000\t.\t" + strings.Repeat("A", 200) + "\tA\t.\tPASS\t.\tGT\t0/1",
		// a gVCF-style reference block, length carried per sample
		"chr1\t4000\t.\tA\t<NON_REF>\t.\t.\t.\tGT:LEN\t0/0:5000",
		// a plain SNV well past all of them, as a control
		"chr1\t20000\t.\tG\tC\t.\tPASS\t.\tGT\t0/1",
		// A 100kb block, whose span crosses many 16kb bin windows. This is the
		// case the WRITER's binning is needed for: the records above all sit in
		// the same window as the queries that find them, so a correct overlap
		// test alone suffices for those. Here the query's bin cannot hold a
		// record binned at its start, so the chunk is never even examined.
		"chr1\t100000\t.\tA\t<NON_REF>\t.\t.\tEND=200000\tGT\t0/0",
	}
	for _, l := range lines {
		if err := w.Write(l); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	r, err := NewReader(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()

	// 0-based half-open query windows, each landing strictly inside a record and
	// well past its first base.
	for _, tc := range []struct {
		name     string
		beg, end int
		wantPos  string
	}{
		{"inside a symbolic <DEL> via INFO/END", 1400, 1600, "1000"},
		{"inside a 200bp explicit-REF deletion", 3100, 3150, "3000"},
		{"inside a gVCF block via FORMAT/LEN", 8000, 8100, "4000"},
		{"a plain SNV still works", 19999, 20100, "20000"},
		// Past the end of every span: the derivation must not widen records
		// without bound, or queries would return everything.
		{"inside a 100kb block, crossing bin windows", 150000, 150100, "100000"},
		{"beyond every span returns nothing", 500000, 500100, ""},
	} {
		t.Run(tc.name, func(t *testing.T) {
			it, err := r.Query("chr1", tc.beg, tc.end)
			if err != nil {
				t.Fatal(err)
			}
			var got []string
			for rec, err := range it {
				if err != nil {
					t.Fatal(err)
				}
				got = append(got, strings.Split(rec.Line, "\t")[1])
			}
			if tc.wantPos == "" {
				if len(got) != 0 {
					t.Errorf("query [%d,%d) = %v, want no records", tc.beg, tc.end, got)
				}
				return
			}
			if len(got) != 1 || got[0] != tc.wantPos {
				t.Errorf("query [%d,%d) = %v, want [%s]", tc.beg, tc.end, got, tc.wantPos)
			}
		})
	}
}
