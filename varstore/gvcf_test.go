package varstore

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/htsio/tabix"
)

// gvcfLines is a single-sample gVCF covering every case that matters:
//
//	chr1  100..5000   reference block, MIN_DP 28
//	chr1  5001        a real variant, ALT "G,<NON_REF>" as GATK writes it
//	chr1  5002..9000  reference block, MIN_DP 25
//	chr1  9001..11999 NOTHING -- a gap. Never sequenced, must stay unanswerable.
//	chr1  12000..13000 reference block, MIN_DP 3 -- below a typical gate
//	chr1  20000..20500 reference block whose genotype is ./. -- a no-call, not an
//	                   observation of reference
var gvcfLines = []string{
	"##fileformat=VCFv4.2",
	"##contig=<ID=chr1,length=100000>",
	`##ALT=<ID=NON_REF,Description="Represents any possible alternative allele">`,
	`##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">`,
	`##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`,
	`##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">`,
	`##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">`,
	`##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP in the block">`,
	`##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="Reference genotype confidence">`,
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
	"chr1\t100\t.\tA\t<NON_REF>\t.\t.\tEND=5000\tGT:DP:MIN_DP:RGQ\t0/0:40:28:60",
	"chr1\t5001\trs1\tA\tG,<NON_REF>\t500\tPASS\t.\tGT:DP:GQ\t0/1:44:99",
	"chr1\t5002\t.\tC\t<NON_REF>\t.\t.\tEND=9000\tGT:DP:MIN_DP:RGQ\t0/0:30:25:55",
	"chr1\t12000\t.\tT\t<NON_REF>\t.\t.\tEND=13000\tGT:DP:MIN_DP:RGQ\t0/0:9:3:20",
	"chr1\t20000\t.\tG\t<NON_REF>\t.\t.\tEND=20500\tGT:DP:MIN_DP:RGQ\t./.:0:0:0",
}

// writeGvcf writes the fixture twice: plain, and bgzipped with a tabix index. Both
// are returned because correctness must not depend on the index -- an unindexed gVCF
// scans, which is slower and must not be *different*.
func writeGvcf(t *testing.T) (plain, indexed string) {
	t.Helper()
	dir := t.TempDir()
	plain = filepath.Join(dir, "s.g.vcf")
	if err := os.WriteFile(plain, []byte(strings.Join(gvcfLines, "\n")+"\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	indexed = filepath.Join(dir, "s.g.vcf.gz")
	w := tabix.NewWriter(indexed, tabix.NewWriterOpts().VCF().AutoIndex())
	for _, l := range gvcfLines {
		if strings.HasPrefix(l, "#") {
			w.WriteHeader(l)
			continue
		}
		if err := w.Write(l); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	return plain, indexed
}

func rowsOf(t *testing.T, path string, q Query) []Call {
	t.Helper()
	s, err := OpenVcf(path)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	out, err := CollectCalls(s, q)
	if err != nil {
		t.Fatalf("Calls: %v", err)
	}
	return out
}

// TestGvcfBlockAnswersInsideItsSpan is the core of gVCF support: a position no variant
// record mentions is still answerable, because a block asserted coverage there.
// A span query needs the index -- that is pre-existing (seek hands the span to
// scan, which requires it; only the locus path degrades to a full scan), so this
// covers the indexed path and TestGvcfExactLocusInsideBlockMatches covers the
// unindexed one through the same block logic.
func TestGvcfBlockAnswersInsideItsSpan(t *testing.T) {
	_, indexed := writeGvcf(t)
	for _, path := range []string{indexed} {
		t.Run("indexed", func(t *testing.T) {
			got := rowsOf(t, path, Query{
				Spans:      []Span{{Chrom: "chr1", Start: 1999, End: 2100}},
				IncludeRef: true,
			})
			if len(got) != 1 {
				t.Fatalf("want exactly one row for a position inside one block, got %d: %+v", len(got), got)
			}
			c := got[0]
			if c.Alt == "<NON_REF>" {
				t.Errorf("block reported as carrying <NON_REF>; that is a genotype the sample does not have")
			}
			if c.Alt != "." {
				t.Errorf("Alt = %q, want \".\" -- a block names no alternate", c.Alt)
			}
			if c.GT != "0/0" {
				t.Errorf("GT = %q, want 0/0", c.GT)
			}
			if c.Pos != 100 || c.RefEnd != 5000 {
				t.Errorf("Pos/RefEnd = %d/%d, want 100/5000 (the block's extent)", c.Pos, c.RefEnd)
			}
			if c.MinDP != 28 {
				t.Errorf("MinDP = %d, want 28 (the block's MIN_DP)", c.MinDP)
			}
			if c.DP != Missing {
				t.Errorf("DP = %d, want Missing -- a block records no depth at any single base", c.DP)
			}
			if c.GQ != 60 {
				t.Errorf("GQ = %d, want 60 (the block's RGQ)", c.GQ)
			}
		})
	}
}

// TestGvcfGapIsUnanswerable is the boundary. A gVCF says a great deal, but not
// everything: between blocks nothing was reported, and inventing reference there
// would be exactly the error the whole sites-catalog rule exists to prevent.
func TestGvcfGapIsUnanswerable(t *testing.T) {
	_, indexed := writeGvcf(t)
	got := rowsOf(t, indexed, Query{
		Spans:      []Span{{Chrom: "chr1", Start: 9999, End: 10100}},
		IncludeRef: true,
	})
	if len(got) != 0 {
		t.Errorf("a gap between blocks must yield nothing, got %+v", got)
	}
}

// TestGvcfGateUsesBlockMinDP checks the gate is not silently fail-open. A block-derived
// call has no DP, and Missing passes every gate, so gating had to move to MIN_DP -- if
// it did not, a block covered at depth 3 would pass --min-dp 10.
func TestGvcfGateUsesBlockMinDP(t *testing.T) {
	_, indexed := writeGvcf(t)
	span := []Span{{Chrom: "chr1", Start: 12499, End: 12600}}

	ungated := rowsOf(t, indexed, Query{Spans: span, IncludeRef: true})
	if len(ungated) != 1 {
		t.Fatalf("ungated: want the block at 12000, got %+v", ungated)
	}
	if ungated[0].MinDP != 3 {
		t.Fatalf("ungated MinDP = %d, want 3", ungated[0].MinDP)
	}

	gated := rowsOf(t, indexed, Query{Spans: span, IncludeRef: true, Gate: Gate{MinDP: 10}})
	if len(gated) != 0 {
		t.Errorf("a block whose MIN_DP is 3 must not pass --min-dp 10, got %+v", gated)
	}

	// And a well-covered block still passes, so the gate is not simply rejecting
	// every block.
	ok := rowsOf(t, indexed, Query{
		Spans:      []Span{{Chrom: "chr1", Start: 1999, End: 2100}},
		IncludeRef: true, Gate: Gate{MinDP: 10},
	})
	if len(ok) != 1 {
		t.Errorf("a block with MIN_DP 28 must pass --min-dp 10, got %+v", ok)
	}
}

// TestGvcfNoCallBlockIsNotReference: a block whose genotype is ./. states that nothing
// was called, not that the sample is reference.
func TestGvcfNoCallBlockIsNotReference(t *testing.T) {
	_, indexed := writeGvcf(t)
	got := rowsOf(t, indexed, Query{
		Spans:      []Span{{Chrom: "chr1", Start: 20099, End: 20200}},
		IncludeRef: true,
	})
	if len(got) != 0 {
		t.Errorf("a ./. block is not a reference observation, got %+v", got)
	}
}

// TestGvcfVariantKeepsRealAlleleOnly: GATK writes the block allele beside the real one.
// The variant must be reported and the block allele must not.
func TestGvcfVariantKeepsRealAlleleOnly(t *testing.T) {
	_, indexed := writeGvcf(t)
	got := rowsOf(t, indexed, Query{
		Loci: []Locus{{Chrom: "chr1", Pos: 5001, Ref: "A", Alt: "G"}},
	})
	if len(got) != 1 {
		t.Fatalf("want the G call at 5001, got %+v", got)
	}
	if got[0].Alt != "G" || got[0].GT != "0/1" {
		t.Errorf("got Alt=%q GT=%q, want G and 0/1", got[0].Alt, got[0].GT)
	}

	// Nothing anywhere may report <NON_REF> as an allele.
	all := rowsOf(t, indexed, Query{IncludeRef: true})
	for _, c := range all {
		if c.Alt == "<NON_REF>" || c.Alt == "<*>" {
			t.Errorf("row reports a block allele as an alternate: %+v", c)
		}
	}
}

// TestGvcfExactLocusInsideBlockMatches pins the positional-matching rule. A block has
// no ALT, so identity matching cannot find it -- a query naming chr1:2000:A:T is asking
// whether the sample is reference at 2000, and the block is the answer.
// Also the unindexed path: a locus query degrades to a full scan, so this proves the
// block logic does not depend on the index -- only its speed does.
func TestGvcfExactLocusInsideBlockMatches(t *testing.T) {
	plain, indexed := writeGvcf(t)
	for _, path := range []string{plain, indexed} {
		name := "unindexed"
		if strings.HasSuffix(path, ".gz") {
			name = "indexed"
		}
		t.Run(name, func(t *testing.T) {
			got := rowsOf(t, path, Query{
				Loci:       []Locus{{Chrom: "chr1", Pos: 2000, Ref: "A", Alt: "T"}},
				IncludeRef: true,
			})
			if len(got) != 1 || got[0].Pos != 100 {
				t.Fatalf("an exact locus inside a block must be answered by it, got %+v", got)
			}
			if got[0].MinDP != 28 || got[0].Alt != "." {
				t.Errorf("got MinDP=%d Alt=%q, want 28 and \".\"", got[0].MinDP, got[0].Alt)
			}
		})
	}
}

// TestGvcfClassifyUsesBlockCoverage is the state-level counterpart: a position inside a
// block is an *observed* non-carrier, which is the whole point of a gVCF. A plain VCF
// can only say not-assayed there.
func TestGvcfClassifyUsesBlockCoverage(t *testing.T) {
	_, indexed := writeGvcf(t)
	s, err := OpenVcf(indexed)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	for _, tc := range []struct {
		name  string
		locus Locus
		gate  Gate
		want  State
	}{
		{"inside a well-covered block", Locus{"chr1", 2000, "A", "T"}, Gate{MinDP: 10}, StateNonCarrier},
		{"in a gap between blocks", Locus{"chr1", 10000, "A", "T"}, Gate{MinDP: 10}, StateNotAssayed},
		{"inside a block below the gate", Locus{"chr1", 12500, "A", "T"}, Gate{MinDP: 10}, StateNotAssayed},
		{"inside that block, ungated", Locus{"chr1", 12500, "A", "T"}, Gate{}, StateNonCarrier},
		{"inside a no-call block", Locus{"chr1", 20100, "A", "T"}, Gate{}, StateNotAssayed},
		// An explicit record outranks a block: the sample really does carry G here,
		// and the block starting at 5002 must not turn that into a non-carrier.
		{"at a real variant", Locus{"chr1", 5001, "A", "G"}, Gate{MinDP: 10}, StateCarrier},
	} {
		t.Run(tc.name, func(t *testing.T) {
			states, err := s.Classify(tc.locus, tc.gate)
			if err != nil {
				t.Fatal(err)
			}
			if len(states) != 1 {
				t.Fatalf("want one sample, got %d", len(states))
			}
			if states[0].State != tc.want {
				t.Errorf("Classify(%v) = %s, want %s", tc.locus, states[0].State, tc.want)
			}
		})
	}
}

// TestGvcfSiteKnownCountsBlockCoverage: SiteKnown asks "did the source look here",
// and a block covering the position means it did. Getting this wrong makes
// vcf-varquery warn that a position is absent from a file that in fact reported it.
func TestGvcfSiteKnownCountsBlockCoverage(t *testing.T) {
	_, indexed := writeGvcf(t)
	s, err := OpenVcf(indexed)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	for _, tc := range []struct {
		name  string
		locus Locus
		want  bool
	}{
		{"inside a block", Locus{"chr1", 2000, "A", "T"}, true},
		{"in a gap", Locus{"chr1", 10000, "A", "T"}, false},
		{"at a real variant", Locus{"chr1", 5001, "A", "G"}, true},
		{"past every record", Locus{"chr1", 90000, "A", "T"}, false},
	} {
		t.Run(tc.name, func(t *testing.T) {
			got, err := s.SiteKnown(tc.locus)
			if err != nil {
				t.Fatal(err)
			}
			if got != tc.want {
				t.Errorf("SiteKnown(%v) = %v, want %v", tc.locus, got, tc.want)
			}
		})
	}
}
