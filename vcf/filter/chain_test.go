package filter

import (
	"errors"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/vcf"
)

// vhdr carries the INFO/FORMAT declarations the value filters read. The base
// hdr in filter_test.go declares only GT.
const vhdr = "##fileformat=VCFv4.2\n" +
	`##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">` + "\n" +
	`##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">` + "\n" +
	`##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">` + "\n" +
	`##INFO=<ID=CSQ,Number=1,Type=String,Description="Consequence">` + "\n" +
	`##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">` + "\n" +
	`##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">` + "\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"

func valueRecs(t *testing.T, lines ...string) (*vcf.VcfHeader, []*vcf.VcfRecord) {
	t.Helper()
	r, err := vcf.NewVcfReader(strings.NewReader(vhdr + strings.Join(lines, "\n") + "\n"))
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
		if err != nil {
			break
		}
		out = append(out, rec)
	}
	return h, out
}

// Chain is the composition mechanism every consumer goes through, and it had no
// test: filter_test.go exercised only the individual filters.
func TestChainAppliesEveryFilterInOrder(t *testing.T) {
	h, recs := valueRecs(t,
		"chr1\t100\t.\tA\tG\t50\t.\tDP=5;AF=0.9\tGT:GQ\t0/1:15",
	)
	rec := recs[0]

	// "INFO" addresses the INFO column; "" would mean "any sample", and a
	// sample name addresses that one. See valueFilter's doc.
	fs := []Filter{
		NewLessThan("DP", 10, "INFO", ""),
		NewGreaterThan("AF", 0.5, "INFO", ""),
		NewLessThan("GQ", 20, "S1", ""),
	}
	c := NewChain()
	if c.Len() != 0 {
		t.Fatalf("a new chain has %d filters, want 0", c.Len())
	}
	for _, f := range fs {
		c.Add(f)
	}
	if c.Len() != 3 {
		t.Fatalf("Len = %d, want 3", c.Len())
	}

	if err := c.SetupHeaders(h); err != nil {
		t.Fatal(err)
	}
	// Every filter declared itself, so the output VCF can explain its own codes.
	for _, f := range fs {
		if _, ok := h.FilterDef(f.ID()); !ok {
			t.Errorf("no ##FILTER declared for %q", f.ID())
		}
	}

	if err := c.Apply(rec); err != nil {
		t.Fatal(err)
	}
	// All three fail this record, and all three stamp. A chain does not stop at
	// the first failure -- the FILTER column is meant to accumulate reasons.
	if got := rec.Filters(); len(got) != 3 {
		t.Errorf("Filters() = %v, want all three codes", got)
	}
	if err := c.Close(); err != nil {
		t.Fatal(err)
	}
}

func TestChainLeavesPassingRecordsAlone(t *testing.T) {
	h, recs := valueRecs(t, "chr1\t100\t.\tA\tG\t50\t.\tDP=50;AF=0.1\tGT:GQ\t0/1:99")

	c := NewChain()
	c.Add(NewLessThan("DP", 10, "INFO", ""))
	c.Add(NewGreaterThan("AF", 0.5, "INFO", ""))
	if err := c.SetupHeaders(h); err != nil {
		t.Fatal(err)
	}
	if err := c.Apply(recs[0]); err != nil {
		t.Fatal(err)
	}
	if got := recs[0].Filters(); len(got) != 0 {
		t.Errorf("a passing record was stamped %v", got)
	}
}

// An empty chain is a no-op rather than an error, which is what lets a CLI build
// one unconditionally and only then decide whether any filters were requested.
func TestEmptyChainIsANoOp(t *testing.T) {
	h, recs := valueRecs(t, "chr1\t100\t.\tA\tG\t50\t.\tDP=50\tGT:GQ\t0/1:99")
	c := NewChain()
	if err := c.SetupHeaders(h); err != nil {
		t.Fatal(err)
	}
	if err := c.Apply(recs[0]); err != nil {
		t.Fatal(err)
	}
	if got := recs[0].Filters(); len(got) != 0 {
		t.Errorf("empty chain stamped %v", got)
	}
	if err := c.Close(); err != nil {
		t.Fatal(err)
	}
}

// errFilter fails at each stage in turn, so the chain's propagation can be
// checked without a real filter that happens to be able to fail.
type errFilter struct {
	id                            string
	setupErr, filterErr, closeErr error
	filtered                      *bool
}

func (e *errFilter) ID() string                       { return e.id }
func (e *errFilter) SetupHeader(*vcf.VcfHeader) error { return e.setupErr }
func (e *errFilter) Close() error                     { return e.closeErr }
func (e *errFilter) Filter(*vcf.VcfRecord) error {
	if e.filtered != nil {
		*e.filtered = true
	}
	return e.filterErr
}

func TestChainPropagatesErrors(t *testing.T) {
	boom := errors.New("boom")
	h, recs := valueRecs(t, "chr1\t100\t.\tA\tG\t50\t.\tDP=50\tGT:GQ\t0/1:99")

	t.Run("SetupHeaders", func(t *testing.T) {
		c := NewChain()
		c.Add(&errFilter{id: "bad", setupErr: boom})
		if err := c.SetupHeaders(h); !errors.Is(err, boom) {
			t.Errorf("SetupHeaders err = %v, want boom", err)
		}
	})

	t.Run("Apply stops at the first failure", func(t *testing.T) {
		reached := false
		c := NewChain()
		c.Add(&errFilter{id: "bad", filterErr: boom})
		c.Add(&errFilter{id: "after", filtered: &reached})
		if err := c.Apply(recs[0]); !errors.Is(err, boom) {
			t.Errorf("Apply err = %v, want boom", err)
		}
		if reached {
			t.Error("a filter after the failing one still ran")
		}
	})

	t.Run("Close closes all and reports the first error", func(t *testing.T) {
		second := errors.New("second")
		closed := false
		c := NewChain()
		c.Add(&errFilter{id: "bad", closeErr: boom})
		c.Add(&errFilter{id: "also", closeErr: second, filtered: &closed})
		if err := c.Close(); !errors.Is(err, boom) {
			t.Errorf("Close err = %v, want the first error", err)
		}
	})
}

// One constructor per comparison family in value.go, which had no test file at
// all. Each case names the record it should and should not stamp; a filter that
// silently never fires is the failure mode worth catching.
func TestValueFilterFamilies(t *testing.T) {
	cases := []struct {
		name   string
		filter Filter
		line   string
		stamp  bool
	}{
		{"less than fires", NewLessThan("DP", 10, "INFO", ""), "DP=5", true},
		{"less than passes", NewLessThan("DP", 10, "INFO", ""), "DP=50", false},
		{"less than is strict at the threshold", NewLessThan("DP", 10, "INFO", ""), "DP=10", false},
		{"less-or-equal at the threshold", NewLessThanEqual("DP", 10, "INFO", ""), "DP=10", true},

		{"greater than fires", NewGreaterThan("DP", 10, "INFO", ""), "DP=50", true},
		{"greater than passes", NewGreaterThan("DP", 10, "INFO", ""), "DP=5", false},
		{"greater-or-equal at the threshold", NewGreaterThanEqual("DP", 10, "INFO", ""), "DP=10", true},

		{"equals fires", NewEquals("CSQ", "missense", "INFO", ""), "CSQ=missense", true},
		{"equals passes", NewEquals("CSQ", "missense", "INFO", ""), "CSQ=synonymous", false},
		{"not-equals fires", NewNotEquals("CSQ", "missense", "INFO", ""), "CSQ=synonymous", true},
		{"not-equals passes", NewNotEquals("CSQ", "missense", "INFO", ""), "CSQ=missense", false},

		{"contains fires", NewContains("CSQ", "mis", "INFO", ""), "CSQ=missense", true},
		{"contains passes", NewContains("CSQ", "mis", "INFO", ""), "CSQ=synonymous", false},
		{"not-contains fires", NewNotContains("CSQ", "mis", "INFO", ""), "CSQ=synonymous", true},

		{"in-list fires", NewInList("CSQ", []string{"a", "missense"}, "INFO", ""), "CSQ=missense", true},
		{"in-list passes", NewInList("CSQ", []string{"a", "b"}, "INFO", ""), "CSQ=missense", false},
		{"not-in-list fires", NewNotInList("CSQ", []string{"a", "b"}, "INFO", ""), "CSQ=missense", true},

		{"flag present fires", NewFlagPresent("DB"), "DB", true},
		{"flag present passes", NewFlagPresent("DB"), "DP=5", false},
		{"flag absent fires", NewFlagAbsent("DB"), "DP=5", true},
		{"flag absent passes", NewFlagAbsent("DB"), "DB", false},

		{"value missing fires", NewValueMissing("DP", "INFO"), "CSQ=x", true},
		{"value missing passes", NewValueMissing("DP", "INFO"), "DP=5", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			h, recs := valueRecs(t, "chr1\t100\t.\tA\tG\t50\t.\t"+c.line+"\tGT:GQ\t0/1:99")
			got := apply(t, c.filter, h, recs[0])
			if stamped := has(got, c.filter.ID()); stamped != c.stamp {
				t.Errorf("%s stamped=%v (codes %v), want %v", c.filter.ID(), stamped, got, c.stamp)
			}
		})
	}
}

// A per-sample filter reads that sample's FORMAT value, not the INFO column.
func TestValueFilterReadsNamedSample(t *testing.T) {
	h, recs := valueRecs(t, "chr1\t100\t.\tA\tG\t50\t.\tDP=99\tGT:GQ\t0/1:5")
	f := NewLessThan("GQ", 20, "S1", "")
	got := apply(t, f, h, recs[0])
	if !has(got, f.ID()) {
		t.Errorf("a sample GQ of 5 should fail --min GQ 20; codes %v", got)
	}
}
