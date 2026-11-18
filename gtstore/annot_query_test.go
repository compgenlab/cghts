package gtstore

import (
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/vcf"
)

// A cohort with typed and imputed records, a Number=A INFO field missing for
// one allele, and a FORMAT field captured on calls. Every record carries DP so
// the two backends' reference calls agree under the gate.
const annotVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=IMPUTED,Number=0,Type=Flag,Description="imputed">
##INFO=<ID=R2,Number=1,Type=Float,Description="r2">
##INFO=<ID=AF,Number=A,Type=Float,Description="af">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="DS">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	100	.	A	G	.	PASS	R2=0.870;AF=0.25	GT:DP:DS	0/1:30:1.000	0/0:30:0	0/0:5:0.01
chr1	200	.	C	T	.	PASS	IMPUTED;R2=0.5;AF=0.4	GT:DP:DS	1/1:30:2	0/0:30:0	0/1:30:1
chr1	300	.	G	A,C	.	PASS	AF=0.1,.;R2=0.95	GT:DP:DS	1/2:30:0.9,0.6	0/0:30:0,0	0/1:30:0.8,0
`

var (
	annotInfo = []InfoField{
		InfoFieldFor("IMPUTED", "Flag", "0"),
		InfoFieldFor("R2", "Float", "1"),
		InfoFieldFor("AF", "Float", "A"),
	}
	annotFormat = []FormatField{FormatFieldFor("DS", "Float", "A")}
	annotGate   = Gate{MinDP: 10}
)

func writeAnnotVCF(t *testing.T) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), "cohort.vcf")
	if err := os.WriteFile(p, []byte(annotVCF), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

// convertForTest is a minimal converter: sites with their captured INFO, ALT
// calls with their captured FORMAT, and one catalog run per called site. It
// reads values through InfoValue/FormatValue, the reading conversion shares.
func convertForTest(t *testing.T, vcfPath string, blocks []RefBlock) string {
	t.Helper()
	base := filepath.Join(t.TempDir(), "store")
	r, err := vcf.OpenVcfFile(context.Background(), vcfPath)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	h, err := r.Header()
	if err != nil {
		t.Fatal(err)
	}
	w, err := NewWriter(base, WriterOpts{
		Samples: h.Samples(), MinDP: 10, Info: annotInfo, Format: annotFormat,
		RefBlocks: len(blocks) > 0,
	})
	if err != nil {
		t.Fatal(err)
	}
	for {
		rec, err := r.NextRecord()
		if err != nil {
			break
		}
		sr := NewSampleReader(rec)
		for j, alt := range rec.Alt() {
			info := map[string]any{}
			for _, f := range annotInfo {
				if v, ok := InfoValue(rec, f, j); ok {
					info[f.Name] = v
				}
			}
			if err := w.WriteSiteInfo(Site{Chrom: rec.Chrom, Pos: int32(rec.Pos), Ref: rec.Ref, Alt: alt}, info); err != nil {
				t.Fatal(err)
			}
			for i, name := range h.Samples() {
				c, ok := CallFor(rec, name, sr.Read(i), j+1, alt)
				if !ok {
					continue
				}
				fv := map[string]any{}
				for _, f := range annotFormat {
					if v, ok := FormatValue(rec, i, f, j); ok {
						fv[f.Name] = v
					}
				}
				if err := w.WriteCallFormat(c, fv); err != nil {
					t.Fatal(err)
				}
			}
		}
		for i, name := range h.Samples() {
			sf := sr.Read(i)
			if HasCall(sf.GT) && sf.DP != Missing && sf.DP >= 10 {
				if err := w.WriteRegion(CalledSiteRun{SampleID: name, Chrom: rec.Chrom,
					Start: int32(rec.Pos), End: int32(rec.Pos), NSites: 1, MinDP: sf.DP}); err != nil {
					t.Fatal(err)
				}
			}
		}
	}
	for _, b := range blocks {
		if err := w.WriteRefBlock(b); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
	return base
}

// rowText renders what both backends must agree on. Reference rows leave out
// FORMAT and the quality fields: a store keeps none for a genotype it never
// wrote down, an asymmetry that predates captured values.
func rowText(c Call) string {
	s := fmt.Sprintf("%s:%d:%s:%s %s %s", c.Chrom, c.Pos, c.Ref, c.Alt, c.SampleID, c.GT)
	if c.Extra != nil {
		s += " info{" + mapText(c.Extra.Info) + "}"
		if c.GT != "0/0" {
			s += " fmt{" + mapText(c.Extra.Format) + "}"
		}
	}
	return s
}

func mapText(m map[string]string) string {
	keys := make([]string, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	var b strings.Builder
	for _, k := range keys {
		fmt.Fprintf(&b, "%s=%s;", k, m[k])
	}
	return b.String()
}

func queryRows(t *testing.T, s Store, q Query) []string {
	t.Helper()
	calls, err := CollectCalls(s, q)
	if err != nil {
		t.Fatalf("Calls(%+v): %v", q, err)
	}
	out := make([]string, len(calls))
	for i, c := range calls {
		out[i] = rowText(c)
	}
	return out
}

func TestBackendsAgreeOnInfoFilterAndValues(t *testing.T) {
	path := writeAnnotVCF(t)
	vs, err := OpenVcf(path)
	if err != nil {
		t.Fatal(err)
	}
	ps, err := OpenParquet(convertForTest(t, path, nil))
	if err != nil {
		t.Fatal(err)
	}
	defer ps.Close()

	cases := map[string]Query{
		"plain":            {},
		"exclude imputed":  {ExcludeInfo: []string{"IMPUTED"}},
		"require imputed":  {RequireInfo: []string{"IMPUTED"}},
		"require AF":       {RequireInfo: []string{"AF"}},
		"show info":        {WithInfo: []string{"*"}},
		"show format":      {WithFormat: []string{"DS"}},
		"hom-ref + filter": {IncludeRef: true, ExcludeInfo: []string{"IMPUTED"}},
		"hom-ref + values": {IncludeRef: true, WithInfo: []string{"R2", "AF"}, WithFormat: []string{"*"}},
		"locus + exclude":  {Loci: []Locus{{Chrom: "chr1", Pos: 200, Ref: "C", Alt: "T"}}, IncludeRef: true, ExcludeInfo: []string{"IMPUTED"}},
	}
	for name, q := range cases {
		t.Run(name, func(t *testing.T) {
			q.Gate = annotGate
			got, want := queryRows(t, ps, q), queryRows(t, vs, q)
			if strings.Join(got, "\n") != strings.Join(want, "\n") {
				t.Errorf("store and VCF disagree\nstore:\n  %s\nvcf:\n  %s",
					strings.Join(got, "\n  "), strings.Join(want, "\n  "))
			}
		})
	}

	// The values themselves: "0.870" in the VCF and 0.87 in the store are one
	// value, rendered one way; a Number=A field takes its allele's element, and
	// "." for the second allele is absent rather than "." or zero.
	rows := queryRows(t, ps, Query{WithInfo: []string{"*"}, WithFormat: []string{"DS"}, Gate: annotGate})
	for _, want := range []string{
		"chr1:100:A:G S1 0/1 info{AF=0.25;R2=0.87;} fmt{DS=1;}",
		"chr1:200:C:T S1 1/1 info{AF=0.4;IMPUTED=1;R2=0.5;} fmt{DS=2;}",
		"chr1:300:G:A S1 1/. info{AF=0.1;R2=0.95;} fmt{DS=0.9;}",
		"chr1:300:G:C S1 ./1 info{R2=0.95;} fmt{DS=0.6;}",
	} {
		if !contains(rows, want) {
			t.Errorf("missing %q in:\n  %s", want, strings.Join(rows, "\n  "))
		}
	}
	// Excluding IMPUTED removes the site entirely, ALT calls included.
	for _, r := range queryRows(t, ps, Query{ExcludeInfo: []string{"IMPUTED"}, IncludeRef: true, Gate: annotGate}) {
		if strings.HasPrefix(r, "chr1:200:") {
			t.Errorf("excluded site still answered: %s", r)
		}
	}
}

func contains(rows []string, want string) bool {
	for _, r := range rows {
		if r == want {
			return true
		}
	}
	return false
}

// A site the filter excludes is off the catalog, so an exact locus can still
// be answered 0/0 from a reference block -- except for a sample the store has
// as a carrier there, where the store and the gVCF would disagree.
func TestExcludedSiteAnsweredByRefBlocks(t *testing.T) {
	path := writeAnnotVCF(t)
	blocks := []RefBlock{
		{SampleID: "S1", Chrom: "chr1", Start: 150, End: 250, MinDP: 20}, // a carrier in the store
		{SampleID: "S2", Chrom: "chr1", Start: 150, End: 250, MinDP: 20},
		{SampleID: "S3", Chrom: "chr1", Start: 150, End: 199, MinDP: 20}, // stops short
	}
	ps, err := OpenParquet(convertForTest(t, path, blocks))
	if err != nil {
		t.Fatal(err)
	}
	defer ps.Close()
	locus := Locus{Chrom: "chr1", Pos: 200, Ref: "C", Alt: "T"}

	got := queryRows(t, ps, Query{Loci: []Locus{locus}, IncludeRef: true, Gate: annotGate,
		ExcludeInfo: []string{"IMPUTED"}, WithInfo: []string{"R2"}})
	want := []string{"chr1:200:C:T S2 0/0 info{R2=0.5;}"}
	if strings.Join(got, "\n") != strings.Join(want, "\n") {
		t.Errorf("excluded locus answered as\n  %s\nwant\n  %s", strings.Join(got, "\n  "), strings.Join(want, "\n  "))
	}

	// A span never enumerates off-catalog positions, so there it is simply gone.
	for _, r := range queryRows(t, ps, Query{Spans: []Span{{Chrom: "chr1", Start: 150, End: 250}},
		IncludeRef: true, Gate: annotGate, ExcludeInfo: []string{"IMPUTED"}}) {
		t.Errorf("span query answered an excluded site: %s", r)
	}
}

// A field the source never recorded cannot filter or be shown: an error, not
// "absent everywhere", which would keep or drop every site in silence.
func TestUncapturedFieldIsAnError(t *testing.T) {
	path := writeAnnotVCF(t)
	vs, _ := OpenVcf(path)
	ps, err := OpenParquet(convertForTest(t, path, nil))
	if err != nil {
		t.Fatal(err)
	}
	defer ps.Close()
	for name, s := range map[string]Store{"store": ps, "vcf": vs} {
		for _, q := range []Query{
			{ExcludeInfo: []string{"NOPE"}},
			{RequireInfo: []string{"imputed"}}, // keys are case-sensitive, as VCF is
			{WithInfo: []string{"NOPE"}},
			{WithFormat: []string{"PL"}},
		} {
			if _, err := s.Calls(q); !errors.Is(err, ErrNotCaptured) {
				t.Errorf("%s: Calls(%+v) = %v, want ErrNotCaptured", name, q, err)
			}
		}
	}
}

func TestValueText(t *testing.T) {
	for v, want := range map[any]string{
		int32(7): "7", 0.87: "0.87", 1.0: "1", 1e-7: "1e-07", "a,b": "a,b", true: "1",
	} {
		if got := ValueText(v); got != want {
			t.Errorf("ValueText(%v) = %q, want %q", v, got, want)
		}
	}
}

// A panel read from a filtered callset can be limited to the variants that
// passed.
func TestParseTargetsPassingOnly(t *testing.T) {
	p := filepath.Join(t.TempDir(), "panel.vcf")
	body := "##fileformat=VCFv4.2\n##FILTER=<ID=LowQual,Description=\"x\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t100\t.\tA\tG\t.\tPASS\t.\n" +
		"chr1\t200\t.\tC\tT\t.\tLowQual\t.\n" +
		"chr1\t300\t.\tG\tA\t.\t.\t.\n"
	if err := os.WriteFile(p, []byte(body), 0o644); err != nil {
		t.Fatal(err)
	}
	all, err := ParseTargets(context.Background(), []string{p})
	if err != nil {
		t.Fatal(err)
	}
	passing, err := ParseTargetsWith(context.Background(), []string{p}, TargetOptions{Passing: true})
	if err != nil {
		t.Fatal(err)
	}
	if len(all.Loci) != 3 || all.Filtered != 0 {
		t.Errorf("all: %d loci, %d filtered; want 3, 0", len(all.Loci), all.Filtered)
	}
	if len(passing.Loci) != 2 || passing.Filtered != 1 {
		t.Errorf("passing: %d loci, %d filtered; want 2, 1", len(passing.Loci), passing.Filtered)
	}
	for _, l := range passing.Loci {
		if l.Pos == 200 {
			t.Error("the LowQual record was kept")
		}
	}
}
