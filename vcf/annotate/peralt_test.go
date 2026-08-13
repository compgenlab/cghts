package annotate

import (
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/htsio/tabix"
	"github.com/compgenlab/cghts/vcf"
)

// multiAltSource writes a tabix-indexed source VCF from raw data lines.
//
// Built per test rather than shared, because what these assert is how the
// *order* of the source's records affects the output — so each needs to control
// that order itself.
func multiAltSource(t *testing.T, lines ...string) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), "src.vcf.gz")
	w := tabix.NewWriter(path, tabix.NewWriterOpts().VCF().AutoIndex())
	w.WriteHeader("##fileformat=VCFv4.2")
	w.WriteHeader(`##INFO=<ID=AF,Number=1,Type=Float,Description="af">`)
	w.WriteHeader("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
	for _, l := range lines {
		if err := w.Write(l); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	return path
}

// annotatePerAlt runs one exact-match annotation over a single record and returns
// the INFO value it wrote.
func annotatePerAlt(t *testing.T, src string, rec *vcf.VcfRecord) string {
	t.Helper()
	a, err := NewVcfAnnotation(VcfOptions{
		Name: "KAF", Field: "AF", Filename: src, Exact: true, Type: TypeFloat,
	})
	if err != nil {
		t.Fatal(err)
	}
	defer a.Close()
	if err := a.SetupHeader(vcf.NewVcfHeader()); err != nil {
		t.Fatal(err)
	}
	if err := a.Annotate(rec); err != nil {
		t.Fatal(err)
	}
	v, ok := rec.Info().Get("KAF")
	if !ok {
		return ""
	}
	return v.String()
}

// A value lands on the allele it describes, whatever order the source is in.
//
// This is the bug these tests exist for. Values were appended in the order the
// source returned them and joined with commas — the separator VCF uses for
// per-allele lists — so the field read as a per-allele vector and was not one.
// With the source listing C before A, a reader assigned C's value to A.
func TestAValueLandsOnItsOwnAllele(t *testing.T) {
	rec := func() *vcf.VcfRecord { return vcf.NewRecord("chr1", 500, "G", []string{"A", "C"}) }

	// The record's ALTs are A,C. The source lists them the other way round.
	reversed := multiAltSource(t,
		"chr1\t500\t.\tG\tC\t.\t.\tAF=0.9",
		"chr1\t500\t.\tG\tA\t.\t.\tAF=0.1")
	if got := annotatePerAlt(t, reversed, rec()); got != "0.1,0.9" {
		t.Errorf("source order C,A gave %q, want \"0.1,0.9\" — the values follow "+
			"the source's order, not the record's ALT order", got)
	}

	// And in the same order, which used to be right only by luck.
	inOrder := multiAltSource(t,
		"chr1\t500\t.\tG\tA\t.\t.\tAF=0.1",
		"chr1\t500\t.\tG\tC\t.\t.\tAF=0.9")
	if got := annotatePerAlt(t, inOrder, rec()); got != "0.1,0.9" {
		t.Errorf("source order A,C gave %q, want \"0.1,0.9\"", got)
	}
}

// An allele with no match is a "." holding its place.
//
// Without the placeholder a single value is written bare, and a positional
// reader assigns it to the first allele — so a value belonging to ALT 2 is
// reported as ALT 1's. Well-formed, undetectable, wrong.
func TestAnUnmatchedAlleleHoldsItsPlace(t *testing.T) {
	// Only the SECOND ALT is in the source.
	src := multiAltSource(t, "chr1\t500\t.\tG\tC\t.\t.\tAF=0.9")
	rec := vcf.NewRecord("chr1", 500, "G", []string{"A", "C"})
	if got := annotatePerAlt(t, src, rec); got != ".,0.9" {
		t.Errorf("got %q, want \".,0.9\" — a bare value reads as the first allele's", got)
	}

	// Only the first.
	src = multiAltSource(t, "chr1\t500\t.\tG\tA\t.\t.\tAF=0.1")
	rec = vcf.NewRecord("chr1", 500, "G", []string{"A", "C"})
	if got := annotatePerAlt(t, src, rec); got != "0.1,." {
		t.Errorf("got %q, want \"0.1,.\"", got)
	}

	// The middle of three.
	src = multiAltSource(t, "chr1\t600\t.\tT\tA\t.\t.\tAF=0.5")
	rec = vcf.NewRecord("chr1", 600, "T", []string{"G", "A", "C"})
	if got := annotatePerAlt(t, src, rec); got != ".,0.5,." {
		t.Errorf("got %q, want \".,0.5,.\"", got)
	}
}

// Nothing matching writes no field, rather than a row of dots.
func TestNoMatchWritesNoField(t *testing.T) {
	src := multiAltSource(t, "chr1\t500\t.\tG\tT\t.\t.\tAF=0.3")
	rec := vcf.NewRecord("chr1", 500, "G", []string{"A", "C"})
	if got := annotatePerAlt(t, src, rec); got != "" {
		t.Errorf("got %q, want no field at all", got)
	}
}

// A single-allele record is unchanged by any of this.
func TestASingleAlleleRecordIsUnaffected(t *testing.T) {
	src := multiAltSource(t, "chr1\t100\t.\tA\tG\t.\t.\tAF=0.2")
	rec := vcf.NewRecord("chr1", 100, "A", []string{"G"})
	if got := annotatePerAlt(t, src, rec); got != "0.2" {
		t.Errorf("got %q, want \"0.2\"", got)
	}
}

// Two source records answering for one allele keep both values, joined with a
// separator that is not the one holding the alleles apart.
//
// Dropping one would be the alternative, and it loses data silently — with
// which value survived depending on the order the source returned them, which is
// the same instability this whole change removes.
func TestTwoValuesForOneAlleleAreBothKept(t *testing.T) {
	src := multiAltSource(t,
		"chr1\t500\t.\tG\tA\t.\t.\tAF=0.1",
		"chr1\t500\t.\tG\tA\t.\t.\tAF=0.2",
		"chr1\t500\t.\tG\tC\t.\t.\tAF=0.9")
	rec := vcf.NewRecord("chr1", 500, "G", []string{"A", "C"})
	got := annotatePerAlt(t, src, rec)
	if !strings.HasPrefix(got, "0.1&0.2") {
		t.Errorf("got %q, want the allele's two values joined with & and kept "+
			"apart from the next allele's", got)
	}
	if !strings.HasSuffix(got, ",0.9") {
		t.Errorf("got %q, want the second allele's value after the comma", got)
	}
}

// The header says Number=A for an exact match and Number=. for a position one.
//
// It said 1 for both while writing a comma-separated list, which is invalid on
// its own terms quite apart from the misattribution — a strict reader is
// entitled to reject it.
func TestTheNumberMatchesWhatIsWritten(t *testing.T) {
	for _, tc := range []struct {
		exact bool
		want  string
	}{
		{true, "A"},  // one value per ALT
		{false, "."}, // a list for the locus; the count is not known in advance
	} {
		a, err := NewVcfAnnotation(VcfOptions{
			Name: "KAF", Field: "AF", Exact: tc.exact,
			Filename: multiAltSource(t, "chr1\t500\t.\tG\tA\t.\t.\tAF=0.1"),
		})
		if err != nil {
			t.Fatal(err)
		}
		h := vcf.NewVcfHeader()
		if err := a.SetupHeader(h); err != nil {
			t.Fatal(err)
		}
		a.Close()
		def, ok := h.InfoDef("KAF")
		if !ok {
			t.Fatal("no ##INFO def")
		}
		if def.Number != tc.want {
			t.Errorf("exact=%v: Number = %q, want %q", tc.exact, def.Number, tc.want)
		}
	}
}

// A position match is not per-allele — it asked about the locus — so it keeps
// list semantics and does not gain placeholders.
func TestAPositionMatchIsNotPerAllele(t *testing.T) {
	src := multiAltSource(t, "chr1\t500\t.\tG\tT\t.\t.\tAF=0.3")
	a, err := NewVcfAnnotation(VcfOptions{
		Name: "KAF", Field: "AF", Filename: src, Exact: false,
	})
	if err != nil {
		t.Fatal(err)
	}
	defer a.Close()
	if err := a.SetupHeader(vcf.NewVcfHeader()); err != nil {
		t.Fatal(err)
	}
	// The source's ALT (T) is not one of the record's, but a position match does
	// not care.
	rec := vcf.NewRecord("chr1", 500, "G", []string{"A", "C"})
	if err := a.Annotate(rec); err != nil {
		t.Fatal(err)
	}
	v, _ := rec.Info().Get("KAF")
	if v.String() != "0.3" {
		t.Errorf("got %q, want \"0.3\" with no per-allele padding", v.String())
	}
}

// The grouped annotator attributes per allele too.
//
// It is a separate implementation of the same rules — one reader and one query
// shared by N fields — so it had the same bug and needs the same proof. A fix
// applied to one and not the other is the kind that looks complete.
func TestTheGroupAttributesPerAlleleToo(t *testing.T) {
	src := multiAltSource(t,
		"chr1\t500\t.\tG\tC\t.\t.\tAF=0.9", // the source lists C first
		"chr1\t500\t.\tG\tA\t.\t.\tAF=0.1")

	g, err := NewVcfAnnotationGroup(VcfGroupOptions{
		Filename: src,
		Fields: []VcfFieldOptions{
			{Name: "KAF", Field: "AF", Exact: true, Type: TypeFloat},
		},
	})
	if err != nil {
		t.Fatal(err)
	}
	defer g.Close()
	h := vcf.NewVcfHeader()
	if err := g.SetupHeader(h); err != nil {
		t.Fatal(err)
	}
	if def, _ := h.InfoDef("KAF"); def.Number != "A" {
		t.Errorf("Number = %q, want A", def.Number)
	}

	rec := vcf.NewRecord("chr1", 500, "G", []string{"A", "C"})
	if err := g.Annotate(rec); err != nil {
		t.Fatal(err)
	}
	v, _ := rec.Info().Get("KAF")
	if v.String() != "0.1,0.9" {
		t.Errorf("got %q, want \"0.1,0.9\" — the group followed the source's order", v.String())
	}

	// And a half-matched record pads rather than shifting.
	only := multiAltSource(t, "chr1\t500\t.\tG\tC\t.\t.\tAF=0.9")
	g2, err := NewVcfAnnotationGroup(VcfGroupOptions{
		Filename: only,
		Fields:   []VcfFieldOptions{{Name: "KAF", Field: "AF", Exact: true}},
	})
	if err != nil {
		t.Fatal(err)
	}
	defer g2.Close()
	if err := g2.SetupHeader(vcf.NewVcfHeader()); err != nil {
		t.Fatal(err)
	}
	rec2 := vcf.NewRecord("chr1", 500, "G", []string{"A", "C"})
	if err := g2.Annotate(rec2); err != nil {
		t.Fatal(err)
	}
	v2, _ := rec2.Info().Get("KAF")
	if v2.String() != ".,0.9" {
		t.Errorf("got %q, want \".,0.9\"", v2.String())
	}
}
