package varstore

import (
	"context"
	"math"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"
)

// Parser-level tests for ParseTargets. The behaviour through a query -- that a
// BED selects 0-based, that an unmatched target is reported -- is tested where it
// is used, in the consuming CLIs.

func writeTargetFile(t *testing.T, name, content string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

// TestTargetFormatDetection pins the rule that decides what a target file is.
//
// A VCF announces itself. Otherwise column 3 decides, because a BED's third
// column is an end coordinate and always numeric where a site list's is a REF
// allele and never is.
func TestTargetFormatDetection(t *testing.T) {
	cases := []struct {
		name    string
		content string
		want    TargetFormat
	}{
		{"bed.txt", "chr1\t99\t300\n", TargetBED},
		{"bed-extra-cols.txt", "chr1 99 300 name 0 +\n", TargetBED},
		{"bed-with-track.txt", "track name=x\nchr1 99 300\n", TargetBED},
		{"sites.txt", "chr1 100 A G\n", TargetSiteList},
		{"sites-extra-cols.txt", "chr1 100 A G rs123\n", TargetSiteList},
		{"sites-pos-only.txt", "chr1 100\n", TargetSiteList},
		{"sites-colon.txt", "chr1:100:A:G\n", TargetSiteList},
		{"sites-commented.txt", "# a comment\n\nchr1 100 A G\n", TargetSiteList},
		{"panel.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\n", TargetVCF},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got, err := DetectTargetFile(context.Background(), writeTargetFile(t, tc.name, tc.content))
			if err != nil {
				t.Fatalf("detect: %v", err)
			}
			if got != tc.want {
				t.Errorf("detected %q, want %q", got, tc.want)
			}
		})
	}
}

// The inline grammar, including the one ambiguity it has to settle: a contig
// name carrying colons splits like a locus, and only the numeric tail tells
// them apart.
func TestParseSelector(t *testing.T) {
	cases := []struct {
		in    string
		locus *Locus
		span  *Span
	}{
		{"chr1:1000:A:T", &Locus{Chrom: "chr1", Pos: 1000, Ref: "A", Alt: "T"}, nil},
		{"chr1:1000", nil, &Span{Chrom: "chr1", Start: 999, End: 1000}},
		{"chr1:1000-2000", nil, &Span{Chrom: "chr1", Start: 999, End: 2000}},
		{"chr1", nil, &Span{Chrom: "chr1", Start: 0, End: math.MaxInt32}},
		{"HLA-A*01:01:01:01", nil, &Span{Chrom: "HLA-A*01:01:01:01", Start: 0, End: math.MaxInt32}},
		{"HLA-B*07:02:01", nil, &Span{Chrom: "HLA-B*07:02:01", Start: 0, End: math.MaxInt32}},
	}
	for _, c := range cases {
		l, sp, err := ParseSelector(c.in)
		if err != nil {
			t.Errorf("%s: %v", c.in, err)
			continue
		}
		if !reflect.DeepEqual(l, c.locus) || !reflect.DeepEqual(sp, c.span) {
			t.Errorf("%s: got %+v %+v, want %+v %+v", c.in, l, sp, c.locus, c.span)
		}
	}
}

// Order records how targets were given, interleaving loci and spans, which is
// what a command reporting per target in turn relies on.
func TestParseTargetsKeepsOrderAcrossFilesAndInline(t *testing.T) {
	list := writeTargetFile(t, "sites.txt", "# panel\nchr1 100 A G\nchr2 50\n")
	got, err := ParseTargets(context.Background(), []string{"chr3:7:C:T", list, "chr4:1-10"})
	if err != nil {
		t.Fatal(err)
	}
	var order []string
	for _, r := range got.Order {
		if r.Locus {
			order = append(order, got.Loci[r.Index].String())
		} else {
			sp := got.Spans[r.Index]
			order = append(order, sp.Chrom)
		}
	}
	want := []string{"chr3:7:C:T", "chr1:100:A:G", "chr2", "chr4"}
	if !reflect.DeepEqual(order, want) {
		t.Errorf("order %v, want %v", order, want)
	}
	if got.Counts[TargetInline] != 2 || got.Counts[TargetSiteList] != 2 {
		t.Errorf("counts %v", got.Counts)
	}
	var b strings.Builder
	got.Report(&b)
	if !strings.Contains(b.String(), "read as site list") {
		t.Errorf("report does not name the detected format:\n%s", b.String())
	}
}

func TestParseTargetsCountsVcfRecordsWithNoAlt(t *testing.T) {
	vcfPath := writeTargetFile(t, "panel.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"+
		"chr1\t100\t.\tA\tG,T\t.\t.\t.\nchr1\t200\t.\tC\t.\t.\t.\t.\n")
	got, err := ParseTargets(context.Background(), []string{vcfPath})
	if err != nil {
		t.Fatal(err)
	}
	if len(got.Loci) != 2 || got.Skipped != 1 {
		t.Errorf("loci %v, skipped %d; want one locus per ALT and the '.' record counted", got.Loci, got.Skipped)
	}
}
