package varstore

import (
	"context"
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestStoreKind(t *testing.T) {
	dir := t.TempDir()
	store := filepath.Join(dir, "cohort")
	buildCensusStore(t, store, WriterOpts{MinDP: 10})

	ctx := context.Background()
	for _, tc := range []struct {
		in, want string
	}{
		{"x.vcf", KindVcf},
		{"x.vcf.gz", KindVcf},
		{"x.VCF.GZ", KindVcf},
		{"x.vcf.bgz", KindVcf},
		{"x.bcf", KindVcf},
		{"https://host/cohort.vcf.gz", KindVcf},
		{"s3://bucket/cohort.vcf.gz", KindVcf},

		{"cohort/calls.parquet", KindParquet},
		{"s3://bucket/cohort/sites.parquet", KindParquet},
		{"s3://bucket/cohort/" + VolumeManifestFile, KindParquet},

		// A bare directory is resolved by finding a manifest in it.
		{store, KindParquet},
		{store + "/", KindParquet},
	} {
		got, err := StoreKind(ctx, tc.in)
		if err != nil {
			t.Errorf("StoreKind(%q): %v", tc.in, err)
			continue
		}
		if got != tc.want {
			t.Errorf("StoreKind(%q) = %q, want %q", tc.in, got, tc.want)
		}
	}
}

// An unlinked transport is not an unrecognized file type. Reporting it as one
// sends a user to recheck their path when the fix is to import a package.
func TestStoreKindSeparatesTransportFromFormat(t *testing.T) {
	_, err := StoreKind(context.Background(), "gs://bucket/cohort")
	if err == nil {
		t.Fatal("gs:// resolved to a store kind")
	}
	if !strings.Contains(err.Error(), "no transport") {
		t.Errorf("error blames the format rather than the transport: %v", err)
	}

	_, err = StoreKind(context.Background(), filepath.Join(t.TempDir(), "mystery"))
	if err == nil {
		t.Fatal("a directory with no manifest resolved to a store kind")
	}
	if !strings.Contains(err.Error(), "cannot tell") {
		t.Errorf("unhelpful error for an unrecognizable target: %v", err)
	}
}

func TestOpenStoreInfersAndForces(t *testing.T) {
	dir := t.TempDir()
	store := filepath.Join(dir, "cohort")
	buildCensusStore(t, store, WriterOpts{MinDP: 10})
	ctx := context.Background()

	s, err := Open(ctx, store, "")
	if err != nil {
		t.Fatal(err)
	}
	if _, ok := s.(*ParquetVolume); !ok {
		t.Errorf("inferred backend is %T, want *ParquetVolume", s)
	}
	s.Close()

	if _, err := Open(ctx, store, "nonsense"); err == nil {
		t.Error("an unknown --store kind was accepted")
	}
}

// The sites a VCF defines, computed on the way past. Counts have to mean the
// same thing they do in a converted store, or the two backends answer the same
// question differently.
func TestVcfStoreSites(t *testing.T) {
	path := filepath.Join(t.TempDir(), "in.vcf")
	body := strings.Join([]string{
		"##fileformat=VCFv4.2",
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
		"chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\t1/1\t0/0",
		"chr1\t200\t.\tC\tG,T\t.\t.\t.\tGT\t0/1\t0/2\t./.",
		"",
	}, "\n")
	if err := os.WriteFile(path, []byte(body), 0o644); err != nil {
		t.Fatal(err)
	}

	s, err := OpenVcf(path)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	var got []Site
	if err := s.Sites(func(site Site) bool {
		got = append(got, site)
		return true
	}); err != nil {
		t.Fatal(err)
	}

	// One entry per alternate, so the multiallelic record contributes two.
	if len(got) != 3 {
		t.Fatalf("got %d sites, want 3: %+v", len(got), got)
	}

	// chr1:100 A>T -- one het and one hom-alt is AC 3 over AN 6, two carriers.
	if g := got[0]; g.Pos != 100 || g.Alt != "T" || g.AC != 3 || g.AN != 6 || g.NCarriers != 2 {
		t.Errorf("chr1:100 = %+v, want AC 3 / AN 6 / 2 carriers", g)
	}
	// chr1:200 -- S3 is ./. so AN counts only the four called alleles, and each
	// alternate is carried once.
	if g := got[1]; g.Alt != "G" || g.AC != 1 || g.AN != 4 || g.NCarriers != 1 {
		t.Errorf("chr1:200 G = %+v, want AC 1 / AN 4 / 1 carrier", g)
	}
	if g := got[2]; g.Alt != "T" || g.AC != 1 || g.AN != 4 || g.NCarriers != 1 {
		t.Errorf("chr1:200 T = %+v, want AC 1 / AN 4 / 1 carrier", g)
	}

	// A --min-dp threshold is not recorded anywhere in a plain VCF, so these
	// must stay zero rather than imply one was applied.
	for _, g := range got {
		if g.NCalled != 0 || g.NLowDP != 0 {
			t.Errorf("%s:%d reports depth-gated counts a VCF cannot know: %+v", g.Chrom, g.Pos, g)
		}
	}
}

// The two failure modes must stay distinguishable: a caller that suggests
// --store for a transport failure sends the user somewhere useless.
func TestUnknownKindIsDistinctFromTransportFailure(t *testing.T) {
	_, err := StoreKind(context.Background(), filepath.Join(t.TempDir(), "mystery"))
	if !errors.Is(err, ErrUnknownStoreKind) {
		t.Errorf("an unrecognizable target is not ErrUnknownStoreKind: %v", err)
	}
	_, err = StoreKind(context.Background(), "gs://bucket/cohort")
	if errors.Is(err, ErrUnknownStoreKind) {
		t.Errorf("an unlinked transport was reported as an unrecognized store: %v", err)
	}
}

// An interrupted conversion leaves tables but no manifest. That has to resolve
// as a store so opening it produces the actionable error, rather than as an
// unrecognizable path -- the manifest is both what identifies a store and the
// thing this one lacks, so treating its absence as "not a store" would make the
// most likely real failure undiagnosable.
func TestUnfinishedStoreIsStillRecognizedAsOne(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	buildCensusStore(t, base, WriterOpts{MinDP: 10})
	if err := os.Remove(VolumeManifestPath(base)); err != nil {
		t.Fatal(err)
	}

	kind, err := StoreKind(context.Background(), base)
	if err != nil {
		t.Fatalf("an unfinished store was not recognized: %v", err)
	}
	if kind != KindParquet {
		t.Errorf("kind = %q, want %q", kind, KindParquet)
	}

	_, err = Open(context.Background(), base, "")
	if err == nil {
		t.Fatal("an unfinished store opened")
	}
	if !strings.Contains(err.Error(), VolumeManifestFile) {
		t.Errorf("the error does not name the missing manifest: %v", err)
	}
}
