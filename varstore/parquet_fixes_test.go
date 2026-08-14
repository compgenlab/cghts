package varstore

import (
	"os"
	"path/filepath"
	"reflect"
	"testing"
)

// The per-chromosome census is keyed on canonical identity, not on spelling. A
// conversion whose inputs mix "chr1" and "1" -- one per-chromosome VCF from each
// of two sources, say -- was producing two entries for one chromosome in the
// field the manifest's own doc calls the part that earns the file.
func TestCensusFoldsChromosomeSpellings(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	w, err := NewWriter(base, WriterOpts{Samples: []string{"S1"}, MinDP: 10})
	if err != nil {
		t.Fatal(err)
	}
	for _, chrom := range []string{"chr1", "1", "NC_000001.11"} {
		if err := w.WriteSite(Site{Chrom: chrom, Pos: 100, Ref: "A", Alt: "T", AC: 1, AN: 2}); err != nil {
			t.Fatal(err)
		}
		if err := w.WriteCall(Call{
			SampleID: "S1", Chrom: chrom, Pos: 100, Ref: "A", Alt: "T", GT: "0/1", DP: 30,
		}); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}

	m, err := ReadVolumeManifest(base)
	if err != nil {
		t.Fatal(err)
	}
	if len(m.Chromosomes) != 1 {
		t.Fatalf("three spellings of chr1 produced %d census entries: %+v",
			len(m.Chromosomes), m.Chromosomes)
	}
	c := m.Chromosomes[0]
	// The first spelling seen is kept: a store records the naming its source
	// used rather than rewriting it.
	if c.Name != "chr1" {
		t.Errorf("census name = %q, want the first spelling seen (chr1)", c.Name)
	}
	if c.Sites != 3 || c.Calls != 3 {
		t.Errorf("census = %d sites / %d calls, want 3 and 3", c.Sites, c.Calls)
	}
}

// Rows for a call whose site is absent from the catalog are emitted last -- and
// in a stable order. They used to be drained straight out of a Go map, so the
// same query against the same store returned them differently on each run, and
// anyone diffing two runs saw a difference that was not in the data.
func TestIncludeRefTailOrderIsDeterministic(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	samples := []string{"S1", "S2", "S3"}
	w, err := NewWriter(base, WriterOpts{Samples: samples, MinDP: 10})
	if err != nil {
		t.Fatal(err)
	}
	// One site in the catalog, and calls at several loci that are not. The
	// off-catalog ones are what lands in the tail.
	if err := w.WriteSite(Site{Chrom: "chr1", Pos: 100, Ref: "A", Alt: "T", AC: 1, AN: 6}); err != nil {
		t.Fatal(err)
	}
	for _, pos := range []int32{100, 200, 300, 400, 500} {
		for _, s := range samples {
			if err := w.WriteCall(Call{
				SampleID: s, Chrom: "chr1", Pos: pos, Ref: "A", Alt: "T", GT: "0/1", DP: 30,
			}); err != nil {
				t.Fatal(err)
			}
		}
	}
	if err := w.WriteRegion(CalledSiteRun{
		SampleID: "S1", Chrom: "chr1", Start: 100, End: 500, NSites: 5,
	}); err != nil {
		t.Fatal(err)
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}

	var first []Call
	for i := 0; i < 12; i++ {
		s, err := OpenParquet(base)
		if err != nil {
			t.Fatal(err)
		}
		got, err := CollectCalls(s, Query{IncludeRef: true})
		s.Close()
		if err != nil {
			t.Fatal(err)
		}
		if i == 0 {
			first = got
			continue
		}
		if !reflect.DeepEqual(got, first) {
			t.Fatalf("run %d returned a different row order:\n got %+v\nwant %+v", i, got, first)
		}
	}
	if len(first) == 0 {
		t.Fatal("no rows; the test is not exercising the tail")
	}
}

// Finish leaves nothing when it cannot write the manifest. By that point Close
// has given all three tables valid footers, so what is on disk is
// indistinguishable from a finished conversion missing only its marker -- and a
// reader would refuse it forever, with the overwrite guard blocking the retry.
func TestFinishDiscardsWhenTheManifestCannotBeWritten(t *testing.T) {
	dir := t.TempDir()
	base := filepath.Join(dir, "cohort")
	w, err := NewWriter(base, WriterOpts{Samples: []string{"S1"}, MinDP: 10})
	if err != nil {
		t.Fatal(err)
	}
	if err := w.WriteSite(Site{Chrom: "chr1", Pos: 100, Ref: "A", Alt: "T", AC: 1, AN: 2}); err != nil {
		t.Fatal(err)
	}

	// Remove a table out from under the writer. The handle stays open, so Close
	// still writes its footer, but manifest()'s os.Stat then fails -- which is
	// the failure mode this guards.
	if err := os.Remove(SitesPath(base)); err != nil {
		t.Fatal(err)
	}

	if err := w.Finish(); err == nil {
		t.Fatal("Finish succeeded despite being unable to size a table")
	}

	for _, p := range []string{CallsPath(base), SitesPath(base), RegionsPath(base), VolumeManifestPath(base)} {
		if _, err := os.Stat(p); err == nil {
			t.Errorf("%s survived a failed Finish", filepath.Base(p))
		}
	}
}
