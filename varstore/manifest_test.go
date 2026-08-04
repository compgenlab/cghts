package varstore

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// buildCensusStore writes a small two-chromosome store and finishes it.
func buildCensusStore(t *testing.T, base string, opts WriterOpts) {
	t.Helper()
	if opts.Samples == nil {
		opts.Samples = []string{"S1", "S2"}
	}
	if opts.RowGroupSize == 0 {
		opts.RowGroupSize = 4
	}
	w, err := NewWriter(base, opts)
	if err != nil {
		t.Fatal(err)
	}
	for _, chrom := range []string{"chr1", "chr2"} {
		for pos := int32(100); pos < 105; pos++ {
			if err := w.WriteSite(Site{
				Chrom: chrom, Pos: pos, Ref: "A", Alt: "T", AC: 1, AN: 4, NCarriers: 1, NCalled: 2,
			}); err != nil {
				t.Fatal(err)
			}
			if err := w.WriteCall(Call{
				SampleID: "S1", Chrom: chrom, Pos: pos, Ref: "A", Alt: "T", GT: "0/1", DP: 30,
			}); err != nil {
				t.Fatal(err)
			}
		}
		if err := w.WriteRegion(CalledSiteRun{
			SampleID: "S1", Chrom: chrom, Start: 100, End: 104, NSites: 5,
		}); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
}

func TestManifestRecordsWhatWasWritten(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	buildCensusStore(t, base, WriterOpts{
		MinDP:   10,
		Program: "cgkit test",
		Command: "cgkit vcf-toparquet chr1.vcf chr2.vcf --out cohort",
		Sources: []string{"chr1.vcf", "chr2.vcf"},
		Contigs: []string{"##contig=<ID=chr1>", "##contig=<ID=chr2>"},
	})

	m, err := ReadManifest(base)
	if err != nil {
		t.Fatal(err)
	}
	if !m.Complete {
		t.Error("Complete is false on a finished store")
	}
	if m.FormatVersion != ManifestVersion {
		t.Errorf("FormatVersion = %d, want %d", m.FormatVersion, ManifestVersion)
	}
	if m.Created.IsZero() {
		t.Error("Created is zero")
	}
	if got := m.Counts; got.Calls != 10 || got.Sites != 10 || got.Regions != 2 || got.Samples != 2 {
		t.Errorf("counts = %+v, want 10 calls / 10 sites / 2 regions / 2 samples", got)
	}
	if got, want := len(m.Sources), 2; got != want {
		t.Errorf("Sources has %d entries, want %d", got, want)
	}
	if m.Params.MinDP != 10 || m.Params.SpanSemantics != SpansSites {
		t.Errorf("params = %+v, want min_dp 10 and sites semantics", m.Params)
	}

	// The census is the field that can contradict the intent stamped at
	// construction, so it has to describe the rows rather than the inputs.
	if len(m.Chromosomes) != 2 {
		t.Fatalf("census has %d chromosomes, want 2: %+v", len(m.Chromosomes), m.Chromosomes)
	}
	if m.Chromosomes[0].Name != "chr1" || m.Chromosomes[1].Name != "chr2" {
		t.Errorf("census is not in written order: %+v", m.Chromosomes)
	}
	for _, c := range m.Chromosomes {
		if c.Sites != 5 || c.Calls != 5 {
			t.Errorf("%s: %d sites / %d calls, want 5 each", c.Name, c.Sites, c.Calls)
		}
		if c.FirstPos != 100 || c.LastPos != 104 {
			t.Errorf("%s: positions %d-%d, want 100-104", c.Name, c.FirstPos, c.LastPos)
		}
	}

	for _, name := range []string{CallsMember, SitesMember, RegionsMember} {
		info, ok := m.Members[name]
		if !ok {
			t.Errorf("no member entry for %s", name)
			continue
		}
		st, err := os.Stat(MemberPath(base, name))
		if err != nil {
			t.Fatal(err)
		}
		if info.Bytes != st.Size() {
			t.Errorf("%s: manifest says %d bytes, file is %d", name, info.Bytes, st.Size())
		}
	}
}

// A store is readable only once a conversion says it finished. Everything else
// about the members can look perfect while the run covered a fraction of its
// input, and that store answers "not assayed" for the rest -- which is exactly
// what a correct store says about a position the source never reported.
func TestOpenRefusesAStoreWithNoManifest(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	buildCensusStore(t, base, WriterOpts{MinDP: 10})

	if err := os.Remove(ManifestPath(base)); err != nil {
		t.Fatal(err)
	}
	_, err := OpenParquet(base)
	if err == nil {
		t.Fatal("a store with no manifest opened")
	}
	if !strings.Contains(err.Error(), ManifestFile) {
		t.Errorf("error does not name the manifest: %v", err)
	}
}

func TestOpenRefusesAnIncompleteManifest(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	buildCensusStore(t, base, WriterOpts{MinDP: 10})

	m, err := ReadManifest(base)
	if err != nil {
		t.Fatal(err)
	}
	m.Complete = false
	if err := WriteManifest(base, *m); err != nil {
		t.Fatal(err)
	}
	if _, err := OpenParquet(base); err == nil {
		t.Fatal("a store marked incomplete opened")
	}
}

func TestOpenRefusesAFutureManifestVersion(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	buildCensusStore(t, base, WriterOpts{MinDP: 10})

	m, err := ReadManifest(base)
	if err != nil {
		t.Fatal(err)
	}
	m.FormatVersion = ManifestVersion + 1
	if err := WriteManifest(base, *m); err != nil {
		t.Fatal(err)
	}
	if _, err := OpenParquet(base); err == nil {
		t.Fatal("a store from a newer format version opened")
	}
}

// The row-count check is what catches a member that is well-formed but belongs
// to a different conversion. Sites and regions carry no metadata of their own,
// so without the manifest nothing links a member to the store around it.
func TestOpenRefusesAMemberFromAnotherStore(t *testing.T) {
	dir := t.TempDir()
	mine := filepath.Join(dir, "mine")
	theirs := filepath.Join(dir, "theirs")
	buildCensusStore(t, mine, WriterOpts{MinDP: 10})

	// A second store with a different number of sites.
	w, err := NewWriter(theirs, WriterOpts{Samples: []string{"S1", "S2"}, MinDP: 10, RowGroupSize: 4})
	if err != nil {
		t.Fatal(err)
	}
	if err := w.WriteSite(Site{Chrom: "chr1", Pos: 100, Ref: "A", Alt: "T", AC: 1, AN: 4}); err != nil {
		t.Fatal(err)
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}

	swapped, err := os.ReadFile(SitesPath(theirs))
	if err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(SitesPath(mine), swapped, 0o644); err != nil {
		t.Fatal(err)
	}

	_, err = OpenParquet(mine)
	if err == nil {
		t.Fatal("a store opened with a sites member from a different conversion")
	}
	if !strings.Contains(err.Error(), "does not belong") {
		t.Errorf("error does not explain the mismatch: %v", err)
	}
}

// A conversion under way must never carry the previous run's completion marker,
// or a --force retry that dies leaves a manifest vouching for half-written
// members.
func TestNewWriterClearsAStaleManifest(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	buildCensusStore(t, base, WriterOpts{MinDP: 10})
	if _, err := os.Stat(ManifestPath(base)); err != nil {
		t.Fatal(err)
	}

	w, err := NewWriter(base, WriterOpts{Samples: []string{"S1"}, MinDP: 10})
	if err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(ManifestPath(base)); !os.IsNotExist(err) {
		t.Error("the previous manifest survived into a new conversion")
	}
	if _, err := OpenParquet(base); err == nil {
		t.Error("a store under construction was readable")
	}
	if err := w.Discard(); err != nil {
		t.Fatal(err)
	}
}
