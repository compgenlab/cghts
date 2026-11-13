package varstore

import (
	"path/filepath"
	"reflect"
	"testing"
)

// The fixture for every test here.
//
//	catalog  chr1:1000 G>A and chr1:3000 G>A; S1 carries both
//	S1       carrier, no reference block
//	S2       catalog run over both sites at 30x, and a block 900-3100 at 30x
//	S3       NO catalog run -- the caller declined it -- but a block 900-3100 at 40x
//	S4       a block 1500-2500 at 5x
//	S5       a block 1900-2100 that vouches for no depth at all (min_dp 0)
//
// chr1:2000 is off the catalog and inside every block but S1's; chr1:5000 is
// off the catalog and outside them all.
var (
	rbSite    = Locus{Chrom: "chr1", Pos: 1000, Ref: "G", Alt: "A"}
	rbOff     = Locus{Chrom: "chr1", Pos: 2000, Ref: "C", Alt: "T"}
	rbOutside = Locus{Chrom: "chr1", Pos: 5000, Ref: "C", Alt: "T"}
	rbSamples = []string{"S1", "S2", "S3", "S4", "S5"}
)

func buildRefBlockStore(t *testing.T, withBlocks bool, shardSites int64) string {
	t.Helper()
	base := filepath.Join(t.TempDir(), "cohort")
	w, err := NewWriter(base, WriterOpts{
		Samples:      rbSamples,
		MinDP:        10,
		RefBlocks:    withBlocks,
		RowGroupSize: 4,
		ShardSites:   shardSites,
	})
	if err != nil {
		t.Fatal(err)
	}
	for _, pos := range []int32{1000, 3000} {
		if err := w.WriteSite(Site{
			Chrom: "chr1", Pos: pos, Ref: "G", Alt: "A", AC: 1, AN: 10, NCarriers: 1, NCalled: 2,
		}); err != nil {
			t.Fatal(err)
		}
		if err := w.WriteCall(Call{
			SampleID: "S1", Chrom: "chr1", Pos: pos, Ref: "G", Alt: "A",
			GT: "0/1", DP: 60, ADRef: 30, ADAlt: 30, GQ: 99,
		}); err != nil {
			t.Fatal(err)
		}
		// One run per site, so a sharded store never has a run straddling a
		// shard boundary. S3 has none: the caller declined to genotype it.
		for _, r := range []CalledSiteRun{
			{SampleID: "S1", Chrom: "chr1", Start: pos, End: pos, NSites: 1, MinDP: 60},
			{SampleID: "S2", Chrom: "chr1", Start: pos, End: pos, NSites: 1, MinDP: 30},
		} {
			if err := w.WriteRegion(r); err != nil {
				t.Fatal(err)
			}
		}
	}
	if withBlocks {
		for _, b := range []RefBlock{
			{SampleID: "S2", Chrom: "chr1", Start: 900, End: 3100, MinDP: 30, GQ: Missing},
			{SampleID: "S3", Chrom: "chr1", Start: 900, End: 3100, MinDP: 40, GQ: Missing},
			{SampleID: "S4", Chrom: "chr1", Start: 1500, End: 2500, MinDP: 5, GQ: Missing},
			{SampleID: "S5", Chrom: "chr1", Start: 1900, End: 2100, MinDP: 0, GQ: Missing},
		} {
			if err := w.WriteRefBlock(b); err != nil {
				t.Fatal(err)
			}
		}
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
	return base
}

type rbWant struct {
	state  State
	source RefSource
}

// classifyBoth runs Classify and ClassifyMany at one locus, fails if they
// disagree, and returns the states by sample.
func classifyBoth(t *testing.T, s *ParquetVolume, l Locus, g Gate) map[string]rbWant {
	t.Helper()
	one, err := s.Classify(l, g)
	if err != nil {
		t.Fatal(err)
	}
	many, err := s.ClassifyMany([]Locus{l, rbSite}, g)
	if err != nil {
		t.Fatal(err)
	}
	flat := func(states []SampleState) map[string]rbWant {
		out := map[string]rbWant{}
		for _, st := range states {
			out[st.SampleID] = rbWant{st.State, st.RefSource}
		}
		return out
	}
	a, b := flat(one), flat(many[l])
	if !reflect.DeepEqual(a, b) {
		t.Errorf("%s gate %d: Classify %v\n   ClassifyMany %v", l, g.MinDP, a, b)
	}
	return a
}

func TestRefBlocksAnswerOffTheCatalog(t *testing.T) {
	for _, shardSites := range []int64{0, 1} {
		s, err := OpenParquet(buildRefBlockStore(t, true, shardSites))
		if err != nil {
			t.Fatal(err)
		}

		nc := rbWant{StateNonCarrier, RefFromBlock}
		na := rbWant{StateNotAssayed, ""}
		for _, c := range []struct {
			gate int32
			want map[string]rbWant
		}{
			{10, map[string]rbWant{"S1": na, "S2": nc, "S3": nc, "S4": na, "S5": na}},
			// No gate admits the 5x block, and the 0x one too -- exactly as the
			// VCF backend reads a gVCF block under no gate.
			{0, map[string]rbWant{"S1": na, "S2": nc, "S3": nc, "S4": nc, "S5": nc}},
			{35, map[string]rbWant{"S1": na, "S2": na, "S3": nc, "S4": na, "S5": na}},
		} {
			got := classifyBoth(t, s, rbOff, Gate{MinDP: c.gate})
			if !reflect.DeepEqual(got, c.want) {
				t.Errorf("shards %d, off-catalog gate %d:\n got %v\nwant %v", shardSites, c.gate, got, c.want)
			}
		}

		for name, st := range classifyBoth(t, s, rbOutside, Gate{MinDP: 10}) {
			if st != na {
				t.Errorf("shards %d: %s outside every block = %v, want not assayed", shardSites, name, st)
			}
		}
		s.Close()
	}
}

// THE LOAD-BEARING CASE. S3 is in the catalog with no run: in a joint-called
// source that means the caller declined to genotype it there. A reference block
// covering the position says reads were present, which is not the same claim,
// and must not overrule the caller -- those disagreements cluster where calling
// is hard, so letting the block win fabricates reference calls in exactly the
// worst places.
func TestARefBlockDoesNotOverruleADeclinedGenotype(t *testing.T) {
	s, err := OpenParquet(buildRefBlockStore(t, true, 0))
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	got := classifyBoth(t, s, rbSite, Gate{MinDP: 10})
	if got["S3"].state != StateNotAssayed {
		t.Errorf("S3 at a catalog site with no run but a covering block = %v; "+
			"a block answered where the caller had declined", got["S3"])
	}
	if want := (rbWant{StateNonCarrier, RefFromCatalog}); got["S2"] != want {
		t.Errorf("S2 at a catalog site with a run = %v, want %v", got["S2"], want)
	}
	if got["S1"].state != StateCarrier || got["S1"].source != "" {
		t.Errorf("S1 = %v, want a carrier with no reference source", got["S1"])
	}
}

// A store without ref_blocks answers nothing off the catalog -- unchanged from
// before the table existed.
func TestWithoutRefBlocksOffCatalogIsNotAssayed(t *testing.T) {
	s, err := OpenParquet(buildRefBlockStore(t, false, 0))
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	if s.HasRefBlocks() || s.Provenance().RefBlocks {
		t.Fatal("a store written without blocks reports having them")
	}
	for name, st := range classifyBoth(t, s, rbOff, Gate{MinDP: 10}) {
		if st.state != StateNotAssayed {
			t.Errorf("%s off the catalog of a store with no blocks = %v", name, st)
		}
	}
}

// Calls(IncludeRef) answers a NAMED off-catalog locus from the blocks, and
// never enumerates one out of a span.
func TestIncludeRefAnswersNamedOffCatalogLociOnly(t *testing.T) {
	s, err := OpenParquet(buildRefBlockStore(t, true, 0))
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	calls, err := CollectCalls(s, Query{
		Loci: []Locus{rbOff, rbSite}, IncludeRef: true, Gate: Gate{MinDP: 10},
	})
	if err != nil {
		t.Fatal(err)
	}
	var got []string
	for _, c := range calls {
		got = append(got, c.Locus().String()+" "+c.SampleID+" "+c.GT)
		if c.Pos == rbOff.Pos {
			if c.GT != HomRefGT {
				t.Errorf("off-catalog row %+v is not a reference call", c)
			}
			want := map[string]int32{"S2": 30, "S3": 40}[c.SampleID]
			if c.MinDP != want {
				t.Errorf("off-catalog 0/0 for %s carries min_dp %d, want the block's %d",
					c.SampleID, c.MinDP, want)
			}
		}
	}
	// Store order: the catalog site at 1000 before the named position at 2000,
	// whatever order the query named them in.
	want := []string{
		rbSite.String() + " S1 0/1",
		rbSite.String() + " S2 0/0",
		rbOff.String() + " S2 0/0",
		rbOff.String() + " S3 0/0",
	}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("rows:\n got %v\nwant %v", got, want)
	}

	// A span over the same position names no locus, so nothing off the catalog
	// may come back -- a block would otherwise yield a row per base per sample.
	spanned, err := CollectCalls(s, Query{
		Spans: []Span{{Chrom: "chr1", Start: 1500, End: 2500}}, IncludeRef: true,
	})
	if err != nil {
		t.Fatal(err)
	}
	if len(spanned) != 0 {
		t.Errorf("a span with no catalog site in it returned %d rows: %+v", len(spanned), spanned)
	}
}

// The VCF backend labels its reference observations the same way, so a
// consumer reconciling the two never sees provenance appear or vanish with the
// backend.
func TestVcfBackendReportsRefSource(t *testing.T) {
	path := filepath.Join(t.TempDir(), "in.vcf")
	if err := writeFile(path, `##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	1000	.	G	A	50	PASS	.	GT:DP	0/1:30	0/0:30
`); err != nil {
		t.Fatal(err)
	}
	s, err := OpenVcf(path)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	states, err := s.Classify(rbSite, Gate{MinDP: 10})
	if err != nil {
		t.Fatal(err)
	}
	for _, st := range states {
		want := RefSource("")
		if st.State == StateNonCarrier {
			want = RefFromCatalog
		}
		if st.RefSource != want {
			t.Errorf("%s: %s with RefSource %q, want %q", st.SampleID, st.State, st.RefSource, want)
		}
	}
}

func TestARefBlockForASampleOffTheRosterIsRefused(t *testing.T) {
	w, err := NewWriter(filepath.Join(t.TempDir(), "s"), WriterOpts{
		Samples: []string{"S1"}, RefBlocks: true,
	})
	if err != nil {
		t.Fatal(err)
	}
	defer w.Discard()
	if err := w.WriteRefBlock(RefBlock{SampleID: "NOPE", Chrom: "chr1", Start: 1, End: 10, MinDP: 20}); err == nil {
		t.Error("a block for a sample the store does not hold was accepted")
	}
}
