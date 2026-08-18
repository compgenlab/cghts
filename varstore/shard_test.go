package varstore

import (
	"fmt"
	"path/filepath"
	"testing"
)

// Splitting a table by coordinate.
//
// It is a layout change and nothing else, so the test that matters is that a
// split store and an unsplit one built from the same rows answer identically --
// every locus, every sample, every state. A faster store that disagrees with
// the slow one is not an optimisation.

// buildStore writes the same rows either split or whole.
func buildShardedStore(t *testing.T, shardSites int64) string {
	t.Helper()
	dir := filepath.Join(t.TempDir(), "store")
	var curChrom string
	// Open runs, held the way the converter holds them: one per sample, emitted
	// when coverage breaks OR when a shard boundary arrives.
	type open struct {
		start, last int32
		n           int32
	}
	runs := map[string]*open{}
	var w *Writer
	emit := func() error {
		for _, name := range []string{"S1", "S2", "S3"} {
			r := runs[name]
			if r == nil {
				continue
			}
			delete(runs, name)
			if err := w.WriteRegion(CalledSiteRun{
				SampleID: name, Chrom: curChrom, Start: r.start, End: r.last,
				NSites: r.n, MinDP: 30,
			}); err != nil {
				return err
			}
		}
		return nil
	}

	var err error
	w, err = NewWriter(dir, WriterOpts{
		Samples: []string{"S1", "S2", "S3"}, MinDP: 10, ShardSites: shardSites,
		// The contract: flush anything spanning sites into the shard it belongs
		// to, before the boundary closes it.
		BeforeRotate: emit,
	})
	if err != nil {
		t.Fatal(err)
	}

	// Two chromosomes, so a chromosome change forces a boundary regardless of
	// the site count -- a shard spanning chromosomes would make its First and
	// Last meaningless.
	for _, chrom := range []string{"chr1", "chr2"} {
		// A chromosome change also ends every run, exactly as it does in the
		// converter -- a run cannot span chromosomes any more than a shard can.
		if curChrom != "" && curChrom != chrom {
			if err := emit(); err != nil {
				t.Fatal(err)
			}
		}
		curChrom = chrom
		for i := 0; i < 25; i++ {
			pos := int32(100 + i*10)
			if err := w.WriteSite(Site{
				Chrom: chrom, Pos: pos, Ref: "G", Alt: "A", AN: 6, NCalled: 3,
			}); err != nil {
				t.Fatal(err)
			}
			// S1 and S2 are covered across the whole chromosome; S3 only the
			// first half, so it is reference early and NOT ASSAYED late.
			covered := []string{"S1", "S2"}
			if i < 13 {
				covered = append(covered, "S3")
			}
			for _, name := range covered {
				if r := runs[name]; r != nil {
					r.last, r.n = pos, r.n+1
				} else {
					runs[name] = &open{start: pos, last: pos, n: 1}
				}
			}

			// Every third site has a carrier, and every seventh a gated one.
			if i%3 == 0 {
				dp := int32(40)
				if i%7 == 0 {
					dp = 4 // below the gate: uncertain, not reference
				}
				if err := w.WriteCall(Call{
					SampleID: "S1", Chrom: chrom, Pos: pos, Ref: "G", Alt: "A",
					GT: "0/1", DP: dp, ADRef: dp / 2, ADAlt: dp - dp/2, GQ: Missing,
				}); err != nil {
					t.Fatal(err)
				}
			}
		}
	}
	if err := emit(); err != nil {
		t.Fatal(err)
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
	return dir
}

func allLoci() []Locus {
	var out []Locus
	for _, chrom := range []string{"chr1", "chr2"} {
		for i := 0; i < 25; i++ {
			out = append(out, Locus{Chrom: chrom, Pos: int32(100 + i*10), Ref: "G", Alt: "A"})
		}
	}
	// And one the store never interrogated.
	out = append(out, Locus{Chrom: "chr1", Pos: 999, Ref: "T", Alt: "C"})
	return out
}

func TestSplitStoreAnswersIdentically(t *testing.T) {
	whole, err := OpenParquet(buildShardedStore(t, 0))
	if err != nil {
		t.Fatal(err)
	}
	defer whole.Close()

	// 10 sites a shard over 50 sites on two chromosomes: several shards, and
	// boundaries that fall both on and off the chromosome change.
	split, err := OpenParquet(buildShardedStore(t, 10))
	if err != nil {
		t.Fatal(err)
	}
	defer split.Close()

	gate := Gate{MinDP: 10}
	for _, l := range allLoci() {
		a, err := whole.Classify(l, gate)
		if err != nil {
			t.Fatalf("%s whole: %v", l, err)
		}
		b, err := split.Classify(l, gate)
		if err != nil {
			t.Fatalf("%s split: %v", l, err)
		}
		if len(a) != len(b) {
			t.Fatalf("%s: whole %d states, split %d", l, len(a), len(b))
		}
		byName := map[string]State{}
		for _, st := range b {
			byName[st.SampleID] = st.State
		}
		for _, st := range a {
			if got := byName[st.SampleID]; got != st.State {
				t.Errorf("%s %s: whole %q, split %q", l, st.SampleID, st.State, got)
			}
		}
	}
}

// The batched read must agree too, since that is the path a real query takes.
func TestSplitStoreClassifyManyAgrees(t *testing.T) {
	whole, err := OpenParquet(buildShardedStore(t, 0))
	if err != nil {
		t.Fatal(err)
	}
	defer whole.Close()
	split, err := OpenParquet(buildShardedStore(t, 10))
	if err != nil {
		t.Fatal(err)
	}
	defer split.Close()

	loci, gate := allLoci(), Gate{MinDP: 10}
	a, err := whole.ClassifyMany(loci, gate)
	if err != nil {
		t.Fatal(err)
	}
	b, err := split.ClassifyMany(loci, gate)
	if err != nil {
		t.Fatal(err)
	}
	for _, l := range loci {
		x, y := a[l], b[l]
		if len(x) != len(y) {
			t.Fatalf("%s: whole %d, split %d", l, len(x), len(y))
		}
		byName := map[string]State{}
		for _, st := range y {
			byName[st.SampleID] = st.State
		}
		for _, st := range x {
			if got := byName[st.SampleID]; got != st.State {
				t.Errorf("%s %s: whole %q, split %q", l, st.SampleID, st.State, got)
			}
		}
	}
}

// A shard must never span chromosomes: its First and Last are compared without
// a chromosome test on every row, which is only sound if the shard has one.
func TestShardsNeverSpanChromosomes(t *testing.T) {
	s, err := OpenParquet(buildShardedStore(t, 1000)) // far more than the 25 per chromosome
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	m := s.VolumeManifest()
	shards := m.Tables[SitesTable].Shards
	if len(shards) < 2 {
		t.Fatalf("got %d shards, want at least 2 -- the chromosome change must cut one "+
			"even when the site count is nowhere near the limit", len(shards))
	}
	seen := map[string]bool{}
	for _, si := range shards {
		if seen[si.Chrom] {
			t.Errorf("chromosome %s appears in more than one shard, which is fine, "+
				"but check the ordering", si.Chrom)
		}
		seen[si.Chrom] = true
		if si.First > si.Last {
			t.Errorf("shard %s has First %d past Last %d", si.Name, si.First, si.Last)
		}
	}
}

// The index must describe what is actually there, or a reader prunes on a lie.
func TestShardIndexMatchesTheFiles(t *testing.T) {
	dir := buildShardedStore(t, 10)
	s, err := OpenParquet(dir)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	m := s.VolumeManifest()
	for _, table := range []string{CallsTable, SitesTable, RegionsTable} {
		info := m.Tables[table]
		if len(info.Shards) == 0 {
			t.Errorf("%s was not split", table)
			continue
		}
		var rows int64
		for i, si := range info.Shards {
			if si.Name != ShardFile(table, i) {
				t.Errorf("%s shard %d named %q, want %q", table, i, si.Name, ShardFile(table, i))
			}
			if si.Bytes <= 0 {
				t.Errorf("%s has no recorded size", si.Name)
			}
			rows += si.Rows
		}
		if rows != info.Rows {
			t.Errorf("%s shards total %d rows, table records %d", table, rows, info.Rows)
		}
	}
}

// Every run must lie wholly inside one shard. Without that a query would have
// to read every earlier shard in case a run started there and reached in --
// which is the whole reason the shard index can be trusted to prune.
func TestRunsDoNotCrossShards(t *testing.T) {
	s, err := OpenParquet(buildShardedStore(t, 10))
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	bounds := map[string]ShardInfo{}
	for _, si := range s.VolumeManifest().Tables[RegionsTable].Shards {
		bounds[si.Name] = si
	}
	if len(bounds) < 2 {
		t.Fatal("fixture produced fewer than two region shards; it proves nothing")
	}
	// Read every run and check it against the shard it was filed under. The
	// reader does not expose which shard a row came from, so this asserts the
	// weaker but sufficient property: no run spans a recorded boundary.
	var crossings int
	err = s.Regions(func(r CalledSiteRun) bool {
		for _, si := range bounds {
			if !SameChrom(si.Chrom, r.Chrom) {
				continue
			}
			startsIn := r.Start >= si.First && r.Start <= si.Last
			endsAfter := r.End > si.Last
			if startsIn && endsAfter {
				crossings++
				t.Errorf("run %s %s:%d-%d starts in shard %s (%d-%d) and reaches past it",
					r.SampleID, r.Chrom, r.Start, r.End, si.Name, si.First, si.Last)
			}
		}
		return crossings < 5
	})
	if err != nil {
		t.Fatal(err)
	}
}

func TestUnsplitStoreHasNoShardIndex(t *testing.T) {
	s, err := OpenParquet(buildShardedStore(t, 0))
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	for _, table := range []string{CallsTable, SitesTable, RegionsTable} {
		if n := len(s.VolumeManifest().Tables[table].Shards); n != 0 {
			t.Errorf("%s recorded %d shards in an unsplit store", table, n)
		}
	}
	fmt.Fprint(devNull{}, "")
}

type devNull struct{}

func (devNull) Write(p []byte) (int, error) { return len(p), nil }

// A caller that lets a run land outside the open shard is REFUSED, in both
// directions.
//
// This is the bug the differential comparison found, and it is worth a test
// because neither mistake announces itself. A run beginning EARLIER than the
// shard was never broken at the boundary. One beginning LATER was created after
// the boundary passed and filed in the shard that is closing -- which is what
// happens when a caller extends a run into the next shard's first site and only
// then writes that site, letting the writer discover the rotation too late.
//
// Either way the run lands in a shard whose range does not contain it, so no
// query for those positions ever sees it, and every locus it covered reads as
// never assayed. On a 200-sample fixture that was 1,796 wrong states out of
// 80,000 -- nine boundaries times two hundred people -- and nothing errored.
func TestRunOutsideTheOpenShardIsRefused(t *testing.T) {
	for _, c := range []struct {
		name string
		run  CalledSiteRun
	}{
		{"begins before the shard", CalledSiteRun{
			SampleID: "S1", Chrom: "chr1", Start: 100, End: 150, NSites: 5, MinDP: 30,
		}},
		{"begins after the shard", CalledSiteRun{
			SampleID: "S1", Chrom: "chr1", Start: 9000, End: 9100, NSites: 5, MinDP: 30,
		}},
	} {
		t.Run(c.name, func(t *testing.T) {
			dir := filepath.Join(t.TempDir(), "store")
			w, err := NewWriter(dir, WriterOpts{
				Samples: []string{"S1"}, MinDP: 10, ShardSites: 3,
			})
			if err != nil {
				t.Fatal(err)
			}
			// Sites 200..400 land in the open shard; the runs above sit outside it.
			for _, pos := range []int32{200, 300, 400, 500} {
				if err := w.WriteSite(Site{
					Chrom: "chr1", Pos: pos, Ref: "G", Alt: "A", AN: 2, NCalled: 1,
				}); err != nil {
					t.Fatal(err)
				}
			}
			if err := w.WriteRegion(c.run); err == nil {
				t.Fatal("a run outside the open shard was accepted; every locus it covers " +
					"would read as never assayed and nothing would say so")
			}
			_ = w.Close()
		})
	}
}

// WouldRotate is what lets a caller get the ordering right, so it has to be
// answerable before the site is written rather than after.
func TestWouldRotateAnticipatesTheBoundary(t *testing.T) {
	dir := filepath.Join(t.TempDir(), "store")
	w, err := NewWriter(dir, WriterOpts{
		Samples: []string{"S1"}, MinDP: 10, ShardSites: 2,
	})
	if err != nil {
		t.Fatal(err)
	}
	defer w.Close()

	if w.WouldRotate("chr1") {
		t.Error("an empty writer wants to rotate before writing anything")
	}
	for i, pos := range []int32{100, 200} {
		if err := w.WriteSite(Site{Chrom: "chr1", Pos: pos, Ref: "G", Alt: "A"}); err != nil {
			t.Fatal(err)
		}
		if want := i == 1; w.WouldRotate("chr1") != want {
			t.Errorf("after %d sites WouldRotate = %v, want %v", i+1, !want, want)
		}
	}
	// A chromosome change rotates whatever the count.
	if !w.WouldRotate("chr2") {
		t.Error("a chromosome change must rotate: a shard's First and Last are compared " +
			"without a chromosome test, which is only sound if it holds one")
	}
}
