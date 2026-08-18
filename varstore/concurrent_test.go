package varstore

import (
	"context"
	"path/filepath"
	"reflect"
	"sync"
	"testing"
)

// One open store serves concurrent queries.
//
// THIS IS WHAT LETS A SERVER HOLD ONE HANDLE instead of a pool of them. The
// expensive part of opening is the parsed parquet footers, and a pool would
// carry one copy per handle -- duplicating the very thing it exists to keep,
// and paying a fresh parse for each handle a burst forces open. That is only
// avoidable if a single handle is safe to share.
//
// The read path is built for it: shards already scan in parallel within one
// call, the lazily parsed footer is guarded by a sync.Once, and reads go
// through ReadAt, which io.SectionReader documents as safe for concurrent use.
// This asserts the property end to end rather than by inspection, because the
// failure mode is a data race -- silent, load-dependent, and invisible to every
// other test in this package.
//
// Run under -race or it proves very little.
func TestOneStoreServesConcurrentQueries(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	buildCensusStore(t, base, WriterOpts{MinDP: 10})

	s, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	gate := Gate{MinDP: 10}
	loci := []Locus{
		{Chrom: "chr1", Pos: 100, Ref: "A", Alt: "T"},
		{Chrom: "chr1", Pos: 102, Ref: "A", Alt: "T"},
		{Chrom: "chr2", Pos: 101, Ref: "A", Alt: "T"},
	}

	// The expected answers, computed serially first, so the concurrent run is
	// checked against a value rather than merely for absence of a crash.
	want := map[Locus][]SampleState{}
	for _, l := range loci {
		st, err := s.Classify(l, gate)
		if err != nil {
			t.Fatal(err)
		}
		want[l] = st
	}

	const goroutines = 16
	const rounds = 8
	var wg sync.WaitGroup
	for g := 0; g < goroutines; g++ {
		wg.Add(1)
		go func(g int) {
			defer wg.Done()
			for r := 0; r < rounds; r++ {
				// Deliberately mixed: the three entry points take different
				// paths through the same tables, and a race between two
				// DIFFERENT readers is the one inspection is worst at finding.
				switch (g + r) % 3 {
				case 0:
					l := loci[r%len(loci)]
					got, err := s.Classify(l, gate)
					if err != nil {
						t.Error(err)
						return
					}
					if !reflect.DeepEqual(got, want[l]) {
						t.Errorf("Classify(%v) disagrees with the serial answer under concurrency", l)
						return
					}
				case 1:
					got, err := s.ClassifyMany(loci, gate)
					if err != nil {
						t.Error(err)
						return
					}
					for _, l := range loci {
						if !reflect.DeepEqual(got[l], want[l]) {
							t.Errorf("ClassifyMany(%v) disagrees with the serial answer under concurrency", l)
							return
						}
					}
				case 2:
					n := 0
					if err := s.Sites(func(Site) bool { n++; return true }); err != nil {
						t.Error(err)
						return
					}
					if n != 10 {
						t.Errorf("a concurrent catalog walk saw %d sites, want 10", n)
						return
					}
				}
			}
		}(g)
	}
	wg.Wait()
}

// And the same for a multi-volume archive, whose lazy per-volume open is a
// second piece of mutable state on the read path.
func TestOneArchiveServesConcurrentQueries(t *testing.T) {
	base := buildSet(t, []string{"chr1", "chr2"}, []string{"S1", "S2", "S3"})

	s, err := OpenStore(context.Background(), base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	gate := Gate{MinDP: 10}
	loci := []Locus{
		{Chrom: "chr1", Pos: 100, Ref: "G", Alt: "A"},
		{Chrom: "chr1", Pos: 130, Ref: "G", Alt: "A"},
		{Chrom: "chr2", Pos: 100, Ref: "G", Alt: "A"},
	}
	want, err := s.ClassifyMany(loci, gate)
	if err != nil {
		t.Fatal(err)
	}

	var wg sync.WaitGroup
	for g := 0; g < 16; g++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for r := 0; r < 8; r++ {
				got, err := s.ClassifyMany(loci, gate)
				if err != nil {
					t.Error(err)
					return
				}
				for _, l := range loci {
					if !reflect.DeepEqual(got[l], want[l]) {
						t.Errorf("ClassifyMany(%v) disagrees with the serial answer under concurrency", l)
						return
					}
				}
			}
		}()
	}
	wg.Wait()
}

// And with the shards that make a query parallel in the first place.
//
// scanShardsParallel fans a single query across shards, so a shared handle puts
// TWO layers of concurrency on the same tables: several queries at once, each
// fanning out internally. That is the combination a pool would have hidden and
// the one most likely to race, so it is worth its own case rather than trusting
// the unsharded result to cover it.
func TestOneShardedStoreServesConcurrentQueries(t *testing.T) {
	base := buildShardedStore(t, 3)

	s, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	if n := len(s.VolumeManifest().Tables[CallsTable].Shards); n < 2 {
		t.Fatalf("the fixture produced %d shards, so it does not exercise the parallel scan", n)
	}

	gate := Gate{MinDP: 10}
	var loci []Locus
	for _, chrom := range []string{"chr1", "chr2"} {
		for i := 0; i < 25; i++ {
			loci = append(loci, Locus{Chrom: chrom, Pos: int32(100 + i*10), Ref: "G", Alt: "A"})
		}
	}
	want, err := s.ClassifyMany(loci, gate)
	if err != nil {
		t.Fatal(err)
	}

	var wg sync.WaitGroup
	for g := 0; g < 16; g++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for r := 0; r < 8; r++ {
				got, err := s.ClassifyMany(loci, gate)
				if err != nil {
					t.Error(err)
					return
				}
				for _, l := range loci {
					if !reflect.DeepEqual(got[l], want[l]) {
						t.Errorf("ClassifyMany(%v) disagrees with the serial answer under concurrency", l)
						return
					}
				}
			}
		}()
	}
	wg.Wait()
}
