package varstore

import (
	"os"
	"runtime"
	"strconv"
	"sync"
)

// Classifying a whole locus set in one pass.
//
// WHY THIS EXISTS. Classify answers one locus, and answering N of them meant
// calling it N times -- which re-scanned calls.parquet and RELOADED EVERY
// CALLABLE RUN, from scratch, once per locus. Measured on a 500-sample store:
// 298ms for one locus, 889ms for three, exactly linear. A gene's two hundred
// qualifying sites would have cost a minute, and the sites are not even the
// expensive part -- the runs are, and they are identical every time.
//
// So the work is inverted: one pruned pass over calls, one over catalog_runs, both
// bounded by the whole locus set at once, and then the same per-sample decision
// Classify makes. Cost becomes a function of the DATA rather than of the number
// of loci asked about.
//
// WHY NOT BUILD IT ON Calls, which already takes a locus set. Calls drops a
// gated-out ALT call from its output while still recording that the sample
// carries something there, so the sample falls through every branch and is
// simply absent -- which a caller reads as "never assayed". Classify calls the
// same case StateUncertain. Those are different claims about a person, and the
// distinction is the one this package exists to preserve, so the batched
// version reproduces Classify's decision rather than Calls'.

// ClassifyMany returns every sample's state at every requested locus.
//
// The result is keyed by locus; a locus the store never interrogated maps to
// every sample as StateNotAssayed, exactly as Classify reports it.
func (s *ParquetVolume) ClassifyMany(loci []Locus, g Gate) (map[Locus][]SampleState, error) {
	if err := s.classifiable(); err != nil {
		return nil, err
	}
	out := make(map[Locus][]SampleState, len(loci))
	if len(loci) == 0 {
		return out, nil
	}
	samples, err := s.Samples()
	if err != nil {
		return nil, err
	}

	// THREE MEMBERS, THREE GOROUTINES, and each internally parallel across the
	// shards it must visit.
	//
	// The scans are independent -- sites says what was interrogated, calls says
	// who carried something, catalog_runs says who was covered -- so the query costs
	// the SLOWEST of them rather than their sum. Measured sequentially at
	// 175ms for calls and 123ms for catalog_runs on a 500-sample store, which is
	// 298ms of waiting for 175ms of work.
	//
	// Each table is a distinct set of files, and within a table each shard is
	// a distinct file, so no two goroutines ever touch one reader. That is the
	// rule the whole scheme rests on: an io.SectionReader is not safe for
	// concurrent use, and the unit of parallelism is therefore the file.
	workers := runtime.GOMAXPROCS(0)
	// An escape hatch for measuring, and for a deployment that would rather
	// spend its cores elsewhere. Setting it to 1 restores the sequential scan
	// exactly.
	if n := os.Getenv("VARSTORE_SCAN_WORKERS"); n != "" {
		if v, err := strconv.Atoi(n); err == nil && v > 0 {
			workers = v
		}
	}

	// Which of the requested loci the source actually reported. Asked once for
	// the whole set: an off-catalog locus must not consult the run intervals,
	// because a run bracketing a position the source never mentioned does not
	// mean the position was called.
	want := make(map[Locus]bool, len(loci))
	for _, l := range loci {
		want[canonLocus(l)] = true
	}
	interrogated := make(map[Locus]bool, len(loci))
	calls := make(map[Locus]map[string]Call, len(loci))
	called := make(map[Locus]map[string]bool, len(loci))
	q := Query{Loci: loci}

	// Positions per chromosome, so a run tests against a slice rather than the
	// whole set. Built before the scans start, since all three read it.
	byChrom := map[string][]Locus{}
	for l := range want {
		byChrom[CanonKey(l.Chrom)] = append(byChrom[CanonKey(l.Chrom)], l)
	}

	var wg sync.WaitGroup
	errs := make([]error, 4)
	blocked := make(map[Locus]map[string]bool, len(loci))

	// 1. Which of the requested loci the source actually reported. An
	//    off-catalog locus must not consult the run intervals, because a run
	//    bracketing a position the source never mentioned does not mean the
	//    position was examined.
	wg.Add(1)
	go func() {
		defer wg.Done()
		errs[0] = scanShardsParallel(s.sites, callsFilter(q), workers,
			func() map[Locus]bool { return map[Locus]bool{} },
			func(acc map[Locus]bool, site Site) bool {
				if k := canonLocus(site.Locus()); want[k] {
					acc[k] = true
				}
				return true
			},
			func(acc map[Locus]bool) {
				for k := range acc {
					interrogated[k] = true
				}
			})
	}()

	// 2. The ALT calls, for every locus in the set at once.
	wg.Add(1)
	go func() {
		defer wg.Done()
		errs[1] = scanShardsParallel(s.calls, callsFilter(q), workers,
			func() map[Locus]map[string]Call { return map[Locus]map[string]Call{} },
			func(acc map[Locus]map[string]Call, c Call) bool {
				k := canonLocus(c.Locus())
				if !want[k] {
					return true
				}
				if acc[k] == nil {
					acc[k] = map[string]Call{}
				}
				acc[k][c.SampleID] = c
				return true
			},
			func(acc map[Locus]map[string]Call) {
				for k, at := range acc {
					if calls[k] == nil {
						calls[k] = at
						continue
					}
					for name, c := range at {
						calls[k][name] = c
					}
				}
			})
	}()

	// 3. The catalog runs, and the expensive one. A run is kept if it covers any
	//    requested locus, and recorded against the loci it brackets.
	wg.Add(1)
	go func() {
		defer wg.Done()
		errs[2] = scanShardsParallel(s.catalogRuns, runScanFilter(q), workers,
			func() map[Locus]map[string]bool { return map[Locus]map[string]bool{} },
			func(acc map[Locus]map[string]bool, r CalledSiteRun) bool {
				// A run the gate refuses supports no reference call, so it is
				// not recorded at all: the sample stays not-assayed here.
				if !s.runAdmits(r, g) {
					return true
				}
				for _, l := range byChrom[CanonKey(r.Chrom)] {
					if l.Pos < r.Start || l.Pos > r.End {
						continue
					}
					if acc[l] == nil {
						acc[l] = map[string]bool{}
					}
					acc[l][r.SampleID] = true
				}
				return true
			},
			func(acc map[Locus]map[string]bool) {
				for k, cov := range acc {
					if called[k] == nil {
						called[k] = cov
						continue
					}
					for name := range cov {
						called[k][name] = true
					}
				}
			})
	}()

	// 4. The reference blocks, on a store that has them. Read for every locus
	//    while the catalog scan is still deciding which loci are off it, so the
	//    batch costs the slowest scan rather than an extra round; the decision
	//    below consults a locus's blocks only once the catalog has said it is
	//    not a site. Reading them for a catalog site is wasted work, never a
	//    wrong answer.
	if s.hasRefBlocks {
		wg.Add(1)
		go func() {
			defer wg.Done()
			errs[3] = scanShardsParallel(s.refBlocks, refBlockFilter(loci), workers,
				func() map[Locus]map[string]bool { return map[Locus]map[string]bool{} },
				func(acc map[Locus]map[string]bool, b RefBlock) bool {
					if !blockAdmits(b, g) {
						return true
					}
					for _, l := range byChrom[CanonKey(b.Chrom)] {
						if l.Pos < b.Start || l.Pos > b.End {
							continue
						}
						if acc[l] == nil {
							acc[l] = map[string]bool{}
						}
						acc[l][b.SampleID] = true
					}
					return true
				},
				func(acc map[Locus]map[string]bool) {
					for k, cov := range acc {
						if blocked[k] == nil {
							blocked[k] = cov
							continue
						}
						for name := range cov {
							blocked[k][name] = true
						}
					}
				})
		}()
	}

	wg.Wait()
	for _, err := range errs {
		if err != nil {
			return nil, err
		}
	}

	// The same decision Classify makes, per sample per locus -- including which
	// reference source the catalog selects, so a block never answers at a site.
	for _, raw := range loci {
		l := canonLocus(raw)
		states := make([]SampleState, 0, len(samples))
		if !interrogated[l] && !s.hasRefBlocks {
			for _, name := range samples {
				states = append(states, SampleState{SampleID: name, State: StateNotAssayed})
			}
			out[raw] = states
			continue
		}
		ref, source := called[l], RefFromCatalog
		at := calls[l]
		if !interrogated[l] {
			// Calls only exist at catalog sites; a stray match here would be a
			// store disagreeing with itself, and the catalog is the authority.
			ref, source, at = blocked[l], RefFromBlock, nil
		}
		for _, name := range samples {
			var c *Call
			if x, ok := at[name]; ok {
				c = &x
			}
			states = append(states, resolveState(name, c, ref[name], source, g))
		}
		out[raw] = states
	}
	return out, nil
}
