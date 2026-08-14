package varstore

// Classifying a whole locus set in one pass.
//
// WHY THIS EXISTS. Classify answers one locus, and answering N of them meant
// calling it N times -- which re-scanned calls.parquet and RELOADED EVERY
// CALLABLE RUN, from scratch, once per locus. Measured on a 500-sample store:
// 298ms for one locus, 889ms for three, exactly linear. A gene's two hundred
// qualifying sites would have cost a minute, and the sites are not even the
// expensive part -- the runs are, and they are identical every time.
//
// So the work is inverted: one pruned pass over calls, one over regions, both
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
func (s *ParquetStore) ClassifyMany(loci []Locus, g Gate) (map[Locus][]SampleState, error) {
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

	// Which of the requested loci the source actually reported. Asked once for
	// the whole set: an off-catalog locus must not consult the run intervals,
	// because a run bracketing a position the source never mentioned does not
	// mean the position was called.
	want := make(map[Locus]bool, len(loci))
	for _, l := range loci {
		want[canonLocus(l)] = true
	}
	interrogated := make(map[Locus]bool, len(loci))
	q := Query{Loci: loci}
	if err := scanParquetPruned(s.sites, callsFilter(q), func(site Site) bool {
		if k := canonLocus(site.Locus()); want[k] {
			interrogated[k] = true
		}
		return true
	}); err != nil {
		return nil, err
	}

	// ONE pass over the calls for every locus in the set.
	calls := make(map[Locus]map[string]Call, len(loci))
	if err := scanParquetPruned(s.calls, callsFilter(q), func(c Call) bool {
		k := canonLocus(c.Locus())
		if !want[k] {
			return true
		}
		if calls[k] == nil {
			calls[k] = map[string]Call{}
		}
		calls[k][c.SampleID] = c
		return true
	}); err != nil {
		return nil, err
	}

	// ONE pass over the runs, and the reason this function is worth having. A
	// run is kept if it covers ANY requested locus, and each is recorded
	// against the loci it actually brackets -- so the file is read once rather
	// than once per locus.
	//
	// Positions are collected per chromosome and sorted, so a run tests against
	// a slice rather than against the whole set.
	byChrom := map[string][]Locus{}
	for l := range want {
		key := CanonKey(l.Chrom)
		byChrom[key] = append(byChrom[key], l)
	}
	called := make(map[Locus]map[string]bool, len(loci))
	if err := scanParquetPruned(s.regions, runScanFilter(q), func(r CalledSiteRun) bool {
		for _, l := range byChrom[CanonKey(r.Chrom)] {
			if l.Pos < r.Start || l.Pos > r.End {
				continue
			}
			if called[l] == nil {
				called[l] = map[string]bool{}
			}
			called[l][r.SampleID] = true
		}
		return true
	}); err != nil {
		return nil, err
	}

	// The same decision Classify makes, per sample per locus.
	for _, raw := range loci {
		l := canonLocus(raw)
		states := make([]SampleState, 0, len(samples))
		if !interrogated[l] && s.spans != SpansBlocks {
			for _, name := range samples {
				states = append(states, SampleState{SampleID: name, State: StateNotAssayed})
			}
			out[raw] = states
			continue
		}
		at, cov := calls[l], called[l]
		for _, name := range samples {
			st := SampleState{SampleID: name}
			if c, ok := at[name]; ok {
				cc := c
				st.Call = &cc
				if g.Admits(c) {
					st.State = StateCarrier
				} else {
					// An ALT call we do not believe. NOT a reference, and not
					// an absence: we looked, and saw something we cannot vouch
					// for.
					st.State = StateUncertain
				}
			} else if cov[name] {
				st.State = StateNonCarrier
			} else {
				st.State = StateNotAssayed
			}
			states = append(states, st)
		}
		out[raw] = states
	}
	return out, nil
}
