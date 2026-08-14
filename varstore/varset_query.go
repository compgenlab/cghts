package varstore

import (
	"fmt"
	"iter"
	"runtime"
	"sort"
	"sync"
)

// Querying a set.
//
// A locus is answered by exactly one member, so dispatch is a lookup rather
// than a merge -- which is the whole reason a set is not a set of parts. What
// is left is concurrency: a question spanning several chromosomes is several
// independent questions, and they are asked at once.

// SiteKnown reports whether the set interrogated a locus.
func (s *VarSet) SiteKnown(l Locus) (bool, error) {
	st, ok, err := s.member(l.Chrom)
	if err != nil || !ok {
		return false, err
	}
	return st.SiteKnown(l)
}

// Classify resolves one locus, through the member that holds its chromosome.
func (s *VarSet) Classify(l Locus, g Gate) ([]SampleState, error) {
	st, ok, err := s.member(l.Chrom)
	if err != nil {
		return nil, err
	}
	if !ok {
		// No member covers this chromosome, so nobody looked. Every sample is
		// NotAssayed -- which is an answer, not a failure.
		out := make([]SampleState, 0, len(s.man.Samples))
		for _, name := range s.man.Samples {
			out = append(out, SampleState{SampleID: name, State: StateNotAssayed})
		}
		return out, nil
	}
	return st.Classify(l, g)
}

// ClassifyMany resolves a whole locus set, asking every member concerned at the
// same time.
//
// THE SHAPE THIS EXISTS FOR. A panel spanning fifteen chromosomes is fifteen
// independent reads of fifteen different files, and doing them one after
// another wastes exactly the parallelism the layout was arranged to permit --
// especially over object storage, where each is dominated by latency rather
// than by work.
func (s *VarSet) ClassifyMany(loci []Locus, g Gate) (map[Locus][]SampleState, error) {
	out := make(map[Locus][]SampleState, len(loci))
	if len(loci) == 0 {
		return out, nil
	}

	// Group by member first, so each is asked once for everything it holds.
	byMember := map[int][]Locus{}
	var unheld []Locus
	for _, l := range loci {
		if i, ok := s.byChrom[CanonKey(l.Chrom)]; ok {
			byMember[i] = append(byMember[i], l)
			continue
		}
		unheld = append(unheld, l)
	}

	// A chromosome no member holds is not assayed for everyone. Answered
	// without opening anything.
	for _, l := range unheld {
		states := make([]SampleState, 0, len(s.man.Samples))
		for _, name := range s.man.Samples {
			states = append(states, SampleState{SampleID: name, State: StateNotAssayed})
		}
		out[l] = states
	}
	if len(byMember) == 0 {
		return out, nil
	}

	idxs := make([]int, 0, len(byMember))
	for i := range byMember {
		idxs = append(idxs, i)
	}
	sort.Ints(idxs)

	// Bounded, because a whole-genome panel would otherwise open twenty-five
	// stores at once and each of those fans out across its own shards. The cap
	// is on members; the shard fan-out inside each is capped separately.
	workers := runtime.GOMAXPROCS(0)
	if workers > len(idxs) {
		workers = len(idxs)
	}

	results := make([]map[Locus][]SampleState, len(idxs))
	errs := make([]error, len(idxs))
	sem := make(chan struct{}, workers)
	var wg sync.WaitGroup
	for n, i := range idxs {
		wg.Add(1)
		go func(n, i int) {
			defer wg.Done()
			sem <- struct{}{}
			defer func() { <-sem }()
			st, err := s.openMemberAt(i)
			if err != nil {
				errs[n] = err
				return
			}
			results[n], errs[n] = st.ClassifyMany(byMember[i], g)
		}(n, i)
	}
	wg.Wait()

	for n, err := range errs {
		if err != nil {
			return nil, fmt.Errorf("querying set member %s: %w", s.man.Members[idxs[n]].Name, err)
		}
	}
	// Merged in member order, so a result never depends on which finished first.
	for _, r := range results {
		for l, states := range r {
			out[l] = states
		}
	}
	return out, nil
}

// Sites walks every member's catalog, in member order.
func (s *VarSet) Sites(fn func(Site) bool) error {
	for i := range s.man.Members {
		st, err := s.openMemberAt(i)
		if err != nil {
			return err
		}
		stop := false
		if err := st.Sites(func(x Site) bool {
			if !fn(x) {
				stop = true
				return false
			}
			return true
		}); err != nil {
			return err
		}
		if stop {
			return nil
		}
	}
	return nil
}

// Regions walks every member's callable runs, in member order.
func (s *VarSet) Regions(fn func(CalledSiteRun) bool) error {
	for i := range s.man.Members {
		st, err := s.openMemberAt(i)
		if err != nil {
			return err
		}
		stop := false
		if err := st.Regions(func(r CalledSiteRun) bool {
			if !fn(r) {
				stop = true
				return false
			}
			return true
		}); err != nil {
			return err
		}
		if stop {
			return nil
		}
	}
	return nil
}

// Site returns one catalog entry, from the member holding its chromosome.
func (s *VarSet) Site(l Locus) (Site, bool, error) {
	st, ok, err := s.member(l.Chrom)
	if err != nil || !ok {
		return Site{}, false, err
	}
	return st.Site(l)
}

// Calls streams the genotypes a query selects, member by member.
//
// SEQUENTIAL AND ORDERED, unlike ClassifyMany. The contract is that rows arrive
// in the store's own order -- contig order, then position, then ALT, then
// sample -- and a set keeps that by visiting its members in the order they were
// declared. Racing them would produce rows faster and in an order no caller
// could rely on, and the callers that stream are the ones writing a VCF.
//
// A query naming no site selector reads every member; one naming loci or spans
// reads only the members holding those chromosomes.
func (s *VarSet) Calls(q Query) (iter.Seq2[Call, error], error) {
	want := s.membersFor(q)

	// Setup failures surface here rather than mid-iteration, so a caller learns
	// about an unusable query before it starts reading. That means opening the
	// members up front -- which is the one place a set is not lazy, and it is
	// the contract that requires it.
	stores := make([]*ParquetStore, 0, len(want))
	for _, i := range want {
		st, err := s.openMemberAt(i)
		if err != nil {
			return nil, err
		}
		stores = append(stores, st)
	}

	seqs := make([]iter.Seq2[Call, error], 0, len(stores))
	for _, st := range stores {
		seq, err := st.Calls(q)
		if err != nil {
			return nil, err
		}
		seqs = append(seqs, seq)
	}

	return func(yield func(Call, error) bool) {
		for _, seq := range seqs {
			for c, err := range seq {
				if !yield(c, err) {
					return
				}
			}
		}
	}, nil
}

// membersFor returns the members a query can touch, in declared order.
func (s *VarSet) membersFor(q Query) []int {
	if len(q.Loci) == 0 && len(q.Spans) == 0 {
		out := make([]int, len(s.man.Members))
		for i := range out {
			out[i] = i
		}
		return out
	}
	hit := map[int]bool{}
	for _, l := range q.Loci {
		if i, ok := s.byChrom[CanonKey(l.Chrom)]; ok {
			hit[i] = true
		}
	}
	for _, sp := range q.Spans {
		if i, ok := s.byChrom[CanonKey(sp.Chrom)]; ok {
			hit[i] = true
		}
	}
	out := make([]int, 0, len(hit))
	for i := range hit {
		out = append(out, i)
	}
	sort.Ints(out)
	return out
}
