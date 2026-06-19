package bed

import (
	"io"
	"iter"
	"sort"
)

// SetOp selects the set-algebra operation performed by [SetOperation].
type SetOp int

const (
	// OpInter keeps bases covered by both A and B.
	OpInter SetOp = iota
	// OpUnion keeps bases covered by A or B (merged).
	OpUnion
	// OpSub keeps bases covered by A but not B (A − B).
	OpSub
	// OpXor keeps bases covered by exactly one of A or B (symmetric difference).
	OpXor
)

// SetOpts configures a [SetOperation].
type SetOpts struct {
	Op SetOp
	// IgnoreStrand collapses all records into a single strand-agnostic group;
	// otherwise the operation is run independently per strand value.
	IgnoreStrand bool
	// FlankX (OpUnion only) merges regions that are within X bases of each other
	// (gap <= FlankX); output coordinates remain the true unpadded extents.
	FlankX int
}

// Contrib identifies an input region that overlaps an output segment.
type Contrib struct {
	Source int // 0 = A, 1 = B
	Name   string
	Score  float64
}

// OutSegment is one emitted, merged, non-overlapping interval plus the input
// regions that cover it. Coordinates are 0-based half-open.
type OutSegment struct {
	Ref      string
	Start    int
	End      int
	Strand   Strand
	Contribs []Contrib
}

// SetOperation streams the result of a coordinate set-algebra operation over two
// BED readers as merged, non-overlapping [OutSegment] values in sorted order.
//
// Both inputs are read fully and grouped by chromosome; the operation is then
// applied per chromosome (and, unless IgnoreStrand is set, per strand). Inputs
// need not be pre-sorted — events are sorted internally — but output is always
// emitted sorted by (start, end). The chromosome output order is the order of
// first appearance in A, then chromosomes seen only in B.
func SetOperation(a, b *BedReader, opts SetOpts) iter.Seq2[*OutSegment, error] {
	return func(yield func(*OutSegment, error) bool) {
		aRecs, aOrder, err := readAllByChrom(a)
		if err != nil {
			yield(nil, err)
			return
		}
		bRecs, bOrder, err := readAllByChrom(b)
		if err != nil {
			yield(nil, err)
			return
		}
		emitSegments(aRecs, aOrder, bRecs, bOrder, opts, func(s *OutSegment) bool {
			return yield(s, nil)
		})
	}
}

// SetOperationRecords is like [SetOperation] but operates on already-read record
// slices (e.g. when the caller has inspected the inputs first). It cannot fail.
func SetOperationRecords(aRecs, bRecs []*BedRecord, opts SetOpts) iter.Seq[*OutSegment] {
	am, aOrder := groupByChrom(aRecs)
	bm, bOrder := groupByChrom(bRecs)
	return func(yield func(*OutSegment) bool) {
		emitSegments(am, aOrder, bm, bOrder, opts, yield)
	}
}

// emitSegments runs the operation per chromosome (A's chrom order, then B-only
// chroms) and passes each output segment to emit; it stops if emit returns false.
func emitSegments(aRecs map[string][]*BedRecord, aOrder []string, bRecs map[string][]*BedRecord, bOrder []string, opts SetOpts, emit func(*OutSegment) bool) {
	order := append([]string(nil), aOrder...)
	seen := make(map[string]bool, len(aOrder))
	for _, c := range aOrder {
		seen[c] = true
	}
	for _, c := range bOrder {
		if !seen[c] {
			seen[c] = true
			order = append(order, c)
		}
	}
	for _, chrom := range order {
		for _, seg := range processChrom(chrom, aRecs[chrom], bRecs[chrom], opts) {
			if !emit(seg) {
				return
			}
		}
	}
}

func groupByChrom(recs []*BedRecord) (map[string][]*BedRecord, []string) {
	m := map[string][]*BedRecord{}
	var order []string
	for _, rec := range recs {
		if _, ok := m[rec.Ref]; !ok {
			order = append(order, rec.Ref)
		}
		m[rec.Ref] = append(m[rec.Ref], rec)
	}
	return m, order
}

func readAllByChrom(r *BedReader) (map[string][]*BedRecord, []string, error) {
	var recs []*BedRecord
	for {
		rec, err := r.NextRecord()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, nil, err
		}
		recs = append(recs, rec)
	}
	m, order := groupByChrom(recs)
	return m, order, nil
}

func processChrom(chrom string, aRecs, bRecs []*BedRecord, opts SetOpts) []*OutSegment {
	var segs []*OutSegment

	if opts.IgnoreStrand {
		segs = processGroup(chrom, StrandNone, aRecs, bRecs, opts)
	} else {
		var strands []Strand
		seen := map[Strand]bool{}
		collect := func(recs []*BedRecord) {
			for _, r := range recs {
				if !seen[r.Strand] {
					seen[r.Strand] = true
					strands = append(strands, r.Strand)
				}
			}
		}
		collect(aRecs)
		collect(bRecs)
		for _, s := range strands {
			segs = append(segs, processGroup(chrom, s, filterStrand(aRecs, s), filterStrand(bRecs, s), opts)...)
		}
	}

	sort.Slice(segs, func(i, j int) bool {
		if segs[i].Start != segs[j].Start {
			return segs[i].Start < segs[j].Start
		}
		if segs[i].End != segs[j].End {
			return segs[i].End < segs[j].End
		}
		return segs[i].Strand < segs[j].Strand
	})
	return segs
}

func filterStrand(recs []*BedRecord, s Strand) []*BedRecord {
	var out []*BedRecord
	for _, r := range recs {
		if r.Strand == s {
			out = append(out, r)
		}
	}
	return out
}

type iv struct {
	start, end int
}

func processGroup(chrom string, strand Strand, aRecs, bRecs []*BedRecord, opts SetOpts) []*OutSegment {
	var ivs []iv
	if opts.Op == OpUnion {
		ivs = unionIntervals(aRecs, bRecs, opts.FlankX)
	} else {
		ivs = sweepIntervals(aRecs, bRecs, opts.Op)
	}

	out := make([]*OutSegment, 0, len(ivs))
	for _, v := range ivs {
		out = append(out, &OutSegment{
			Ref:      chrom,
			Start:    v.start,
			End:      v.end,
			Strand:   strand,
			Contribs: collectContribs(v, aRecs, bRecs),
		})
	}
	return out
}

// unionIntervals merges all (non-empty) A and B intervals, bridging regions
// within flankX bases of each other (gap <= flankX). Output coordinates are the
// true min-start..max-end per cluster.
func unionIntervals(aRecs, bRecs []*BedRecord, flankX int) []iv {
	var list []iv
	add := func(recs []*BedRecord) {
		for _, r := range recs {
			if r.End > r.Start {
				list = append(list, iv{r.Start, r.End})
			}
		}
	}
	add(aRecs)
	add(bRecs)
	if len(list) == 0 {
		return nil
	}
	sort.Slice(list, func(i, j int) bool {
		if list[i].start != list[j].start {
			return list[i].start < list[j].start
		}
		return list[i].end < list[j].end
	})

	out := []iv{list[0]}
	for _, v := range list[1:] {
		cur := &out[len(out)-1]
		if v.start <= cur.end+flankX {
			if v.end > cur.end {
				cur.end = v.end
			}
		} else {
			out = append(out, v)
		}
	}
	return out
}

type event struct {
	pos    int
	dA, dB int
}

// sweepIntervals computes the inter/sub/xor output intervals via a coverage
// sweep with two independent counters.
func sweepIntervals(aRecs, bRecs []*BedRecord, op SetOp) []iv {
	var evs []event
	for _, r := range aRecs {
		if r.End > r.Start {
			evs = append(evs, event{r.Start, 1, 0}, event{r.End, -1, 0})
		}
	}
	for _, r := range bRecs {
		if r.End > r.Start {
			evs = append(evs, event{r.Start, 0, 1}, event{r.End, 0, -1})
		}
	}
	if len(evs) == 0 {
		return nil
	}
	sort.Slice(evs, func(i, j int) bool { return evs[i].pos < evs[j].pos })

	var out []iv
	countA, countB := 0, 0
	i, n := 0, len(evs)
	for i < n {
		p := evs[i].pos
		for i < n && evs[i].pos == p {
			countA += evs[i].dA
			countB += evs[i].dB
			i++
		}
		if i >= n {
			break
		}
		nextP := evs[i].pos
		if nextP > p && setPredicate(op, countA, countB) {
			if len(out) > 0 && out[len(out)-1].end == p {
				out[len(out)-1].end = nextP // coalesce adjacent emitted spans
			} else {
				out = append(out, iv{p, nextP})
			}
		}
	}
	return out
}

func setPredicate(op SetOp, countA, countB int) bool {
	switch op {
	case OpInter:
		return countA > 0 && countB > 0
	case OpSub:
		return countA > 0 && countB == 0
	case OpXor:
		return (countA > 0) != (countB > 0)
	}
	return false
}

// collectContribs returns every non-empty input region overlapping segment v.
// For A−B this naturally yields only A regions (B does not overlap A−B output);
// for inter it yields both sides; for xor only the covering side.
func collectContribs(v iv, aRecs, bRecs []*BedRecord) []Contrib {
	var c []Contrib
	scan := func(recs []*BedRecord, source int) {
		for _, r := range recs {
			if r.End > r.Start && r.Start < v.end && r.End > v.start {
				c = append(c, Contrib{Source: source, Name: r.Name, Score: r.Score})
			}
		}
	}
	scan(aRecs, 0)
	scan(bRecs, 1)
	return c
}
