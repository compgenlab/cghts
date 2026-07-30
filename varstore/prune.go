package varstore

import (
	"io"
	"math"
	"os"

	"github.com/parquet-go/parquet-go"
)

// Row-group pruning. Parquet records min/max statistics per column per row
// group, and vcf-toparquet writes rows in VCF order, which is coordinate
// sorted -- so those ranges are already tight and non-overlapping. Reading them
// turns a whole-file scan into a scan of the few groups that could match.
//
// Every predicate here must be conservative in one direction: it may keep a row
// group that turns out to hold nothing, but it must never skip one that holds a
// match. The per-row checks downstream stay exactly as they were, so pruning can
// only affect how much is read, never what is found.

// rowGroupFilter decides whether a row group could contain a match.
type rowGroupFilter func(rg parquet.RowGroup) bool

// keepAll is the filter for an unrestricted scan.
func keepAll(parquet.RowGroup) bool { return true }

// colBounds returns the min and max of a named column within a row group,
// taken across its pages. ok is false when the column carries no statistics, in
// which case the caller must keep the group.
func colBounds(rg parquet.RowGroup, name string) (min, max parquet.Value, ok bool) {
	idx := -1
	for i, path := range rg.Schema().Columns() {
		if len(path) == 1 && path[0] == name {
			idx = i
			break
		}
	}
	if idx < 0 {
		return min, max, false
	}
	chunks := rg.ColumnChunks()
	if idx >= len(chunks) {
		return min, max, false
	}
	chunk := chunks[idx]
	ci, err := chunk.ColumnIndex()
	if err != nil || ci == nil || ci.NumPages() == 0 {
		return min, max, false
	}
	cmp := chunk.Type().Compare
	for p := 0; p < ci.NumPages(); p++ {
		if ci.NullPage(p) {
			continue
		}
		lo, hi := ci.MinValue(p), ci.MaxValue(p)
		if lo.IsNull() || hi.IsNull() {
			continue
		}
		if !ok {
			min, max, ok = lo, hi, true
			continue
		}
		if cmp(lo, min) < 0 {
			min = lo
		}
		if cmp(hi, max) > 0 {
			max = hi
		}
	}
	return min, max, ok
}

// locusFilter keeps only row groups that could hold a given locus.
//
// Position does the real work: it is numeric, high-cardinality and sorted, so
// its bounds are tight. Chromosome is only used when a group holds exactly one,
// because the statistics are lexicographic while chromosome names are not
// ordered that way -- comparing "chr17" against a chr1..chr9 range would be
// meaningless. Restricting the test to single-chromosome groups keeps it sound
// and still catches the common case, since coordinate-sorted input puts almost
// every group inside one chromosome.
func locusFilter(l Locus) rowGroupFilter {
	return func(rg parquet.RowGroup) bool {
		if lo, hi, ok := colBounds(rg, "pos"); ok {
			if int64(l.Pos) < lo.Int64() || int64(l.Pos) > hi.Int64() {
				return false
			}
		}
		if lo, hi, ok := colBounds(rg, "chrom"); ok {
			a, b := lo.String(), hi.String()
			if a == b && !SameChrom(a, l.Chrom) {
				return false
			}
		}
		return true
	}
}

// spanFilter keeps row groups overlapping a 0-based half-open span.
func spanFilter(s *Span) rowGroupFilter {
	if s == nil {
		return keepAll
	}
	return func(rg parquet.RowGroup) bool {
		if lo, hi, ok := colBounds(rg, "pos"); ok {
			// Span is 0-based half-open; stored positions are 1-based.
			if int64(s.End) < lo.Int64() || int64(s.Start)+1 > hi.Int64() {
				return false
			}
		}
		if lo, hi, ok := colBounds(rg, "chrom"); ok {
			a, b := lo.String(), hi.String()
			if a == b && !SameChrom(a, s.Chrom) {
				return false
			}
		}
		return true
	}
}

// coveringFilter keeps row groups whose runs could span a position. A run's
// End is not sorted, so only Start can bound the search: any run covering pos
// must begin at or before it.
func coveringFilter(chrom string, pos int32) rowGroupFilter {
	return func(rg parquet.RowGroup) bool {
		if lo, _, ok := colBounds(rg, "start"); ok {
			if int64(pos) < lo.Int64() {
				return false
			}
		}
		if lo, hi, ok := colBounds(rg, "chrom"); ok {
			a, b := lo.String(), hi.String()
			if a == b && !SameChrom(a, chrom) {
				return false
			}
		}
		return true
	}
}

// spanRunFilter keeps row groups whose runs could overlap a span.
//
// Like coveringFilter, only Start can bound the search: End is not sorted, so a
// group whose runs all start before the span may still hold one reaching into
// it. Only the other side is safe -- a group whose earliest run begins after the
// span's last position cannot overlap it at all.
func spanRunFilter(s *Span) rowGroupFilter {
	if s == nil {
		return keepAll
	}
	return func(rg parquet.RowGroup) bool {
		if lo, _, ok := colBounds(rg, "start"); ok {
			// Span is 0-based half-open, so its last 1-based position is End.
			if int64(s.End) < lo.Int64() {
				return false
			}
		}
		if lo, hi, ok := colBounds(rg, "chrom"); ok {
			a, b := lo.String(), hi.String()
			if a == b && !SameChrom(a, s.Chrom) {
				return false
			}
		}
		return true
	}
}

// callsFilter turns a query's site axis into a row-group predicate for the calls
// and sites files, which share a (chrom, pos) layout.
//
// Sample selection deliberately does not participate, and should not be added:
// every sample appears in nearly every row group, so the sample_id bloom filter
// measured 429ms against 426ms, and a union over several samples is weaker still.
func callsFilter(q Query) rowGroupFilter {
	switch {
	case len(q.Loci) == 0 && len(q.Spans) == 0:
		return keepAll
	case len(q.Loci) == 1 && len(q.Spans) == 0:
		return locusFilter(q.Loci[0])
	case len(q.Loci) == 0 && len(q.Spans) == 1:
		return spanFilter(&q.Spans[0])
	}
	return unionFilter(q)
}

// siteScanFilter prunes the sites file for a query.
func siteScanFilter(q Query) rowGroupFilter { return callsFilter(q) }

// runScanFilter prunes the regions file, whose position column is "start" rather
// than "pos" and whose "end" is unsorted -- so only a lower bound on start is
// usable, and with several selectors there is little left to prune on.
func runScanFilter(q Query) rowGroupFilter {
	switch {
	case len(q.Loci) == 1 && len(q.Spans) == 0:
		return coveringFilter(q.Loci[0].Chrom, q.Loci[0].Pos)
	case len(q.Loci) == 0 && len(q.Spans) == 1:
		return spanRunFilter(&q.Spans[0])
	}
	return keepAll
}

// unionFilter keeps a row group that any of the query's site selectors admits.
//
// It must stay a union: skipping a group because one selector rejects it would
// drop matches for another. Bounds are collapsed across all selectors up front
// rather than evaluated per selector per group, so a million-locus panel costs one
// pass over the panel instead of a million predicate calls for every row group.
// That makes it conservative -- a group between two distant selectors is kept --
// which is the safe direction.
func unionFilter(q Query) rowGroupFilter {
	lo, hi := int32(math.MaxInt32), int32(0)
	chroms := map[string]bool{}
	note := func(chrom string, first, last int32) {
		if first < lo {
			lo = first
		}
		if last > hi {
			hi = last
		}
		chroms[CanonKey(chrom)] = true
	}
	for _, l := range q.Loci {
		note(l.Chrom, l.Pos, l.Pos)
	}
	for _, sp := range q.Spans {
		// Spans are 0-based half-open; stored positions are 1-based.
		note(sp.Chrom, sp.Start+1, sp.End)
	}
	return func(rg parquet.RowGroup) bool {
		if a, b, ok := colBounds(rg, "pos"); ok {
			if int64(hi) < a.Int64() || int64(lo) > b.Int64() {
				return false
			}
		}
		if a, b, ok := colBounds(rg, "chrom"); ok {
			x, y := a.String(), b.String()
			if x == y && !chroms[CanonKey(x)] {
				return false
			}
		}
		return true
	}
}

// eachRowGroup calls fn for every row group in path, without reading any rows.
func eachRowGroup(path string, fn func(parquet.RowGroup)) error {
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()
	st, err := f.Stat()
	if err != nil {
		return err
	}
	pf, err := parquet.OpenFile(f, st.Size())
	if err != nil {
		return err
	}
	for _, rg := range pf.RowGroups() {
		fn(rg)
	}
	return nil
}

// scanParquetPruned walks path, reading only the row groups keep admits.
//
// Rows are addressed by seeking the generic reader to each surviving group's
// offset, so decoding still goes through the same typed path as a full scan.
func scanParquetPruned[T any](path string, keep rowGroupFilter, fn func(T) bool) error {
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()
	st, err := f.Stat()
	if err != nil {
		return err
	}
	pf, err := parquet.OpenFile(f, st.Size())
	if err != nil {
		return err
	}

	r := parquet.NewGenericReader[T](f)
	defer r.Close()

	buf := make([]T, 1024)
	var offset int64
	for _, rg := range pf.RowGroups() {
		n := rg.NumRows()
		if !keep(rg) {
			offset += n
			continue
		}
		if err := r.SeekToRow(offset); err != nil {
			return err
		}
		for remaining := n; remaining > 0; {
			want := int64(len(buf))
			if remaining < want {
				want = remaining
			}
			got, err := r.Read(buf[:want])
			for i := 0; i < got; i++ {
				if !fn(buf[i]) {
					return nil
				}
			}
			remaining -= int64(got)
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if got == 0 {
				break
			}
		}
		offset += n
	}
	return nil
}
