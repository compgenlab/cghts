package varstore

import (
	"io"
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

// sampleFilter keeps row groups that may contain a sample, consulting the
// bloom filter on sample_id when the writer produced one. Sample ids are
// high-cardinality and unsorted, so statistics are useless here but a bloom
// filter is exact about absence.
func sampleFilter(sample string) rowGroupFilter {
	return func(rg parquet.RowGroup) bool {
		for i, path := range rg.Schema().Columns() {
			if len(path) != 1 || path[0] != "sample_id" {
				continue
			}
			chunks := rg.ColumnChunks()
			if i >= len(chunks) {
				return true
			}
			bf := chunks[i].BloomFilter()
			if bf == nil {
				return true
			}
			ok, err := bf.Check(parquet.ByteArrayValue([]byte(sample)))
			if err != nil {
				return true
			}
			return ok
		}
		return true
	}
}

// bothFilters keeps a row group only when both filters do.
func bothFilters(a, b rowGroupFilter) rowGroupFilter {
	return func(rg parquet.RowGroup) bool { return a(rg) && b(rg) }
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
