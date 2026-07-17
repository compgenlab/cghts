package bed

import (
	"iter"

	"github.com/compgenlab/cghts/htsio/tabix"
)

// IndexedBedReader provides random access to a tabix-indexed BED file
// (BGZF-compressed with a companion .tbi or .csi index).
type IndexedBedReader struct {
	tr           *tabix.Reader
	ignoreStrand bool
}

// NewIndexedBedReader opens a tabix-indexed BED file for random access. The
// file must be BGZF-compressed and have a companion .tbi or .csi index. The
// caller should [IndexedBedReader.Close] the reader when done.
func NewIndexedBedReader(filename string) (*IndexedBedReader, error) {
	tr, err := tabix.NewReader(filename)
	if err != nil {
		return nil, err
	}
	return &IndexedBedReader{tr: tr}, nil
}

// IgnoreStrand configures the reader to leave every record's strand
// unspecified ([StrandNone]) regardless of the strand column. It returns the
// reader to allow chaining at construction.
func (r *IndexedBedReader) IgnoreStrand(b bool) *IndexedBedReader {
	r.ignoreStrand = b
	return r
}

// Query returns an iterator over the BED records overlapping the 0-based
// half-open region [start, end) on the given reference. The iterator yields
// (nil, err) and stops if a record line cannot be parsed.
func (r *IndexedBedReader) Query(ref string, start, end int) (iter.Seq2[*BedRecord, error], error) {
	recs, err := r.tr.Query(ref, start, end)
	if err != nil {
		return nil, err
	}
	return func(yield func(*BedRecord, error) bool) {
		for rec, err := range recs {
			if err != nil {
				yield(nil, err)
				return
			}
			bedRec, ok, perr := parseBedLine(rec.Line, r.ignoreStrand)
			if perr != nil {
				yield(nil, perr)
				return
			}
			if !ok {
				continue
			}
			if !yield(bedRec, nil) {
				return
			}
		}
	}, nil
}

// Close releases resources held by the reader.
func (r *IndexedBedReader) Close() error {
	return r.tr.Close()
}
