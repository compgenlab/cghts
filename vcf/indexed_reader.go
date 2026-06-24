package vcf

import (
	"iter"

	"github.com/compgenlab/hts/htsio/tabix"
)

// IndexedVcfReader provides random access to a tabix-indexed VCF file
// (BGZF-compressed with a companion .tbi or .csi index).
type IndexedVcfReader struct {
	tr       *tabix.Reader
	filename string
	header   *VcfHeader
}

// NewIndexedVcfReader opens a tabix-indexed VCF file for random access. The file
// must be BGZF-compressed and have a companion .tbi or .csi index. The caller
// should [IndexedVcfReader.Close] the reader when done.
func NewIndexedVcfReader(filename string) (*IndexedVcfReader, error) {
	tr, err := tabix.NewReader(filename)
	if err != nil {
		return nil, err
	}
	return &IndexedVcfReader{tr: tr, filename: filename}, nil
}

// Header returns the VCF header, reading it from the start of the BGZF stream on
// first call.
func (r *IndexedVcfReader) Header() (*VcfHeader, error) {
	if r.header == nil {
		hr, err := NewVcfFile(r.filename)
		if err != nil {
			return nil, err
		}
		defer hr.Close()
		h, err := hr.Header()
		if err != nil {
			return nil, err
		}
		r.header = h
	}
	return r.header, nil
}

// Query returns an iterator over the VCF records overlapping the 0-based
// half-open region [start, end) on the given reference. The iterator yields
// (nil, err) and stops if a record line cannot be parsed.
func (r *IndexedVcfReader) Query(ref string, start, end int) (iter.Seq2[*VcfRecord, error], error) {
	header, err := r.Header()
	if err != nil {
		return nil, err
	}
	recs, err := r.tr.Query(ref, start, end)
	if err != nil {
		return nil, err
	}
	return func(yield func(*VcfRecord, error) bool) {
		for rec, err := range recs {
			if err != nil {
				yield(nil, err)
				return
			}
			vr, perr := newRecord(rec.Line, header)
			if perr != nil {
				yield(nil, perr)
				return
			}
			if !yield(vr, nil) {
				return
			}
		}
	}, nil
}

// Close releases resources held by the reader.
func (r *IndexedVcfReader) Close() error {
	return r.tr.Close()
}
