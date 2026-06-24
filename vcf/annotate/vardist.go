package annotate

import (
	"io"
	"strconv"

	"github.com/compgenlab/hts/vcf"
)

// VariantDistance annotates each variant with CG_VARDIST, the distance to its
// nearest neighboring variant on the same chromosome (sorted input assumed). It
// is a look-ahead annotator: it must see the next variant before emitting the
// current one, so it wraps the record stream. Isolated variants get -1. Ports
// ngsutilsj VariantDistance.
type VariantDistance struct{ closeNoop }

// NewVariantDistance returns a VariantDistance annotator (--vardist).
func NewVariantDistance() *VariantDistance { return &VariantDistance{} }

// SetupHeader declares the CG_VARDIST INFO field. (ngsutilsj registers this as a
// FORMAT def, a bug; cgio registers it as INFO.)
func (a *VariantDistance) SetupHeader(h *vcf.VcfHeader) error {
	h.AddInfo(infoDef("CG_VARDIST", "1", "Integer", "Distance to the nearest variant (absolute value)"))
	return nil
}

// Wrap returns a source that emits each record annotated with its distance to
// the nearest neighbor. It buffers one record (holding pointers, not copies).
func (a *VariantDistance) Wrap(next Source) Source {
	var last, cur *vcf.VcfRecord
	var lastDist int64 = -1
	done := false

	annotate := func(rec *vcf.VcfRecord, d int64) {
		rec.AddInfo("CG_VARDIST", strconv.FormatInt(d, 10))
	}

	return func() (*vcf.VcfRecord, error) {
		if done {
			return nil, io.EOF
		}
		if cur == nil {
			c, err := pull(next)
			if err != nil {
				return nil, err
			}
			cur = c
			if cur == nil {
				return nil, io.EOF
			}
			if last != nil && cur.Chrom != last.Chrom {
				annotate(last, lastDist)
				lastDist = -1
				return last, nil
			}
		}

		if last == nil {
			lastDist = -1
		}
		last = cur
		c, err := pull(next)
		if err != nil {
			return nil, err
		}
		cur = c
		if cur == nil {
			done = true
		}
		if cur == nil || cur.Chrom != last.Chrom {
			annotate(last, lastDist)
			lastDist = -1
			return last, nil
		}

		dist := int64(cur.Pos - last.Pos)
		if lastDist > -1 && lastDist < dist {
			annotate(last, lastDist)
		} else {
			annotate(last, dist)
		}
		lastDist = dist
		return last, nil
	}
}

// pull reads one record from next, mapping io.EOF to a nil record.
func pull(next Source) (*vcf.VcfRecord, error) {
	rec, err := next()
	if err == io.EOF {
		return nil, nil
	}
	if err != nil {
		return nil, err
	}
	return rec, nil
}
