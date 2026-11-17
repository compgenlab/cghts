package tabix

import (
	"encoding/binary"
	"errors"
	"io"

	"github.com/compgenlab/cghts/htsio/bgzf"
)

// RefStats is what an index records about one reference beyond its bins: the
// optional "pseudo-bin" of the SAM spec's BAI section, which TBI and CSI share.
// htslib writes it for every index it builds, BAM and tabix alike, whatever the
// preset or columns; cghts writes it since v0.12.9.
type RefStats struct {
	// Mapped is the number of records indexed on the reference. For a tabix
	// index that is every data line on it: header (skip) and meta lines are
	// not counted. For a BAI it is the reads with a position that are mapped.
	Mapped uint64
	// Unmapped is a BAI's count of unmapped reads placed on the reference
	// (normally beside a mapped mate). Tabix always writes 0.
	Unmapped uint64
	// Begin is the virtual offset of the reference's first record and End the
	// offset just past its last, so [Begin, End) holds all of its records.
	Begin, End bgzf.VirtualOffset
}

// tbiPseudoBin is the pseudo-bin's number in a BAI or TBI: one past the last
// real bin of the 6-level scheme (37449, the last 16 kb leaf below 2^29).
const tbiPseudoBin = 37450

// tbiMaxPos is the first position a BAI/TBI cannot bin: 2^29.
const tbiMaxPos = 1 << 29

// csiPseudoBin is the pseudo-bin's number in a CSI of the given depth, again
// one past the last real bin. ok is false for a depth too deep to bin at all.
func csiPseudoBin(depth int32) (uint32, bool) {
	if depth < 0 || depth > 9 { // 3*(depth+1) bits must fit a bin number
		return 0, false
	}
	return uint32(((1<<(3*(depth+1)))-1)/7 + 1), true
}

// parsePseudoBin reads a pseudo-bin's two chunks into RefStats. It is not a bin
// -- its number lies outside the binning scheme and its second chunk holds
// counts, not offsets -- so it must never reach a query as one.
func parsePseudoBin(chunks []Chunk) (*RefStats, bool) {
	if len(chunks) != 2 {
		return nil, false
	}
	return &RefStats{
		Begin:    chunks[0].Begin,
		End:      chunks[0].End,
		Mapped:   uint64(chunks[1].Begin),
		Unmapped: uint64(chunks[1].End),
	}, true
}

// readNoCoor reads the optional n_no_coor that may follow the last reference
// in a BAI, TBI or CSI. It is absent from older indexes, so a clean EOF means
// "not recorded", not an error.
func readNoCoor(r io.Reader) (*uint64, error) {
	var n uint64
	if err := binary.Read(r, binary.LittleEndian, &n); err != nil {
		if errors.Is(err, io.EOF) || errors.Is(err, io.ErrUnexpectedEOF) {
			return nil, nil
		}
		return nil, err
	}
	return &n, nil
}

// refStats answers RefStats for one reference of either index type. A
// reference with no bins at all holds no records, which is known even without a
// pseudo-bin; one with bins but no pseudo-bin is unknown.
func refStats(stats *RefStats, nBins int) (RefStats, bool) {
	switch {
	case stats != nil:
		return *stats, true
	case nBins == 0:
		return RefStats{}, true
	}
	return RefStats{}, false
}

// RefStats returns the counts and record extent the index records for the
// reference with the given ID. ok is false when the index does not record them
// -- a cghts index from before v0.12.9, or a tabix older than htslib 1.0 --
// which is "unknown", never zero.
func (idx *BinIndex) RefStats(refID int) (RefStats, bool) {
	if refID < 0 || refID >= len(idx.refs) {
		return RefStats{}, false
	}
	ref := &idx.refs[refID]
	return refStats(ref.stats, len(ref.bins))
}

// NoCoor returns the index's trailing n_no_coor -- a BAI's count of reads with
// no position; 0 in a tabix index -- and whether the index records one.
func (idx *BinIndex) NoCoor() (uint64, bool) {
	if idx.noCoor == nil {
		return 0, false
	}
	return *idx.noCoor, true
}

// RefStats is BinIndex.RefStats for a CSI.
func (idx *CSIIndex) RefStats(refID int) (RefStats, bool) {
	if refID < 0 || refID >= len(idx.refs) {
		return RefStats{}, false
	}
	ref := &idx.refs[refID]
	return refStats(ref.stats, len(ref.bins))
}

// NoCoor is BinIndex.NoCoor for a CSI.
func (idx *CSIIndex) NoCoor() (uint64, bool) {
	if idx.noCoor == nil {
		return 0, false
	}
	return *idx.noCoor, true
}

// RecordCounts returns the number of records the index says each reference
// holds, in RefNames order, without reading the data. ok is false if the index
// does not record a count for every reference that has records (see
// BinIndex.RefStats); the counts are then partial and should not be summed.
func (tr *Reader) RecordCounts() (counts []uint64, ok bool) {
	names := tr.idx.RefNames()
	counts = make([]uint64, len(names))
	ok = true
	for i := range names {
		st, known := tr.idx.RefStats(i)
		if !known {
			ok = false
			continue
		}
		counts[i] = st.Mapped
	}
	return counts, ok
}
