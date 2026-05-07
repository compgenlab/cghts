package tabix

import (
	"encoding/binary"
	"fmt"
	"io"
	"os"
	"sort"

	"github.com/compgen-io/cgkit/htsio/bgzf"
)

// csiBin holds a single bin from a CSI index, including its loffset
// (minimum virtual offset for the bin, replacing the linear index).
type csiBin struct {
	bin     uint32
	loffset bgzf.VirtualOffset
	chunks  []Chunk
}

// csiRefIndex holds the bins for one reference sequence in a CSI index.
type csiRefIndex struct {
	bins map[uint32]*csiBin
}

// CSIIndex is a parsed CSI (coordinate-sorted index) file. CSI supports
// variable binning depth and min_shift, enabling indexing of sequences
// longer than 512 Mb.
type CSIIndex struct {
	minShift int32
	depth    int32
	refs     []csiRefIndex

	// Tabix metadata from the auxiliary data block.
	Format    int32
	ColSeq    int32
	ColBeg    int32
	ColEnd    int32
	Meta      int32
	Skip      int32
	ZeroBased bool
	Names     []string
}

// LoadCSI reads a CSI index file. CSI files are BGZF-compressed.
func LoadCSI(filename string) (*CSIIndex, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	r := bgzf.NewReader(f)

	var magic [4]byte
	if _, err := io.ReadFull(r, magic[:]); err != nil {
		return nil, fmt.Errorf("csi: reading magic: %w", err)
	}
	if magic != [4]byte{'C', 'S', 'I', 1} {
		return nil, fmt.Errorf("csi: invalid magic: %x", magic)
	}

	idx := &CSIIndex{}

	if err := binary.Read(r, binary.LittleEndian, &idx.minShift); err != nil {
		return nil, fmt.Errorf("csi: reading min_shift: %w", err)
	}
	if err := binary.Read(r, binary.LittleEndian, &idx.depth); err != nil {
		return nil, fmt.Errorf("csi: reading depth: %w", err)
	}

	var lAux int32
	if err := binary.Read(r, binary.LittleEndian, &lAux); err != nil {
		return nil, fmt.Errorf("csi: reading l_aux: %w", err)
	}

	if lAux > 0 {
		auxData := make([]byte, lAux)
		if _, err := io.ReadFull(r, auxData); err != nil {
			return nil, fmt.Errorf("csi: reading aux data: %w", err)
		}
		if err := idx.parseTabixAux(auxData); err != nil {
			return nil, fmt.Errorf("csi: parsing aux: %w", err)
		}
	}

	var nRef int32
	if err := binary.Read(r, binary.LittleEndian, &nRef); err != nil {
		return nil, fmt.Errorf("csi: reading n_ref: %w", err)
	}

	idx.refs = make([]csiRefIndex, nRef)
	for i := int32(0); i < nRef; i++ {
		if err := idx.readCSIRef(r, &idx.refs[i]); err != nil {
			return nil, fmt.Errorf("csi: ref %d: %w", i, err)
		}
	}

	return idx, nil
}

func (idx *CSIIndex) parseTabixAux(data []byte) error {
	if len(data) < 28 {
		return nil
	}

	idx.Format = int32(binary.LittleEndian.Uint32(data[0:4]))
	idx.ColSeq = int32(binary.LittleEndian.Uint32(data[4:8]))
	idx.ColBeg = int32(binary.LittleEndian.Uint32(data[8:12]))
	idx.ColEnd = int32(binary.LittleEndian.Uint32(data[12:16]))
	idx.Meta = int32(binary.LittleEndian.Uint32(data[16:20]))
	idx.Skip = int32(binary.LittleEndian.Uint32(data[20:24]))
	namesLen := int32(binary.LittleEndian.Uint32(data[24:28]))

	if idx.Format&0x10000 != 0 {
		idx.ZeroBased = true
		idx.Format &= 0xFFFF
	}

	if int(28+namesLen) <= len(data) {
		idx.Names = splitNulTerminated(data[28 : 28+namesLen])
	}

	return nil
}

func (idx *CSIIndex) readCSIRef(r io.Reader, ref *csiRefIndex) error {
	var nBins int32
	if err := binary.Read(r, binary.LittleEndian, &nBins); err != nil {
		return fmt.Errorf("reading n_bins: %w", err)
	}

	ref.bins = make(map[uint32]*csiBin, nBins)
	for i := int32(0); i < nBins; i++ {
		var binNum uint32
		if err := binary.Read(r, binary.LittleEndian, &binNum); err != nil {
			return fmt.Errorf("reading bin number: %w", err)
		}

		var loffset uint64
		if err := binary.Read(r, binary.LittleEndian, &loffset); err != nil {
			return fmt.Errorf("reading loffset: %w", err)
		}

		var nChunks int32
		if err := binary.Read(r, binary.LittleEndian, &nChunks); err != nil {
			return fmt.Errorf("reading n_chunks: %w", err)
		}

		chunks := make([]Chunk, nChunks)
		for j := int32(0); j < nChunks; j++ {
			var begin, end uint64
			if err := binary.Read(r, binary.LittleEndian, &begin); err != nil {
				return fmt.Errorf("reading chunk begin: %w", err)
			}
			if err := binary.Read(r, binary.LittleEndian, &end); err != nil {
				return fmt.Errorf("reading chunk end: %w", err)
			}
			chunks[j] = Chunk{
				Begin: bgzf.VirtualOffset(begin),
				End:   bgzf.VirtualOffset(end),
			}
		}

		ref.bins[binNum] = &csiBin{
			bin:     binNum,
			loffset: bgzf.VirtualOffset(loffset),
			chunks:  chunks,
		}
	}

	return nil
}

// csiReg2bins returns the list of bins that may overlap the region [beg, end)
// using the CSI variable binning scheme with the given minShift and depth.
func csiReg2bins(beg, end int, minShift, depth int32) []uint32 {
	end--
	var bins []uint32
	bins = append(bins, 0)
	for level := int32(1); level <= depth; level++ {
		shift := minShift + 3*(depth-level)
		offset := ((1 << (3 * level)) - 1) / 7
		lo := offset + (beg >> shift)
		hi := offset + (end >> shift)
		for k := lo; k <= hi; k++ {
			bins = append(bins, uint32(k))
		}
	}
	return bins
}

// Query returns the sorted, merged list of chunks that contain records
// overlapping the 0-based half-open region [start, end) on the given
// reference sequence (by index).
func (idx *CSIIndex) Query(refID int, start, end int) []Chunk {
	if refID < 0 || refID >= len(idx.refs) {
		return nil
	}
	ref := &idx.refs[refID]

	var minOffset bgzf.VirtualOffset
	leafShift := int(idx.minShift)
	leafOffset := ((1 << (3 * int(idx.depth))) - 1) / 7
	startBin := uint32(leafOffset + start>>leafShift)
	if b, ok := ref.bins[startBin]; ok {
		minOffset = b.loffset
	}

	bins := csiReg2bins(start, end, idx.minShift, idx.depth)
	var chunks []Chunk
	for _, binNum := range bins {
		if b, ok := ref.bins[binNum]; ok {
			for _, c := range b.chunks {
				if c.End <= minOffset {
					continue
				}
				chunks = append(chunks, c)
			}
		}
	}

	if len(chunks) == 0 {
		return nil
	}

	sort.Slice(chunks, func(i, j int) bool {
		return chunks[i].Begin < chunks[j].Begin
	})

	merged := []Chunk{chunks[0]}
	for i := 1; i < len(chunks); i++ {
		last := &merged[len(merged)-1]
		if chunks[i].Begin <= last.End {
			if chunks[i].End > last.End {
				last.End = chunks[i].End
			}
		} else {
			merged = append(merged, chunks[i])
		}
	}

	return merged
}

// RefID returns the reference sequence index for the given name.
func (idx *CSIIndex) RefID(name string) int {
	for i, n := range idx.Names {
		if n == name {
			return i
		}
	}
	return -1
}
