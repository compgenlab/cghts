package tabix

import (
	"encoding/binary"
	"fmt"
	"io"
	"os"
	"sort"

	"github.com/compgen-io/cgkit/htsio/bgzf"
)

// Chunk represents a contiguous range of BGZF virtual offsets containing
// alignment records or text lines.
type Chunk struct {
	Begin bgzf.VirtualOffset
	End   bgzf.VirtualOffset
}

// binIndex holds the bins and linear index for one reference sequence.
type binIndex struct {
	bins      map[uint32][]Chunk   // bin number → chunks
	linearIdx []bgzf.VirtualOffset // one entry per 16kb window
}

// BinIndex is a parsed BAI or TBI index.
type BinIndex struct {
	refs []binIndex

	// TBI-specific header fields (zero-valued for BAI).
	Format    int32    // 0=generic, 1=SAM, 2=VCF
	ColSeq    int32    // column for sequence name (1-based)
	ColBeg    int32    // column for region start (1-based)
	ColEnd    int32    // column for region end (1-based); 0 means same as ColBeg
	Meta      int32    // comment character (e.g. '#'), or 0
	Skip      int32    // number of header lines to skip
	ZeroBased bool     // true if coordinates are 0-based
	Names     []string // reference sequence names (TBI only)
}

// LoadBAI reads a BAI index file.
func LoadBAI(filename string) (*BinIndex, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	var magic [4]byte
	if _, err := io.ReadFull(f, magic[:]); err != nil {
		return nil, fmt.Errorf("bai: reading magic: %w", err)
	}
	if magic != [4]byte{'B', 'A', 'I', 1} {
		return nil, fmt.Errorf("bai: invalid magic: %x", magic)
	}

	return readBinIndex(f)
}

// LoadTBI reads a TBI (tabix) index file.
func LoadTBI(filename string) (*BinIndex, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	r := bgzf.NewReader(f)

	var magic [4]byte
	if _, err := io.ReadFull(r, magic[:]); err != nil {
		return nil, fmt.Errorf("tbi: reading magic: %w", err)
	}
	if magic != [4]byte{'T', 'B', 'I', 1} {
		return nil, fmt.Errorf("tbi: invalid magic: %x", magic)
	}

	var nRef int32
	if err := binary.Read(r, binary.LittleEndian, &nRef); err != nil {
		return nil, fmt.Errorf("tbi: reading n_ref: %w", err)
	}

	idx := &BinIndex{}

	if err := binary.Read(r, binary.LittleEndian, &idx.Format); err != nil {
		return nil, fmt.Errorf("tbi: reading format: %w", err)
	}
	if err := binary.Read(r, binary.LittleEndian, &idx.ColSeq); err != nil {
		return nil, fmt.Errorf("tbi: reading col_seq: %w", err)
	}
	if err := binary.Read(r, binary.LittleEndian, &idx.ColBeg); err != nil {
		return nil, fmt.Errorf("tbi: reading col_beg: %w", err)
	}
	if err := binary.Read(r, binary.LittleEndian, &idx.ColEnd); err != nil {
		return nil, fmt.Errorf("tbi: reading col_end: %w", err)
	}
	if err := binary.Read(r, binary.LittleEndian, &idx.Meta); err != nil {
		return nil, fmt.Errorf("tbi: reading meta: %w", err)
	}
	if err := binary.Read(r, binary.LittleEndian, &idx.Skip); err != nil {
		return nil, fmt.Errorf("tbi: reading skip: %w", err)
	}

	if idx.Format&0x10000 != 0 {
		idx.ZeroBased = true
		idx.Format &= 0xFFFF
	}

	var namesLen int32
	if err := binary.Read(r, binary.LittleEndian, &namesLen); err != nil {
		return nil, fmt.Errorf("tbi: reading names length: %w", err)
	}
	namesData := make([]byte, namesLen)
	if _, err := io.ReadFull(r, namesData); err != nil {
		return nil, fmt.Errorf("tbi: reading names: %w", err)
	}

	idx.Names = splitNulTerminated(namesData)

	idx.refs = make([]binIndex, nRef)
	for i := int32(0); i < nRef; i++ {
		if err := readRefBins(r, &idx.refs[i]); err != nil {
			return nil, fmt.Errorf("tbi: ref %d: %w", i, err)
		}
	}

	return idx, nil
}

func readBinIndex(r io.Reader) (*BinIndex, error) {
	var nRef int32
	if err := binary.Read(r, binary.LittleEndian, &nRef); err != nil {
		return nil, fmt.Errorf("reading n_ref: %w", err)
	}

	idx := &BinIndex{
		refs: make([]binIndex, nRef),
	}

	for i := int32(0); i < nRef; i++ {
		if err := readRefBins(r, &idx.refs[i]); err != nil {
			return nil, fmt.Errorf("ref %d: %w", i, err)
		}
	}

	return idx, nil
}

func readRefBins(r io.Reader, ref *binIndex) error {
	var nBins int32
	if err := binary.Read(r, binary.LittleEndian, &nBins); err != nil {
		return fmt.Errorf("reading n_bins: %w", err)
	}

	ref.bins = make(map[uint32][]Chunk, nBins)
	for i := int32(0); i < nBins; i++ {
		var binNum uint32
		if err := binary.Read(r, binary.LittleEndian, &binNum); err != nil {
			return fmt.Errorf("reading bin number: %w", err)
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
		ref.bins[binNum] = chunks
	}

	var nIntervals int32
	if err := binary.Read(r, binary.LittleEndian, &nIntervals); err != nil {
		return fmt.Errorf("reading n_intervals: %w", err)
	}

	ref.linearIdx = make([]bgzf.VirtualOffset, nIntervals)
	for i := int32(0); i < nIntervals; i++ {
		var offset uint64
		if err := binary.Read(r, binary.LittleEndian, &offset); err != nil {
			return fmt.Errorf("reading linear index: %w", err)
		}
		ref.linearIdx[i] = bgzf.VirtualOffset(offset)
	}

	return nil
}

// reg2bins returns the list of bins that may overlap the region [beg, end)
// using the standard BAI/TBI 6-level binning scheme.
func reg2bins(beg, end int) []uint32 {
	end--
	var bins []uint32
	bins = append(bins, 0)
	for k := 1 + (beg >> 26); k <= 1+(end>>26); k++ {
		bins = append(bins, uint32(k))
	}
	for k := 9 + (beg >> 23); k <= 9+(end>>23); k++ {
		bins = append(bins, uint32(k))
	}
	for k := 73 + (beg >> 20); k <= 73+(end>>20); k++ {
		bins = append(bins, uint32(k))
	}
	for k := 585 + (beg >> 17); k <= 585+(end>>17); k++ {
		bins = append(bins, uint32(k))
	}
	for k := 4681 + (beg >> 14); k <= 4681+(end>>14); k++ {
		bins = append(bins, uint32(k))
	}
	return bins
}

// Query returns the sorted, merged list of chunks that contain records
// overlapping the 0-based half-open region [start, end) on the given
// reference sequence (by index).
func (idx *BinIndex) Query(refID int, start, end int) []Chunk {
	if refID < 0 || refID >= len(idx.refs) {
		return nil
	}
	ref := &idx.refs[refID]

	linearMinIdx := start >> 14
	var minOffset bgzf.VirtualOffset
	if linearMinIdx < len(ref.linearIdx) {
		minOffset = ref.linearIdx[linearMinIdx]
	}

	bins := reg2bins(start, end)
	var chunks []Chunk
	for _, b := range bins {
		if cs, ok := ref.bins[b]; ok {
			for _, c := range cs {
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
// Returns -1 if not found.
func (idx *BinIndex) RefID(name string) int {
	for i, n := range idx.Names {
		if n == name {
			return i
		}
	}
	return -1
}

// Reg2Bin calculates the BAM bin for a region [beg, end) using the BAM
// binning scheme (as defined in the SAM specification).
func Reg2Bin(beg, end int) uint16 {
	end--
	// Parentheses around the shifts are redundant in Go (>> binds tighter than
	// +), but are written explicitly here to match the SAM spec formula and to
	// avoid any C-style precedence confusion for future readers.
	if beg>>14 == end>>14 {
		return uint16((((1 << 15) - 1) / 7) + (beg >> 14))
	}
	if beg>>17 == end>>17 {
		return uint16((((1 << 12) - 1) / 7) + (beg >> 17))
	}
	if beg>>20 == end>>20 {
		return uint16((((1 << 9) - 1) / 7) + (beg >> 20))
	}
	if beg>>23 == end>>23 {
		return uint16((((1 << 6) - 1) / 7) + (beg >> 23))
	}
	if beg>>26 == end>>26 {
		return uint16((((1 << 3) - 1) / 7) + (beg >> 26))
	}
	return 0
}

func splitNulTerminated(data []byte) []string {
	var result []string
	start := 0
	for i, b := range data {
		if b == 0 {
			if i > start {
				result = append(result, string(data[start:i]))
			}
			start = i + 1
		}
	}
	return result
}
