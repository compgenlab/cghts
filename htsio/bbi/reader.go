package bbi

import (
	"bytes"
	"compress/zlib"
	"encoding/binary"
	"fmt"
	"io"
	"iter"
	"math"
	"os"
	"strconv"
	"strings"
)

// Magic numbers identifying the two BBI file kinds and the two internal trees.
// The magic is read in the file's native byte order — a BBI file may be written
// little- or big-endian, and the byte order is detected from the magic itself.
const (
	magicBigWig  = 0x888FFC26
	magicBigBed  = 0x8789F2EB
	magicChromBP = 0x78CA8C91 // chromosome B+ tree
	magicCIRTree = 0x2468ACE0 // R-tree data index
)

// Kind is the BBI file kind.
type Kind int

const (
	BigWig Kind = iota
	BigBed
)

// header is the fixed 64-byte BBI common header.
type header struct {
	version           uint16
	zoomLevels        uint16
	chromTreeOffset   uint64
	fullDataOffset    uint64
	fullIndexOffset   uint64
	fieldCount        uint16
	definedFieldCount uint16
	autoSQLOffset     uint64
	totalSummaryOff   uint64
	uncompressBufSize uint32
}

// Record is one item from a BBI query. For bigWig, Value holds the score and
// Line is empty; for bigBed, Line holds a BED-like row "chrom\tstart\tend\t<rest>"
// (start 0-based, matching UCSC bigBed) and Value is unused.
type Record struct {
	Chrom string
	Start int // 0-based
	End   int // 0-based, exclusive
	Value float64
	Line  string
}

// Reader reads a local UCSC BBI file (bigWig or bigBed) with random access by
// genomic region. It mirrors the surface of htsio/tabix.Reader (HasRef/RefNames/
// Query/Close) so annotators can treat the two interchangeably. Zoom-level
// summaries are never read — only the base-resolution data — so queried values
// are exact.
type Reader struct {
	f    *os.File
	bo   binary.ByteOrder
	kind Kind
	hdr  header

	names    []string          // reference names, in id order
	nameToID map[string]uint32 // reference name -> chrom id
	idToName map[uint32]string
}

// Open opens a BBI file (.bw/.bb) and parses its header and chromosome index.
func Open(filename string) (*Reader, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	r := &Reader{f: f, nameToID: map[string]uint32{}, idToName: map[uint32]string{}}
	if err := r.readHeader(); err != nil {
		f.Close()
		return nil, fmt.Errorf("bbi: %s: %w", filename, err)
	}
	if err := r.readChromTree(); err != nil {
		f.Close()
		return nil, fmt.Errorf("bbi: %s: %w", filename, err)
	}
	return r, nil
}

// Close releases the underlying file.
func (r *Reader) Close() error { return r.f.Close() }

// Kind reports whether this is a bigWig or bigBed file.
func (r *Reader) Kind() Kind { return r.kind }

// HasRef reports whether the file contains the given reference name.
func (r *Reader) HasRef(ref string) bool { _, ok := r.nameToID[ref]; return ok }

// RefNames returns the reference names present in the file, in id order — the
// contig list used to build a contig-name converter.
func (r *Reader) RefNames() []string { return r.names }

// readAt reads exactly n bytes at the given file offset.
func (r *Reader) readAt(off int64, n int) ([]byte, error) {
	b := make([]byte, n)
	if _, err := r.f.ReadAt(b, off); err != nil {
		return nil, err
	}
	return b, nil
}

func (r *Reader) readHeader() error {
	b, err := r.readAt(0, 64)
	if err != nil {
		return err
	}
	// Detect byte order from the magic, then re-read the magic in that order.
	switch {
	case binary.LittleEndian.Uint32(b) == magicBigWig || binary.LittleEndian.Uint32(b) == magicBigBed:
		r.bo = binary.LittleEndian
	case binary.BigEndian.Uint32(b) == magicBigWig || binary.BigEndian.Uint32(b) == magicBigBed:
		r.bo = binary.BigEndian
	default:
		return fmt.Errorf("not a bigWig/bigBed file (bad magic)")
	}
	switch r.bo.Uint32(b) {
	case magicBigWig:
		r.kind = BigWig
	case magicBigBed:
		r.kind = BigBed
	}
	r.hdr = header{
		version:           r.bo.Uint16(b[4:]),
		zoomLevels:        r.bo.Uint16(b[6:]),
		chromTreeOffset:   r.bo.Uint64(b[8:]),
		fullDataOffset:    r.bo.Uint64(b[16:]),
		fullIndexOffset:   r.bo.Uint64(b[24:]),
		fieldCount:        r.bo.Uint16(b[32:]),
		definedFieldCount: r.bo.Uint16(b[34:]),
		autoSQLOffset:     r.bo.Uint64(b[36:]),
		totalSummaryOff:   r.bo.Uint64(b[44:]),
		uncompressBufSize: r.bo.Uint32(b[52:]),
	}
	return nil
}

// readChromTree walks the chromosome B+ tree, populating the name<->id maps.
func (r *Reader) readChromTree() error {
	off := int64(r.hdr.chromTreeOffset)
	th, err := r.readAt(off, 32)
	if err != nil {
		return err
	}
	if r.bo.Uint32(th) != magicChromBP {
		return fmt.Errorf("bad chromosome B+ tree magic")
	}
	keySize := int(r.bo.Uint32(th[8:]))
	// root node starts right after the 32-byte tree header
	if err := r.walkChromNode(off+32, keySize); err != nil {
		return err
	}
	// names in id order
	r.names = make([]string, 0, len(r.idToName))
	max := uint32(0)
	for id := range r.idToName {
		if id > max {
			max = id
		}
	}
	tmp := make([]string, max+1)
	for id, n := range r.idToName {
		tmp[id] = n
	}
	for _, n := range tmp {
		if n != "" {
			r.names = append(r.names, n)
		}
	}
	return nil
}

func (r *Reader) walkChromNode(off int64, keySize int) error {
	nh, err := r.readAt(off, 4)
	if err != nil {
		return err
	}
	isLeaf := nh[0] != 0
	count := int(r.bo.Uint16(nh[2:]))
	itemSize := keySize + 8 // leaf: key+chromId+chromSize; internal: key+childOffset
	body, err := r.readAt(off+4, count*itemSize)
	if err != nil {
		return err
	}
	for i := 0; i < count; i++ {
		item := body[i*itemSize : (i+1)*itemSize]
		key := strings.TrimRight(string(item[:keySize]), "\x00")
		rest := item[keySize:]
		if isLeaf {
			id := r.bo.Uint32(rest)
			r.nameToID[key] = id
			r.idToName[id] = key
		} else {
			child := int64(r.bo.Uint64(rest))
			if err := r.walkChromNode(child, keySize); err != nil {
				return err
			}
		}
	}
	return nil
}

// Query returns an iterator over records overlapping the 0-based half-open
// region [start, end) on the given reference.
func (r *Reader) Query(ref string, start, end int) (iter.Seq2[*Record, error], error) {
	id, ok := r.nameToID[ref]
	if !ok {
		return nil, fmt.Errorf("bbi: unknown reference %q", ref)
	}
	if start < 0 {
		start = 0
	}
	if end < start {
		end = start
	}
	blocks, err := r.findBlocks(id, uint32(start), uint32(end))
	if err != nil {
		return nil, err
	}
	return func(yield func(*Record, error) bool) {
		for _, blk := range blocks {
			raw, err := r.readBlock(blk)
			if err != nil {
				yield(nil, err)
				return
			}
			ok := r.emitBlock(raw, id, uint32(start), uint32(end), yield)
			if !ok {
				return
			}
		}
	}, nil
}

type block struct{ offset, size uint64 }

// findBlocks walks the R-tree, collecting data blocks overlapping the query.
func (r *Reader) findBlocks(chromID, start, end uint32) ([]block, error) {
	off := int64(r.hdr.fullIndexOffset)
	rh, err := r.readAt(off, 48)
	if err != nil {
		return nil, err
	}
	if r.bo.Uint32(rh) != magicCIRTree {
		return nil, fmt.Errorf("bad R-tree magic")
	}
	var out []block
	if err := r.walkRTree(off+48, chromID, start, end, &out); err != nil {
		return nil, err
	}
	return out, nil
}

// cirLess reports whether coordinate (c1,b1) < (c2,b2).
func cirLess(c1, b1, c2, b2 uint32) bool {
	if c1 != c2 {
		return c1 < c2
	}
	return b1 < b2
}

// overlaps reports whether span [(sC,sB),(eC,eB)) overlaps query chrom qC [qS,qE).
func overlaps(qC, qS, qE, sC, sB, eC, eB uint32) bool {
	return cirLess(qC, qS, eC, eB) && cirLess(sC, sB, qC, qE)
}

func (r *Reader) walkRTree(off int64, qC, qS, qE uint32, out *[]block) error {
	nh, err := r.readAt(off, 4)
	if err != nil {
		return err
	}
	isLeaf := nh[0] != 0
	count := int(r.bo.Uint16(nh[2:]))
	itemSize := 24
	if isLeaf {
		itemSize = 32
	}
	body, err := r.readAt(off+4, count*itemSize)
	if err != nil {
		return err
	}
	for i := 0; i < count; i++ {
		it := body[i*itemSize : (i+1)*itemSize]
		sC := r.bo.Uint32(it[0:])
		sB := r.bo.Uint32(it[4:])
		eC := r.bo.Uint32(it[8:])
		eB := r.bo.Uint32(it[12:])
		if !overlaps(qC, qS, qE, sC, sB, eC, eB) {
			continue
		}
		if isLeaf {
			*out = append(*out, block{offset: r.bo.Uint64(it[16:]), size: r.bo.Uint64(it[24:])})
		} else {
			if err := r.walkRTree(int64(r.bo.Uint64(it[16:])), qC, qS, qE, out); err != nil {
				return err
			}
		}
	}
	return nil
}

// readBlock reads and (if configured) decompresses a data block.
func (r *Reader) readBlock(blk block) ([]byte, error) {
	raw, err := r.readAt(int64(blk.offset), int(blk.size))
	if err != nil {
		return nil, err
	}
	if r.hdr.uncompressBufSize == 0 {
		return raw, nil // stored uncompressed
	}
	zr, err := zlib.NewReader(bytes.NewReader(raw))
	if err != nil {
		return nil, fmt.Errorf("bbi: inflate block: %w", err)
	}
	defer zr.Close()
	out, err := io.ReadAll(zr)
	if err != nil {
		return nil, fmt.Errorf("bbi: inflate block: %w", err)
	}
	return out, nil
}

// emitBlock parses a decompressed block and yields records overlapping the query.
// It returns false if the consumer stopped iteration.
func (r *Reader) emitBlock(buf []byte, qC, qS, qE uint32, yield func(*Record, error) bool) bool {
	if r.kind == BigWig {
		return r.emitBigWig(buf, qC, qS, qE, yield)
	}
	return r.emitBigBed(buf, qC, qS, qE, yield)
}

func (r *Reader) emitBigWig(buf []byte, qC, qS, qE uint32, yield func(*Record, error) bool) bool {
	for len(buf) >= 24 {
		chromID := r.bo.Uint32(buf[0:])
		chromStart := r.bo.Uint32(buf[4:])
		itemStep := r.bo.Uint32(buf[12:])
		itemSpan := r.bo.Uint32(buf[16:])
		sectType := buf[20]
		itemCount := int(r.bo.Uint16(buf[22:]))
		buf = buf[24:]

		var itemSize int
		switch sectType {
		case 1: // bedGraph: start,end,val
			itemSize = 12
		case 2: // varStep: start,val
			itemSize = 8
		case 3: // fixedStep: val
			itemSize = 4
		default:
			return true // unknown section; skip the rest of this block
		}
		if len(buf) < itemCount*itemSize {
			return true
		}
		for i := 0; i < itemCount; i++ {
			it := buf[i*itemSize : (i+1)*itemSize]
			var s, e uint32
			var val float32
			switch sectType {
			case 1:
				s = r.bo.Uint32(it[0:])
				e = r.bo.Uint32(it[4:])
				val = math.Float32frombits(r.bo.Uint32(it[8:]))
			case 2:
				s = r.bo.Uint32(it[0:])
				e = s + itemSpan
				val = math.Float32frombits(r.bo.Uint32(it[4:]))
			case 3:
				s = chromStart + uint32(i)*itemStep
				e = s + itemSpan
				val = math.Float32frombits(r.bo.Uint32(it[0:]))
			}
			if chromID != qC || e <= qS || s >= qE {
				continue
			}
			rec := &Record{Chrom: r.idToName[chromID], Start: int(s), End: int(e), Value: float64(val)}
			if !yield(rec, nil) {
				return false
			}
		}
		buf = buf[itemCount*itemSize:]
	}
	return true
}

func (r *Reader) emitBigBed(buf []byte, qC, qS, qE uint32, yield func(*Record, error) bool) bool {
	for len(buf) >= 12 {
		chromID := r.bo.Uint32(buf[0:])
		s := r.bo.Uint32(buf[4:])
		e := r.bo.Uint32(buf[8:])
		buf = buf[12:]
		// rest = null-terminated tab-separated remaining fields
		nul := bytes.IndexByte(buf, 0)
		if nul < 0 {
			nul = len(buf)
		}
		rest := string(buf[:nul])
		if nul < len(buf) {
			buf = buf[nul+1:]
		} else {
			buf = buf[nul:]
		}
		if chromID != qC || e <= qS || s >= qE {
			continue
		}
		name := r.idToName[chromID]
		line := name + "\t" + strconv.FormatUint(uint64(s), 10) + "\t" + strconv.FormatUint(uint64(e), 10)
		if rest != "" {
			line += "\t" + rest
		}
		rec := &Record{Chrom: name, Start: int(s), End: int(e), Line: line}
		if !yield(rec, nil) {
			return false
		}
	}
	return true
}
