package bgzf

import (
	"container/list"
	"encoding/binary"
	"fmt"
	"io"
	"os"
	"sort"
)

const (
	// DefaultCacheSize is the default number of decompressed blocks to keep
	// in the LRU cache.
	DefaultCacheSize = 64
)

// cachedBlock holds a decompressed BGZF block and its compressed file offset.
type cachedBlock struct {
	offset int64  // compressed block offset in the file
	data   []byte // decompressed block data (owned by the cache)
	bsize  int    // total compressed block size (for advancing past the block)
}

// blockCache is an LRU cache of decompressed BGZF blocks keyed by compressed
// file offset.
type blockCache struct {
	maxSize int
	ll      *list.List
	items   map[int64]*list.Element
}

func newBlockCache(maxSize int) *blockCache {
	return &blockCache{
		maxSize: maxSize,
		ll:      list.New(),
		items:   make(map[int64]*list.Element, maxSize),
	}
}

// get retrieves a cached block by its compressed file offset.
// Returns nil if not cached.
func (c *blockCache) get(offset int64) *cachedBlock {
	if el, ok := c.items[offset]; ok {
		c.ll.MoveToFront(el)
		return el.Value.(*cachedBlock)
	}
	return nil
}

// put adds a block to the cache. If the cache is full, the least recently
// used block is evicted.
func (c *blockCache) put(b *cachedBlock) {
	if el, ok := c.items[b.offset]; ok {
		c.ll.MoveToFront(el)
		el.Value = b
		return
	}
	if c.ll.Len() >= c.maxSize {
		// Evict LRU
		back := c.ll.Back()
		if back != nil {
			evicted := c.ll.Remove(back).(*cachedBlock)
			delete(c.items, evicted.offset)
		}
	}
	el := c.ll.PushFront(b)
	c.items[b.offset] = el
}

// gziEntry is one entry from a .gzi index: the compressed and uncompressed
// offsets of a BGZF block boundary.
type gziEntry struct {
	compressedOffset   int64
	uncompressedOffset int64
}

// GZIndex is a parsed .gzi index that maps uncompressed byte positions to
// compressed block offsets. It enables Seek(uncompressedOffset) on an
// IndexedReader.
type GZIndex struct {
	entries []gziEntry
}

// LoadGZIndex reads a .gzi index file. The format is: uint64 count, then
// count pairs of (uint64 compressedOffset, uint64 uncompressedOffset).
func LoadGZIndex(filename string) (*GZIndex, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("bgzf: opening gzi: %w", err)
	}
	defer f.Close()

	var count uint64
	if err := binary.Read(f, binary.LittleEndian, &count); err != nil {
		return nil, fmt.Errorf("bgzf: reading gzi count: %w", err)
	}

	// The first block (offset 0,0) is implicit — the gzi stores only
	// subsequent block boundaries.
	entries := make([]gziEntry, 1, count+1)
	entries[0] = gziEntry{0, 0}

	for i := uint64(0); i < count; i++ {
		var cOff, uOff uint64
		if err := binary.Read(f, binary.LittleEndian, &cOff); err != nil {
			return nil, fmt.Errorf("bgzf: reading gzi entry %d: %w", i, err)
		}
		if err := binary.Read(f, binary.LittleEndian, &uOff); err != nil {
			return nil, fmt.Errorf("bgzf: reading gzi entry %d: %w", i, err)
		}
		entries = append(entries, gziEntry{
			compressedOffset:   int64(cOff),
			uncompressedOffset: int64(uOff),
		})
	}

	return &GZIndex{entries: entries}, nil
}

// lookup finds the block containing the given uncompressed offset.
// Returns the compressed block offset and the offset within the
// uncompressed block.
func (idx *GZIndex) lookup(uncompressedOffset int64) (compressedOffset int64, withinBlock int64) {
	// Binary search for the last entry whose uncompressedOffset <= target.
	i := sort.Search(len(idx.entries), func(i int) bool {
		return idx.entries[i].uncompressedOffset > uncompressedOffset
	}) - 1
	if i < 0 {
		i = 0
	}
	e := idx.entries[i]
	return e.compressedOffset, uncompressedOffset - e.uncompressedOffset
}

// IndexedReader reads BGZF-compressed data with random access by virtual
// offset. It wraps an io.ReadSeeker and maintains an LRU cache of
// decompressed blocks to avoid redundant decompression when reading
// nearby regions.
//
// It implements io.Reader and io.ByteReader. If a GZIndex is loaded via
// SetGZIndex, it also supports Seek(offset, whence) by uncompressed position.
type IndexedReader struct {
	rs    io.ReadSeeker
	cache *blockCache
	gzi   *GZIndex // optional; enables Seek by uncompressed offset

	// Current reading state
	block       *cachedBlock // current block (from cache)
	pos         int          // read position within block.data
	blockOffset int64        // compressed offset of current block
	uPos        int64        // current uncompressed position (tracked only when gzi is set)
	err         error        // sticky error
}

// NewIndexedReader creates an IndexedReader with the default cache size.
func NewIndexedReader(rs io.ReadSeeker) *IndexedReader {
	return NewIndexedReaderSize(rs, DefaultCacheSize)
}

// NewIndexedReaderSize creates an IndexedReader with the specified cache size
// (number of decompressed blocks to keep).
func NewIndexedReaderSize(rs io.ReadSeeker, cacheSize int) *IndexedReader {
	if cacheSize < 1 {
		cacheSize = 1
	}
	return &IndexedReader{
		rs:    rs,
		cache: newBlockCache(cacheSize),
	}
}

// OpenIndexedReader opens a BGZF file by name and returns an IndexedReader.
// If a .gzi index file exists at filename.gzi, it is loaded automatically,
// enabling Seek by uncompressed offset.
func OpenIndexedReader(filename string) (*IndexedReader, error) {
	return OpenIndexedReaderSize(filename, DefaultCacheSize)
}

// OpenIndexedReaderSize opens a BGZF file by name with the specified cache size.
// If a .gzi index file exists at filename.gzi, it is loaded automatically.
func OpenIndexedReaderSize(filename string, cacheSize int) (*IndexedReader, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}

	r := NewIndexedReaderSize(f, cacheSize)

	// Auto-load .gzi index if it exists.
	gziPath := filename + ".gzi"
	if _, err := os.Stat(gziPath); err == nil {
		idx, err := LoadGZIndex(gziPath)
		if err != nil {
			f.Close()
			return nil, fmt.Errorf("bgzf: loading gzi index: %w", err)
		}
		r.gzi = idx
	}

	return r, nil
}

// setGZIndex attaches a .gzi index, enabling Seek by uncompressed offset.
func (r *IndexedReader) setGZIndex(idx *GZIndex) {
	r.gzi = idx
}

// Seek positions the reader at the given uncompressed byte offset.
// Requires a GZIndex to be loaded via SetGZIndex. whence follows the
// standard io.Seek* constants (only io.SeekStart is supported).
func (r *IndexedReader) Seek(offset int64, whence int) (int64, error) {
	if r.gzi == nil {
		return 0, fmt.Errorf("bgzf: Seek requires a GZIndex (call SetGZIndex)")
	}
	if whence != io.SeekStart {
		return 0, fmt.Errorf("bgzf: only SeekStart is supported")
	}
	if offset < 0 {
		return 0, fmt.Errorf("bgzf: negative seek offset")
	}

	compOff, withinBlock := r.gzi.lookup(offset)
	vo := NewVirtualOffset(compOff, uint16(withinBlock))
	if err := r.SeekToVirtualOffset(vo); err != nil {
		return 0, err
	}
	r.uPos = offset
	return offset, nil
}

// SeekToVirtualOffset positions the reader at the given virtual offset.
// Subsequent reads will start from the uncompressed position within the
// block identified by the virtual offset.
func (r *IndexedReader) SeekToVirtualOffset(vo VirtualOffset) error {
	r.err = nil // clear sticky error on seek

	blockOff := vo.BlockOffset()
	withinBlock := int(vo.WithinBlock())

	if err := r.loadBlock(blockOff); err != nil {
		r.err = err
		return err
	}

	if withinBlock > len(r.block.data) {
		err := fmt.Errorf("bgzf: within-block offset %d exceeds block size %d", withinBlock, len(r.block.data))
		r.err = err
		return err
	}

	r.pos = withinBlock
	return nil
}

// VirtualTell returns the virtual offset of the next byte that will be read.
func (r *IndexedReader) VirtualTell() VirtualOffset {
	return NewVirtualOffset(r.blockOffset, uint16(r.pos))
}

// Read implements io.Reader. After seeking, reads proceed sequentially
// through blocks.
func (r *IndexedReader) Read(p []byte) (int, error) {
	if r.err != nil {
		return 0, r.err
	}

	total := 0
	for len(p) > 0 {
		if r.block == nil || r.pos >= len(r.block.data) {
			if err := r.advance(); err != nil {
				r.err = err
				if total > 0 {
					return total, nil
				}
				return 0, err
			}
		}

		n := copy(p, r.block.data[r.pos:])
		r.pos += n
		r.uPos += int64(n)
		p = p[n:]
		total += n
	}
	return total, nil
}

// ReadByte implements io.ByteReader.
func (r *IndexedReader) ReadByte() (byte, error) {
	if r.err != nil {
		return 0, r.err
	}
	if r.block == nil || r.pos >= len(r.block.data) {
		if err := r.advance(); err != nil {
			r.err = err
			return 0, err
		}
	}
	b := r.block.data[r.pos]
	r.pos++
	r.uPos++
	return b, nil
}

// loadBlock loads the decompressed block at the given compressed file offset,
// using the cache if available.
func (r *IndexedReader) loadBlock(offset int64) error {
	// Check cache first.
	if b := r.cache.get(offset); b != nil {
		r.block = b
		r.blockOffset = offset
		return nil
	}

	// Cache miss: seek and decompress.
	if _, err := r.rs.Seek(offset, io.SeekStart); err != nil {
		return fmt.Errorf("bgzf: seek to %d: %w", offset, err)
	}

	hdr, err := readBlockHeader(r.rs)
	if err != nil {
		return err
	}

	blockSize := int(hdr.bsize) + 1

	data, err := readBlockData(r.rs, hdr.bsize, nil)
	if err != nil {
		return err
	}

	b := &cachedBlock{
		offset: offset,
		data:   data,
		bsize:  blockSize,
	}

	// Make a copy for the cache so the data is owned by the cache entry
	// and won't be overwritten.
	cached := &cachedBlock{
		offset: offset,
		data:   make([]byte, len(data)),
		bsize:  blockSize,
	}
	copy(cached.data, data)
	r.cache.put(cached)

	r.block = b
	r.blockOffset = offset
	return nil
}

// advance moves to the next sequential block after the current one.
func (r *IndexedReader) advance() error {
	if r.block == nil {
		return io.EOF
	}

	nextOffset := r.blockOffset + int64(r.block.bsize)
	if err := r.loadBlock(nextOffset); err != nil {
		return err
	}

	// EOF block (empty data).
	if len(r.block.data) == 0 {
		return io.EOF
	}

	r.pos = 0
	return nil
}
