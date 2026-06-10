package seqio

import (
	"bufio"
	"container/list"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
	"sync"

	"github.com/compgen-io/cgkit/htsio/bgzf"
)

const (
	faiChunkSize    = 10 * 1024 * 1024 // 10MB of bases per chunk
	faiCacheMaxSize = 1 << 30          // 1GB max cache
)

// FaiEntry represents one sequence in a .fai index file.
type FaiEntry struct {
	Name      string
	Length    int
	Offset    int64 // byte offset of first base (uncompressed offset for bgzip)
	LineBases int   // number of bases per text line
	LineWidth int   // number of bytes per text line (including newline)
}

// IndexedFastaReader provides random access to an indexed FASTA file.
// Supports both plain text and bgzip-compressed FASTA.
// Requires a .fai index file. For bgzip, also uses .gzi if available.
//
// Sequences are loaded in 10MB chunks with an LRU cache (max 1GB).
type IndexedFastaReader struct {
	path  string
	fai   map[string]*FaiEntry
	names []string // ordered sequence names from .fai

	// Exactly one of these is set.
	file *os.File            // plain FASTA
	bgzr *bgzf.IndexedReader // bgzip FASTA

	cache *faiChunkCache
	mu    sync.Mutex
}

// NewIndexedFastaReader opens a FASTA file with its .fai index for random access.
// The .fai file must exist at fastaPath+".fai". For bgzip-compressed FASTA,
// the .gzi index is loaded automatically if present.
func NewIndexedFastaReader(fastaPath string) (*IndexedFastaReader, error) {
	// Parse .fai index.
	faiPath := fastaPath + ".fai"
	fai, names, err := parseFaiIndex(faiPath)
	if err != nil {
		return nil, fmt.Errorf("loading FASTA index: %w", err)
	}

	r := &IndexedFastaReader{
		path:  fastaPath,
		fai:   fai,
		names: names,
		cache: newFaiChunkCache(faiCacheMaxSize),
	}

	// Detect format: bgzip or plain.
	f, err := os.Open(fastaPath)
	if err != nil {
		return nil, err
	}

	magic := make([]byte, 2)
	n, _ := f.Read(magic)
	if n >= 2 && magic[0] == 0x1f && magic[1] == 0x8b {
		// BGZF/gzip magic — open as bgzip.
		f.Close()
		bgzr, err := bgzf.OpenIndexedReader(fastaPath)
		if err != nil {
			return nil, fmt.Errorf("opening bgzip FASTA: %w", err)
		}
		r.bgzr = bgzr
	} else {
		// Plain FASTA — keep file open for ReadAt.
		r.file = f
	}

	return r, nil
}

// Names returns the ordered sequence names from the .fai index.
func (r *IndexedFastaReader) Names() []string {
	return r.names
}

// SequenceLength returns the length of the named sequence, or false if not found.
func (r *IndexedFastaReader) SequenceLength(name string) (int, bool) {
	entry, ok := r.fai[name]
	if !ok {
		return 0, false
	}
	return entry.Length, true
}

// GetSequenceRange returns reference bases for [start, end) (0-based) of the
// named sequence. The returned slice is uppercased.
func (r *IndexedFastaReader) GetSequenceRange(name string, start, end int) ([]byte, error) {
	r.mu.Lock()
	defer r.mu.Unlock()

	entry, ok := r.fai[name]
	if !ok {
		return nil, fmt.Errorf("sequence %q not found in %s.fai", name, r.path)
	}

	if start < 0 {
		start = 0
	}
	if end > entry.Length {
		end = entry.Length
	}
	if start >= end {
		return nil, nil
	}

	// Determine which chunks cover [start, end).
	firstChunk := start / faiChunkSize
	lastChunk := (end - 1) / faiChunkSize

	result := make([]byte, 0, end-start)

	for ci := firstChunk; ci <= lastChunk; ci++ {
		chunk, err := r.loadChunk(name, ci, entry)
		if err != nil {
			return nil, err
		}

		chunkStart := ci * faiChunkSize
		lo := start - chunkStart
		if lo < 0 {
			lo = 0
		}
		hi := end - chunkStart
		if hi > len(chunk) {
			hi = len(chunk)
		}
		result = append(result, chunk[lo:hi]...)
	}

	return result, nil
}

// GetSequence returns the full sequence for the given name, uppercased.
// Prefer GetSequenceRange when you know the range you need.
func (r *IndexedFastaReader) GetSequence(name string) ([]byte, error) {
	entry, ok := r.fai[name]
	if !ok {
		return nil, fmt.Errorf("sequence %q not found in %s.fai", name, r.path)
	}
	return r.GetSequenceRange(name, 0, entry.Length)
}

// Close releases resources.
func (r *IndexedFastaReader) Close() error {
	if r.file != nil {
		return r.file.Close()
	}
	// bgzf.IndexedReader doesn't expose Close, but the underlying file
	// was opened by OpenIndexedReader which keeps the *os.File internally.
	return nil
}

// loadChunk loads a single chunk, using the cache. Must be called with mu held.
func (r *IndexedFastaReader) loadChunk(name string, chunkIdx int, entry *FaiEntry) ([]byte, error) {
	key := faiCacheKey{name: name, chunkIdx: chunkIdx}

	if data := r.cache.get(key); data != nil {
		return data, nil
	}

	// Compute base range for this chunk.
	baseStart := chunkIdx * faiChunkSize
	baseEnd := baseStart + faiChunkSize
	if baseEnd > entry.Length {
		baseEnd = entry.Length
	}

	// Compute byte range in the (uncompressed) FASTA.
	startByte := entry.Offset + int64(baseStart/entry.LineBases)*int64(entry.LineWidth) + int64(baseStart%entry.LineBases)
	lastBase := baseEnd - 1
	endByte := entry.Offset + int64(lastBase/entry.LineBases)*int64(entry.LineWidth) + int64(lastBase%entry.LineBases) + 1

	byteLen := int(endByte - startByte)
	buf := make([]byte, byteLen)

	var n int
	var err error
	if r.file != nil {
		// Plain FASTA: use ReadAt for thread-safe random access.
		n, err = r.file.ReadAt(buf, startByte)
	} else {
		// Bgzip FASTA: seek + read.
		if _, err = r.bgzr.Seek(startByte, io.SeekStart); err != nil {
			return nil, fmt.Errorf("seeking in bgzip FASTA for %s chunk %d: %w", name, chunkIdx, err)
		}
		n, err = io.ReadFull(r.bgzr, buf)
	}
	if err != nil && n == 0 {
		return nil, fmt.Errorf("reading %s chunk %d: %w", name, chunkIdx, err)
	}
	buf = buf[:n]

	// Strip newlines and uppercase.
	bases := make([]byte, 0, baseEnd-baseStart)
	for _, b := range buf {
		if b != '\n' && b != '\r' {
			if b >= 'a' && b <= 'z' {
				b -= 32
			}
			bases = append(bases, b)
		}
	}

	r.cache.put(key, bases)
	return bases, nil
}

// parseFaiIndex reads a .fai file and returns entries + ordered names.
func parseFaiIndex(path string) (map[string]*FaiEntry, []string, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	defer f.Close()

	fai := make(map[string]*FaiEntry)
	var names []string
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := scanner.Text()
		if line == "" {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 5 {
			return nil, nil, fmt.Errorf("malformed .fai line: %s", line)
		}
		length, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, nil, fmt.Errorf("bad length in .fai: %s", fields[1])
		}
		offset, err := strconv.ParseInt(fields[2], 10, 64)
		if err != nil {
			return nil, nil, fmt.Errorf("bad offset in .fai: %s", fields[2])
		}
		lineBases, err := strconv.Atoi(fields[3])
		if err != nil {
			return nil, nil, fmt.Errorf("bad lineBases in .fai: %s", fields[3])
		}
		lineWidth, err := strconv.Atoi(fields[4])
		if err != nil {
			return nil, nil, fmt.Errorf("bad lineWidth in .fai: %s", fields[4])
		}
		// lineBases is used as a divisor when mapping a base offset to a byte
		// offset in loadChunk; a non-positive line geometry would divide by zero
		// or produce nonsensical offsets. Reject malformed entries up front.
		if length < 0 || lineBases <= 0 || lineWidth <= 0 {
			return nil, nil, fmt.Errorf("invalid .fai geometry for %s: length=%d lineBases=%d lineWidth=%d", fields[0], length, lineBases, lineWidth)
		}
		name := fields[0]
		fai[name] = &FaiEntry{
			Name:      name,
			Length:    length,
			Offset:    offset,
			LineBases: lineBases,
			LineWidth: lineWidth,
		}
		names = append(names, name)
	}
	if err := scanner.Err(); err != nil {
		return nil, nil, err
	}
	if len(fai) == 0 {
		return nil, nil, fmt.Errorf("empty .fai file: %s", path)
	}
	return fai, names, nil
}

// --- LRU chunk cache ---

type faiCacheKey struct {
	name     string
	chunkIdx int
}

type faiCacheValue struct {
	key  faiCacheKey
	data []byte
}

type faiChunkCache struct {
	entries   map[faiCacheKey]*list.Element
	lru       *list.List
	totalSize int
	maxSize   int
}

func newFaiChunkCache(maxSize int) *faiChunkCache {
	return &faiChunkCache{
		entries: make(map[faiCacheKey]*list.Element),
		lru:     list.New(),
		maxSize: maxSize,
	}
}

func (c *faiChunkCache) get(key faiCacheKey) []byte {
	if el, ok := c.entries[key]; ok {
		c.lru.MoveToFront(el)
		return el.Value.(*faiCacheValue).data
	}
	return nil
}

func (c *faiChunkCache) put(key faiCacheKey, data []byte) {
	if el, ok := c.entries[key]; ok {
		c.lru.MoveToFront(el)
		return
	}
	for c.totalSize+len(data) > c.maxSize && c.lru.Len() > 0 {
		back := c.lru.Back()
		cv := c.lru.Remove(back).(*faiCacheValue)
		delete(c.entries, cv.key)
		c.totalSize -= len(cv.data)
	}
	el := c.lru.PushFront(&faiCacheValue{key: key, data: data})
	c.entries[key] = el
	c.totalSize += len(data)
}
