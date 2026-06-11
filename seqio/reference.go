package seqio

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"os"
	"strings"
	"sync"
)

// ReferenceReader provides random access to named reference sequences.
// This is the interface for reference genomes — few long sequences accessed
// by name and position range. Implementations include indexed FASTA (local,
// plain or bgzip) and in-memory FASTA.
type ReferenceReader interface {
	// Names returns the ordered sequence names.
	Names() []string

	// SequenceLength returns the length of the named sequence.
	// Returns 0, false if the sequence is not found.
	SequenceLength(name string) (int, bool)

	// GetSequenceRange returns bases for [start, end) (0-based) of the
	// named sequence. The returned slice is uppercased. Out-of-range
	// coordinates are clamped to the sequence boundaries.
	GetSequenceRange(name string, start, end int) ([]byte, error)

	// GetSequence returns the full sequence for the given name, uppercased.
	// Prefer GetSequenceRange when you know the range you need.
	GetSequence(name string) ([]byte, error)

	// Close releases resources held by the reader.
	Close() error
}

// OpenReference opens a reference sequence, auto-detecting the format.
//
// Supported formats:
//   - HTTP/HTTPS URL: fetches .fai index, uses Range requests for chunks
//   - Indexed FASTA (.fa/.fasta with .fai index, plain or bgzip-compressed)
//   - Plain FASTA without index (loaded fully into memory as compressed chunks)
//
// For indexed FASTA and remote references, sequences are loaded in 10MB
// chunks with LRU caching (1GB max).
func OpenReference(path string) (ReferenceReader, error) {
	// Detect HTTP/HTTPS URLs.
	if strings.HasPrefix(path, "http://") || strings.HasPrefix(path, "https://") {
		return NewRemoteFastaReader(path)
	}

	if _, err := os.Stat(path); err != nil {
		return nil, fmt.Errorf("reference not found: %s", path)
	}

	// Try indexed mode first (.fai must exist).
	r, err := NewIndexedFastaReader(path)
	if err == nil {
		return r, nil
	}

	// Fall back to in-memory FASTA.
	return newInMemoryFasta(path)
}

// inMemoryFasta reads a FASTA file and stores sequences as compressed chunks.
// Chunks are decompressed on demand with an LRU cache, keeping memory usage
// bounded even for large genomes without a .fai index.
type inMemoryFasta struct {
	path    string
	names   []string
	lengths map[string]int
	chunks  map[faiCacheKey][]byte // gzip-compressed chunks
	cache   *faiChunkCache         // decompressed chunk LRU
	mu      sync.Mutex
}

func newInMemoryFasta(path string) (*inMemoryFasta, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("opening reference FASTA: %w", err)
	}
	defer f.Close()

	r := &inMemoryFasta{
		path:    path,
		lengths: make(map[string]int),
		chunks:  make(map[faiCacheKey][]byte),
		cache:   newFaiChunkCache(faiCacheMaxSize),
	}

	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 0, 1024*1024), 10*1024*1024)

	var name string
	var seq bytes.Buffer

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			if name != "" {
				if err := r.storeSequence(name, bytes.ToUpper(seq.Bytes())); err != nil {
					return nil, err
				}
			}
			fields := strings.Fields(line[1:])
			if len(fields) > 0 {
				name = fields[0]
			}
			r.names = append(r.names, name)
			seq.Reset()
		} else {
			seq.WriteString(strings.TrimSpace(line))
		}
	}
	if name != "" {
		if err := r.storeSequence(name, bytes.ToUpper(seq.Bytes())); err != nil {
			return nil, err
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return r, nil
}

// storeSequence breaks a sequence into chunks and compresses each one.
func (r *inMemoryFasta) storeSequence(name string, seq []byte) error {
	r.lengths[name] = len(seq)
	for i := 0; i < len(seq); i += faiChunkSize {
		end := i + faiChunkSize
		if end > len(seq) {
			end = len(seq)
		}
		chunkIdx := i / faiChunkSize
		compressed, err := gzipCompress(seq[i:end])
		if err != nil {
			return fmt.Errorf("compressing chunk %d of %q: %w", chunkIdx, name, err)
		}
		r.chunks[faiCacheKey{name: name, chunkIdx: chunkIdx}] = compressed
	}
	return nil
}

func (r *inMemoryFasta) Names() []string { return r.names }

func (r *inMemoryFasta) SequenceLength(name string) (int, bool) {
	l, ok := r.lengths[name]
	return l, ok
}

func (r *inMemoryFasta) GetSequenceRange(name string, start, end int) ([]byte, error) {
	r.mu.Lock()
	defer r.mu.Unlock()

	length, ok := r.lengths[name]
	if !ok {
		return nil, fmt.Errorf("reference %q not found in %s", name, r.path)
	}
	if start < 0 {
		start = 0
	}
	if end > length {
		end = length
	}
	if start >= end {
		return nil, nil
	}

	firstChunk := start / faiChunkSize
	lastChunk := (end - 1) / faiChunkSize

	result := make([]byte, 0, end-start)
	for ci := firstChunk; ci <= lastChunk; ci++ {
		chunk, err := r.decompressChunk(name, ci)
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

func (r *inMemoryFasta) GetSequence(name string) ([]byte, error) {
	length, ok := r.lengths[name]
	if !ok {
		return nil, fmt.Errorf("reference %q not found in %s", name, r.path)
	}
	return r.GetSequenceRange(name, 0, length)
}

func (r *inMemoryFasta) Close() error { return nil }

// decompressChunk returns decompressed chunk data, using the LRU cache.
// Must be called with r.mu held.
func (r *inMemoryFasta) decompressChunk(name string, chunkIdx int) ([]byte, error) {
	key := faiCacheKey{name: name, chunkIdx: chunkIdx}

	if data := r.cache.get(key); data != nil {
		return data, nil
	}

	compressed, ok := r.chunks[key]
	if !ok {
		return nil, fmt.Errorf("chunk %s:%d not found", name, chunkIdx)
	}

	gz, err := gzip.NewReader(bytes.NewReader(compressed))
	if err != nil {
		return nil, fmt.Errorf("decompressing chunk %s:%d: %w", name, chunkIdx, err)
	}
	var buf bytes.Buffer
	if _, err := buf.ReadFrom(gz); err != nil {
		gz.Close()
		return nil, fmt.Errorf("decompressing chunk %s:%d: %w", name, chunkIdx, err)
	}
	gz.Close()

	data := buf.Bytes()
	r.cache.put(key, data)
	return data, nil
}

// gzipCompress compresses data with gzip at BestSpeed.
func gzipCompress(data []byte) ([]byte, error) {
	var buf bytes.Buffer
	gz, err := gzip.NewWriterLevel(&buf, gzip.BestSpeed)
	if err != nil {
		return nil, err
	}
	if _, err := gz.Write(data); err != nil {
		return nil, err
	}
	if err := gz.Close(); err != nil {
		return nil, err
	}
	return buf.Bytes(), nil
}
