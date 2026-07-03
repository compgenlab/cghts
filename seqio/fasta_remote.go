package seqio

import (
	"bufio"
	"fmt"
	"io"
	"net/http"
	"strconv"
	"strings"
	"sync"

	"github.com/compgenlab/hts/iosource"
)

// RemoteFastaReader provides random access to a remote indexed FASTA file
// via HTTP/HTTPS. The .fai index is downloaded once on open. Sequence data
// is fetched in 10MB chunks using HTTP Range requests and cached with the
// same LRU strategy as IndexedFastaReader.
//
// Requires the server to support Range requests (most HTTP servers do).
// Supports both plain FASTA and bgzip-compressed FASTA (bgzip uses the
// same .fai byte offsets as uncompressed for the Range calculation, but
// the fetched bytes need decompression — not yet supported; plain FASTA only).
type RemoteFastaReader struct {
	url   string
	src   *iosource.HTTPRange
	fai   map[string]*FaiEntry
	names []string

	cache *faiChunkCache
	mu    sync.Mutex
}

// NewRemoteFastaReader opens a remote FASTA file for random access.
// The .fai index is fetched from url+".fai" and downloaded fully.
// Sequence chunks are fetched on demand via HTTP Range requests.
func NewRemoteFastaReader(url string) (*RemoteFastaReader, error) {
	// Fetch .fai index.
	faiURL := url + ".fai"
	fai, names, err := fetchFaiIndex(faiURL)
	if err != nil {
		return nil, fmt.Errorf("fetching remote .fai index: %w", err)
	}

	src, err := iosource.NewHTTPRange(url)
	if err != nil {
		return nil, fmt.Errorf("opening remote FASTA: %w", err)
	}

	return &RemoteFastaReader{
		url:   url,
		src:   src,
		fai:   fai,
		names: names,
		cache: newFaiChunkCache(faiCacheMaxSize),
	}, nil
}

// Names returns the ordered sequence names from the remote .fai index.
func (r *RemoteFastaReader) Names() []string { return r.names }

// SequenceLength returns the length of the named sequence, or false if it is
// not present in the .fai index.
func (r *RemoteFastaReader) SequenceLength(name string) (int, bool) {
	entry, ok := r.fai[name]
	if !ok {
		return 0, false
	}
	return entry.Length, true
}

// GetSequenceRange returns the bases for [start, end) (0-based, half-open) of
// the named sequence, fetching any needed 10MB chunks via HTTP Range requests
// and caching them. The result is uppercased; coordinates are clamped to the
// sequence bounds. It is safe for concurrent use.
func (r *RemoteFastaReader) GetSequenceRange(name string, start, end int) ([]byte, error) {
	r.mu.Lock()
	defer r.mu.Unlock()

	entry, ok := r.fai[name]
	if !ok {
		return nil, fmt.Errorf("sequence %q not found in remote .fai", name)
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

// GetSequence returns the full named sequence, uppercased. Prefer
// [RemoteFastaReader.GetSequenceRange] when only a sub-range is needed.
func (r *RemoteFastaReader) GetSequence(name string) ([]byte, error) {
	entry, ok := r.fai[name]
	if !ok {
		return nil, fmt.Errorf("sequence %q not found in remote .fai", name)
	}
	return r.GetSequenceRange(name, 0, entry.Length)
}

// Close releases resources held by the reader.
func (r *RemoteFastaReader) Close() error { return r.src.Close() }

// loadChunk fetches a single chunk via HTTP Range request, using the cache.
// Must be called with r.mu held.
func (r *RemoteFastaReader) loadChunk(name string, chunkIdx int, entry *FaiEntry) ([]byte, error) {
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

	// Compute byte range in the FASTA file.
	startByte := entry.Offset + int64(baseStart/entry.LineBases)*int64(entry.LineWidth) + int64(baseStart%entry.LineBases)
	lastBase := baseEnd - 1
	endByte := entry.Offset + int64(lastBase/entry.LineBases)*int64(entry.LineWidth) + int64(lastBase%entry.LineBases) + 1

	// Fetch via HTTP Range request. A read that reaches EOF returns the bytes
	// available (a trailing partial line) rather than an error.
	buf := make([]byte, endByte-startByte)
	n, err := r.src.ReadAt(buf, startByte)
	if err != nil && err != io.EOF {
		return nil, fmt.Errorf("fetching %s chunk %d: %w", name, chunkIdx, err)
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

// fetchFaiIndex downloads and parses a .fai index from a URL.
func fetchFaiIndex(url string) (map[string]*FaiEntry, []string, error) {
	resp, err := httpClient.Get(url)
	if err != nil {
		return nil, nil, err
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return nil, nil, fmt.Errorf("HTTP %d fetching %s", resp.StatusCode, url)
	}

	fai := make(map[string]*FaiEntry)
	var names []string
	scanner := bufio.NewScanner(resp.Body)
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
		return nil, nil, fmt.Errorf("empty .fai from %s", url)
	}
	return fai, names, nil
}
