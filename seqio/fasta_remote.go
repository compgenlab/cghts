package seqio

import (
	"bufio"
	"fmt"
	"io"
	"net/http"
	"strconv"
	"strings"
	"sync"
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

	return &RemoteFastaReader{
		url:   url,
		fai:   fai,
		names: names,
		cache: newFaiChunkCache(faiCacheMaxSize),
	}, nil
}

func (r *RemoteFastaReader) Names() []string { return r.names }

func (r *RemoteFastaReader) SequenceLength(name string) (int, bool) {
	entry, ok := r.fai[name]
	if !ok {
		return 0, false
	}
	return entry.Length, true
}

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

func (r *RemoteFastaReader) GetSequence(name string) ([]byte, error) {
	entry, ok := r.fai[name]
	if !ok {
		return nil, fmt.Errorf("sequence %q not found in remote .fai", name)
	}
	return r.GetSequenceRange(name, 0, entry.Length)
}

func (r *RemoteFastaReader) Close() error { return nil }

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

	// Fetch via HTTP Range request.
	buf, err := httpRangeGet(r.url, startByte, endByte-1) // Range is inclusive
	if err != nil {
		return nil, fmt.Errorf("fetching %s chunk %d: %w", name, chunkIdx, err)
	}

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
	resp, err := http.Get(url)
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

// httpRangeGet fetches bytes [start, end] (inclusive) from a URL using an HTTP Range request.
func httpRangeGet(url string, start, end int64) ([]byte, error) {
	req, err := http.NewRequest("GET", url, nil)
	if err != nil {
		return nil, err
	}
	req.Header.Set("Range", fmt.Sprintf("bytes=%d-%d", start, end))

	resp, err := http.DefaultClient.Do(req)
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusPartialContent && resp.StatusCode != http.StatusOK {
		return nil, fmt.Errorf("HTTP %d (expected 206 Partial Content)", resp.StatusCode)
	}

	return io.ReadAll(resp.Body)
}
