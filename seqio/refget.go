package seqio

import (
	"fmt"
	"io"
	"net/http"
	"strconv"
	"strings"
	"sync"
)

const (
	// Default EBI refget endpoint for CRAM reference sequences.
	DefaultRefgetServer = "https://www.ebi.ac.uk/ena/cram/md5"
)

// RefgetReader provides random access to reference sequences via the GA4GH
// refget API. Sequences are identified by their MD5 checksum, which comes
// from the M5 tag in SAM @SQ header lines.
//
// Sequences are fetched in chunks and cached with the same LRU strategy as
// other ReferenceReader implementations.
type RefgetReader struct {
	server  string            // base URL, e.g. "https://www.ebi.ac.uk/ena/cram/md5"
	md5s    map[string]string // sequence name → MD5
	lengths map[string]int    // sequence name → length (from header or metadata)
	names   []string          // ordered sequence names

	cache *faiChunkCache
	mu    sync.Mutex
}

// RefgetOption configures a RefgetReader.
type RefgetOption func(*RefgetReader)

// RefgetServer sets a custom refget server URL.
func RefgetServer(url string) RefgetOption {
	return func(r *RefgetReader) { r.server = strings.TrimRight(url, "/") }
}

// NewRefgetReader creates a ReferenceReader that fetches sequences from a
// GA4GH refget server. The md5s map provides sequence name → MD5 mappings
// (typically from SamHeader.ReferenceMD5s()). The lengths map provides
// sequence name → length (typically from @SQ LN tags).
//
// If no server is specified via RefgetServer option, the default EBI
// endpoint is used.
func NewRefgetReader(md5s map[string]string, lengths map[string]int, names []string, opts ...RefgetOption) (*RefgetReader, error) {
	if len(md5s) == 0 {
		return nil, fmt.Errorf("refget: no MD5 mappings provided (need M5 tags in @SQ header lines)")
	}

	r := &RefgetReader{
		server:  DefaultRefgetServer,
		md5s:    md5s,
		lengths: lengths,
		names:   names,
		cache:   newFaiChunkCache(faiCacheMaxSize),
	}
	for _, opt := range opts {
		opt(r)
	}

	return r, nil
}

func (r *RefgetReader) Names() []string { return r.names }

func (r *RefgetReader) SequenceLength(name string) (int, bool) {
	l, ok := r.lengths[name]
	return l, ok
}

func (r *RefgetReader) GetSequenceRange(name string, start, end int) ([]byte, error) {
	r.mu.Lock()
	defer r.mu.Unlock()

	md5, ok := r.md5s[name]
	if !ok {
		return nil, fmt.Errorf("refget: no MD5 for sequence %q", name)
	}

	length, ok := r.lengths[name]
	if !ok {
		// Try to get length from metadata endpoint.
		l, err := r.fetchMetadataLength(md5)
		if err != nil {
			return nil, fmt.Errorf("refget: cannot determine length of %q: %w", name, err)
		}
		r.lengths[name] = l
		length = l
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
		chunk, err := r.loadChunk(name, md5, ci, length)
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

func (r *RefgetReader) GetSequence(name string) ([]byte, error) {
	length, ok := r.lengths[name]
	if !ok {
		return nil, fmt.Errorf("refget: sequence %q not found", name)
	}
	return r.GetSequenceRange(name, 0, length)
}

func (r *RefgetReader) Close() error { return nil }

// loadChunk fetches a chunk from the refget server, using the cache.
// Must be called with r.mu held.
func (r *RefgetReader) loadChunk(name, md5 string, chunkIdx, seqLen int) ([]byte, error) {
	key := faiCacheKey{name: name, chunkIdx: chunkIdx}

	if data := r.cache.get(key); data != nil {
		return data, nil
	}

	baseStart := chunkIdx * faiChunkSize
	baseEnd := baseStart + faiChunkSize
	if baseEnd > seqLen {
		baseEnd = seqLen
	}

	// Refget API: GET /sequence/{md5}?start=X&end=Y (0-based, end exclusive)
	url := fmt.Sprintf("%s/%s?start=%d&end=%d", r.server, md5, baseStart, baseEnd)

	req, err := http.NewRequest("GET", url, nil)
	if err != nil {
		return nil, err
	}
	req.Header.Set("Accept", "text/plain")

	resp, err := httpClient.Do(req)
	if err != nil {
		return nil, fmt.Errorf("refget fetch %s [%d,%d): %w", name, baseStart, baseEnd, err)
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK && resp.StatusCode != http.StatusPartialContent {
		return nil, fmt.Errorf("refget fetch %s [%d,%d): HTTP %d", name, baseStart, baseEnd, resp.StatusCode)
	}

	// Cap the read to the number of bases requested; the server should return
	// exactly this subsequence, so this bounds memory if it returns more.
	data, err := io.ReadAll(io.LimitReader(resp.Body, int64(baseEnd-baseStart)))
	if err != nil {
		return nil, fmt.Errorf("refget fetch %s [%d,%d): %w", name, baseStart, baseEnd, err)
	}

	// Refget returns raw uppercase bases (no newlines), but uppercase just in case.
	for i, b := range data {
		if b >= 'a' && b <= 'z' {
			data[i] = b - 32
		}
	}

	r.cache.put(key, data)
	return data, nil
}

// fetchMetadataLength retrieves sequence length from the refget metadata endpoint.
func (r *RefgetReader) fetchMetadataLength(md5 string) (int, error) {
	url := fmt.Sprintf("%s/%s/metadata", r.server, md5)

	resp, err := httpClient.Get(url)
	if err != nil {
		return 0, err
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return 0, fmt.Errorf("HTTP %d from %s", resp.StatusCode, url)
	}

	// The metadata response is JSON with a "metadata.length" field.
	// We do a simple parse to avoid importing encoding/json for one field.
	// Cap the read — this is a small JSON document, so a multi-megabyte body
	// indicates a misbehaving server, not legitimate metadata.
	body, err := io.ReadAll(io.LimitReader(resp.Body, 1<<20))
	if err != nil {
		return 0, err
	}
	return parseRefgetLength(string(body))
}

// parseRefgetLength extracts the length from a refget metadata JSON response.
// Looks for "length": NNN pattern.
func parseRefgetLength(body string) (int, error) {
	idx := strings.Index(body, `"length"`)
	if idx < 0 {
		return 0, fmt.Errorf("no length field in refget metadata")
	}
	rest := body[idx+len(`"length"`):]
	// Skip whitespace and colon.
	rest = strings.TrimLeft(rest, " \t\n\r:")
	// Read digits.
	end := 0
	for end < len(rest) && rest[end] >= '0' && rest[end] <= '9' {
		end++
	}
	if end == 0 {
		return 0, fmt.Errorf("could not parse length from refget metadata")
	}
	return strconv.Atoi(rest[:end])
}
