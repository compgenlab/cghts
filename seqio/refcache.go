package seqio

import (
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"strings"
	"sync"
)

// RefCacheReader provides random access to reference sequences using htslib's
// REF_PATH and REF_CACHE environment variable conventions.
//
// Resolution order for each sequence (by MD5):
//  1. Search REF_PATH directories for the file named by MD5
//  2. Search REF_CACHE directory (with pattern expansion)
//  3. Fetch from refget server and cache locally
//
// REF_PATH is a colon-separated list of directories and/or URLs where %s
// expands to the MD5 hash.
//
// REF_CACHE is a pattern like "/path/%2s/%2s/%s" where %2s takes the first 2
// chars of the MD5 and %s is the full MD5, creating a directory hierarchy.
//
// Cached reference files contain raw uppercase bases only (no FASTA headers,
// no newlines).
type RefCacheReader struct {
	refPath  []string          // REF_PATH entries (dirs or URL templates)
	cacheDir string            // REF_CACHE pattern (e.g., "/tmp/ref_cache/%2s/%2s/%s")
	md5s     map[string]string // sequence name → MD5
	lengths  map[string]int    // sequence name → length
	names    []string

	// In-memory chunk cache for recently accessed regions.
	cache *faiChunkCache
	mu    sync.Mutex

	// Full sequences loaded from cache files (raw bases).
	seqCache map[string][]byte
}

// NewRefCacheReader creates a ReferenceReader that resolves sequences through
// the htslib REF_PATH/REF_CACHE directory convention. The md5s map provides
// sequence name → MD5 mappings (typically from SamHeader.ReferenceMD5s()).
//
// If REF_PATH and REF_CACHE environment variables are not set, returns an error.
func NewRefCacheReader(md5s map[string]string, lengths map[string]int, names []string) (*RefCacheReader, error) {
	if len(md5s) == 0 {
		return nil, fmt.Errorf("refcache: no MD5 mappings provided (need M5 tags in @SQ header lines)")
	}

	refPath := os.Getenv("REF_PATH")
	refCache := os.Getenv("REF_CACHE")

	if refPath == "" && refCache == "" {
		return nil, fmt.Errorf("refcache: neither REF_PATH nor REF_CACHE environment variables are set")
	}

	var paths []string
	if refPath != "" {
		paths = strings.Split(refPath, ":")
	}

	r := &RefCacheReader{
		refPath:  paths,
		cacheDir: refCache,
		md5s:     md5s,
		lengths:  lengths,
		names:    names,
		cache:    newFaiChunkCache(faiCacheMaxSize),
		seqCache: make(map[string][]byte),
	}

	return r, nil
}

func (r *RefCacheReader) Names() []string { return r.names }

func (r *RefCacheReader) SequenceLength(name string) (int, bool) {
	l, ok := r.lengths[name]
	return l, ok
}

func (r *RefCacheReader) GetSequenceRange(name string, start, end int) ([]byte, error) {
	r.mu.Lock()
	defer r.mu.Unlock()

	length, ok := r.lengths[name]
	if !ok {
		return nil, fmt.Errorf("refcache: sequence %q not found", name)
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

	// Load full sequence if not cached.
	seq, err := r.loadSequence(name)
	if err != nil {
		return nil, err
	}

	if end > len(seq) {
		end = len(seq)
	}
	if start >= end {
		return nil, nil
	}

	return seq[start:end], nil
}

func (r *RefCacheReader) GetSequence(name string) ([]byte, error) {
	r.mu.Lock()
	defer r.mu.Unlock()

	seq, err := r.loadSequence(name)
	if err != nil {
		return nil, err
	}
	result := make([]byte, len(seq))
	copy(result, seq)
	return result, nil
}

func (r *RefCacheReader) Close() error { return nil }

// loadSequence resolves and loads a sequence by name. Must be called with mu held.
func (r *RefCacheReader) loadSequence(name string) ([]byte, error) {
	if seq, ok := r.seqCache[name]; ok {
		return seq, nil
	}

	md5, ok := r.md5s[name]
	if !ok {
		return nil, fmt.Errorf("refcache: no MD5 for sequence %q", name)
	}

	// 1. Search REF_PATH entries.
	for _, entry := range r.refPath {
		data, err := r.resolveRefPathEntry(entry, md5)
		if err == nil && len(data) > 0 {
			r.seqCache[name] = data
			return data, nil
		}
	}

	// 2. Search REF_CACHE directory.
	if r.cacheDir != "" {
		cachePath := expandCachePattern(r.cacheDir, md5)
		data, err := os.ReadFile(cachePath)
		if err == nil && len(data) > 0 {
			data = uppercaseBases(data)
			r.seqCache[name] = data
			return data, nil
		}
	}

	// 3. Try fetching from refget and caching.
	data, err := r.fetchAndCache(name, md5)
	if err != nil {
		return nil, fmt.Errorf("refcache: could not resolve %q (MD5 %s): %w", name, md5, err)
	}

	r.seqCache[name] = data
	return data, nil
}

// resolveRefPathEntry resolves a single REF_PATH entry. The entry may be a
// directory, a URL template (with %s for MD5), or a plain directory path.
func (r *RefCacheReader) resolveRefPathEntry(entry, md5 string) ([]byte, error) {
	if strings.HasPrefix(entry, "http://") || strings.HasPrefix(entry, "https://") {
		// URL template: expand %s to MD5 and fetch.
		url := expandPathTemplate(entry, md5)
		return fetchURL(url)
	}

	// Local directory: look for file named by MD5.
	path := expandPathTemplate(entry, md5)
	if !strings.Contains(entry, "%s") {
		// Plain directory — append MD5 as filename.
		path = filepath.Join(entry, md5)
	}

	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	return uppercaseBases(data), nil
}

// fetchAndCache fetches a sequence from the default refget server and writes
// it to the REF_CACHE directory.
func (r *RefCacheReader) fetchAndCache(name, md5 string) ([]byte, error) {
	// Fetch from EBI refget endpoint.
	url := fmt.Sprintf("%s/%s", DefaultRefgetServer, md5)
	req, err := http.NewRequest("GET", url, nil)
	if err != nil {
		return nil, err
	}
	req.Header.Set("Accept", "text/plain")

	resp, err := httpClient.Do(req)
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return nil, fmt.Errorf("HTTP %d from refget for %s", resp.StatusCode, name)
	}

	data, err := readCappedBody(resp.Body)
	if err != nil {
		return nil, fmt.Errorf("refget fetch %s: %w", name, err)
	}

	data = uppercaseBases(data)

	// Cache to REF_CACHE if configured.
	if r.cacheDir != "" {
		cachePath := expandCachePattern(r.cacheDir, md5)
		dir := filepath.Dir(cachePath)
		// Best-effort cache write — a failure here must not fail the fetch, but
		// we surface it on stderr rather than swallowing it silently so a
		// misconfigured or unwritable REF_CACHE is diagnosable.
		if err := os.MkdirAll(dir, 0755); err != nil {
			fmt.Fprintf(os.Stderr, "warning: could not create REF_CACHE dir %q: %v\n", dir, err)
		} else if err := os.WriteFile(cachePath, data, 0644); err != nil {
			fmt.Fprintf(os.Stderr, "warning: could not write REF_CACHE file %q: %v\n", cachePath, err)
		}
	}

	return data, nil
}

// expandCachePattern expands a REF_CACHE pattern like "/path/%2s/%2s/%s"
// using the MD5 hash. Each %Ns takes the next N chars of the MD5 (advancing
// a cursor), and %s takes the full MD5. This matches htslib's behavior where
// %2s/%2s/%s with MD5 "abcdef..." produces "ab/cd/abcdef...".
func expandCachePattern(pattern, md5 string) string {
	var result strings.Builder
	cursor := 0 // position in md5 for %Ns expansions
	i := 0
	for i < len(pattern) {
		if pattern[i] == '%' && i+1 < len(pattern) {
			if pattern[i+1] == 's' {
				// %s → full MD5
				result.WriteString(md5)
				i += 2
			} else if i+2 < len(pattern) && pattern[i+2] == 's' {
				// %Ns → next N chars of MD5
				n := int(pattern[i+1] - '0')
				if n > 0 && cursor+n <= len(md5) {
					result.WriteString(md5[cursor : cursor+n])
					cursor += n
				}
				i += 3
			} else {
				result.WriteByte(pattern[i])
				i++
			}
		} else {
			result.WriteByte(pattern[i])
			i++
		}
	}
	return result.String()
}

// expandPathTemplate expands %s in a path template to the MD5 hash.
func expandPathTemplate(tmpl, md5 string) string {
	return strings.ReplaceAll(tmpl, "%s", md5)
}

// fetchURL fetches a URL and returns the body as bytes.
func fetchURL(url string) ([]byte, error) {
	resp, err := httpClient.Get(url)
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return nil, fmt.Errorf("HTTP %d from %s", resp.StatusCode, url)
	}

	data, err := readCappedBody(resp.Body)
	if err != nil {
		return nil, fmt.Errorf("fetch %s: %w", url, err)
	}
	return uppercaseBases(data), nil
}

// maxRefSequenceBytes caps a single reference-sequence download. A refget
// server (or REF_PATH source) returns one sequence — a contig/chromosome — per
// request, so this is a generous ceiling above any realistic single sequence;
// its purpose is to stop a malicious or buggy server from streaming an
// unbounded body into memory, not to enforce a tight size limit.
const maxRefSequenceBytes = 8 << 30 // 8 GiB

// readCappedBody reads an HTTP body for a single reference sequence, capped at
// maxRefSequenceBytes. It reads one byte past the cap so an over-large body is
// reported as an error rather than silently truncated.
func readCappedBody(r io.Reader) ([]byte, error) {
	data, err := io.ReadAll(io.LimitReader(r, maxRefSequenceBytes+1))
	if err != nil {
		return nil, err
	}
	if int64(len(data)) > maxRefSequenceBytes {
		return nil, fmt.Errorf("reference sequence exceeds %d bytes", maxRefSequenceBytes)
	}
	return data, nil
}

// uppercaseBases strips whitespace and uppercases DNA bases.
func uppercaseBases(data []byte) []byte {
	result := make([]byte, 0, len(data))
	for _, b := range data {
		if b == '\n' || b == '\r' || b == ' ' {
			continue
		}
		if b >= 'a' && b <= 'z' {
			b -= 32
		}
		result = append(result, b)
	}
	return result
}
