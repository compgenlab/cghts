package htsio

import (
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"os"
)

// ReaderRegistration holds a format detector and constructors for use with
// the NewSamReader auto-detection registry.
type ReaderRegistration struct {
	// Detect returns true if the magic bytes match this format.
	// For raw (uncompressed) formats like CRAM, these are the first bytes
	// of the file. For compressed formats like BAM, these are the first
	// bytes after gzip decompression.
	Detect func(magic []byte) bool

	// NewFromFile creates a SamReader from a filename. The reader opens
	// its own file handle(s). This is the fast path — no nested readers.
	NewFromFile func(filename string, opts *SamReaderOpts) (SamReader, error)

	// NewFromStream creates a SamReader from a stream (e.g., stdin).
	// The io.ReadCloser has peeked bytes prepended back via MultiReader.
	// Query() is not supported on stream-based readers.
	NewFromStream func(r io.ReadCloser, opts *SamReaderOpts) (SamReader, error)
}

var readerRegistry []ReaderRegistration
var fallbackReader *ReaderRegistration

// RegisterReader adds a reader factory to the auto-detection registry.
func RegisterReader(r ReaderRegistration) {
	readerRegistry = append(readerRegistry, r)
}

// RegisterFallbackReader sets the reader used when no other reader's Detect
// matches. This is typically the SAM text reader.
func RegisterFallbackReader(r ReaderRegistration) {
	fallbackReader = &r
}

// peekSize is the number of bytes read from the file for format detection.
// Must be large enough to contain a gzip header + enough compressed data
// to decompress the first 4 bytes (for BAM magic detection).
const peekSize = 65536

// detectFormat reads the first bytes of a file/stream and returns the
// matching registry entry (or fallback). Returns the raw peek buffer
// for stream-based callers that need to prepend it back.
func detectFormat(peek []byte, n int) *ReaderRegistration {
	// Phase 1: Check raw bytes against registered readers (e.g., CRAM).
	for i := range readerRegistry {
		if readerRegistry[i].Detect(peek) {
			return &readerRegistry[i]
		}
	}

	// Phase 2: If gzip-compressed, decompress and check again (e.g., BAM).
	if n >= 2 && peek[0] == 0x1f && peek[1] == 0x8b {
		gz, gzErr := gzip.NewReader(bytes.NewReader(peek))
		if gzErr == nil {
			decompBuf := make([]byte, 16)
			dn, _ := io.ReadAtLeast(gz, decompBuf, 4)
			gz.Close()
			if dn >= 4 {
				for i := range readerRegistry {
					if readerRegistry[i].Detect(decompBuf[:dn]) {
						return &readerRegistry[i]
					}
				}
			}
		}
	}

	// Phase 3: Fallback (typically SAM text).
	return fallbackReader
}

// detectFromFile peeks at a file's magic bytes, closes it, and returns
// the matching reader registration.
func detectFromFile(filename string) (*ReaderRegistration, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	buf := make([]byte, peekSize)
	n, err := f.Read(buf)
	if err != nil && n < 2 {
		return nil, fmt.Errorf("htsio: failed to read magic bytes: %w", err)
	}

	reg := detectFormat(buf[:n], n)
	if reg == nil {
		return nil, fmt.Errorf("htsio: no registered reader for file format (peek: %x)", buf[:min(8, n)])
	}
	return reg, nil
}

// prependedReadCloser wraps an io.Reader (with peeked bytes prepended) and
// an io.Closer for the underlying stream.
type prependedReadCloser struct {
	io.Reader
	closer io.Closer
}

func (p *prependedReadCloser) Close() error {
	if p.closer != nil {
		return p.closer.Close()
	}
	return nil
}

// detectFromStream peeks at a stream's magic bytes, prepends them back,
// and returns the matching reader registration plus the reconstructed stream.
func detectFromStream(r io.ReadCloser) (*ReaderRegistration, io.ReadCloser, error) {
	buf := make([]byte, peekSize)
	n, err := io.ReadAtLeast(r, buf, 4)
	if err != nil && n < 2 {
		r.Close()
		return nil, nil, fmt.Errorf("htsio: failed to read magic bytes: %w", err)
	}
	buf = buf[:n]

	reg := detectFormat(buf, n)
	if reg == nil {
		r.Close()
		return nil, nil, fmt.Errorf("htsio: no registered reader for file format (peek: %x)", buf[:min(8, n)])
	}

	// Prepend peeked bytes back to create a full stream.
	fullReader := &prependedReadCloser{
		Reader: io.MultiReader(bytes.NewReader(buf), r),
		closer: r,
	}

	return reg, fullReader, nil
}
