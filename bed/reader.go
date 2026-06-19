package bed

import (
	"bufio"
	"compress/gzip"
	"io"
	"os"
)

// BedReader is a streaming, forward-only BED parser. Records are read one at a
// time with [BedReader.NextRecord]. A BedReader is not safe for concurrent use.
type BedReader struct {
	file         *os.File
	reader       *bufio.Reader
	closed       bool
	ignoreStrand bool
}

// NewBedReader returns a streaming BED parser over rd. The input is not
// inspected for gzip compression; wrap rd in a gzip reader yourself if needed.
// Returns [io.ErrUnexpectedEOF] if rd is nil.
func NewBedReader(rd io.Reader) (*BedReader, error) {
	if rd == nil {
		return nil, io.ErrUnexpectedEOF
	}
	return &BedReader{reader: bufio.NewReader(rd)}, nil
}

// NewBedFile opens the named BED file for streaming reads. If the file begins
// with the gzip magic bytes it is transparently decompressed. The caller
// should [BedReader.Close] the reader when done. For random access on a
// tabix-indexed file, use [NewIndexedBedReader] instead.
func NewBedFile(filename string) (*BedReader, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}

	r := bufio.NewReader(f)
	// Tolerate files shorter than two bytes (e.g. empty files): only treat the
	// input as gzip when the magic bytes are actually present.
	if magic, err := r.Peek(2); err == nil && magic[0] == 0x1f && magic[1] == 0x8b {
		gz, err := gzip.NewReader(r)
		if err != nil {
			f.Close()
			return nil, err
		}
		r = bufio.NewReader(gz)
	}

	return &BedReader{file: f, reader: r}, nil
}

// IgnoreStrand configures the reader to leave every record's strand
// unspecified ([StrandNone]) regardless of the strand column. It returns the
// reader to allow chaining at construction.
func (r *BedReader) IgnoreStrand(b bool) *BedReader {
	r.ignoreStrand = b
	return r
}

// NextRecord returns the next BED record, skipping blank lines, comment lines
// ("#"), and lines with fewer than three columns. It returns [io.EOF] when the
// input is exhausted.
func (r *BedReader) NextRecord() (*BedRecord, error) {
	if r.closed {
		return nil, io.EOF
	}
	for {
		line, err := r.reader.ReadString('\n')
		if len(line) > 0 {
			rec, ok, perr := parseBedLine(line, r.ignoreStrand)
			if perr != nil {
				return nil, perr
			}
			if ok {
				return rec, nil
			}
		}
		if err != nil {
			// On EOF (or any read error) return the error once the final line
			// has been processed above.
			return nil, err
		}
	}
}

// Close closes the underlying file (if the reader was opened from a file) and
// marks the reader closed.
func (r *BedReader) Close() {
	if r.file != nil {
		r.file.Close()
	}
	r.closed = true
}
