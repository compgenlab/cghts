package htsio

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"iter"
	"os"
	"strings"
)

// SamTextReader reads plain SAM text files (.sam or .sam.gz).
// It implements the SamReader interface.
type SamTextReader struct {
	src       io.ReadCloser
	gz        *gzip.Reader
	scanner   *bufio.Scanner
	header    *SamHeader
	opts      *SamReaderOpts
	nextLine  string // buffered first data line from header parsing
	err       error
}

// NewSamTextReader creates a SAM text reader for the given file.
func NewSamTextReader(filename string, opts ...*SamReaderOpts) (*SamTextReader, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}

	var o *SamReaderOpts
	if len(opts) > 0 {
		o = opts[0]
	} else {
		o = NewSamReaderOpts()
	}

	r := &SamTextReader{
		src:  f,
		opts: o,
	}

	var reader io.Reader = f
	if strings.HasSuffix(filename, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			f.Close()
			return nil, fmt.Errorf("sam: gzip reader: %w", err)
		}
		r.gz = gz
		reader = gz
	}

	r.scanner = bufio.NewScanner(reader)
	r.scanner.Buffer(make([]byte, 0, 1024*1024), 10*1024*1024)
	r.populateHeader()

	return r, nil
}

func (r *SamTextReader) populateHeader() {
	for r.scanner.Scan() {
		line := r.scanner.Text()
		if strings.HasPrefix(line, "@") {
			if r.header == nil {
				r.header = NewSamHeader()
			}
			r.header.AddLine(line)
			continue
		}
		r.nextLine = line
		return
	}
}

// Header returns the parsed SAM header.
func (r *SamTextReader) Header() (*SamHeader, error) {
	return r.header, nil
}

// Next returns the next SamRecord. Returns nil, io.EOF when done.
func (r *SamTextReader) Next() (*SamRecord, error) {
	if r.err != nil {
		return nil, r.err
	}

	// Return the buffered first data line.
	if r.nextLine != "" {
		rec, err := parseSamLine(r.nextLine)
		r.nextLine = ""
		if err != nil {
			return nil, fmt.Errorf("parse SAM: %w", err)
		}
		if r.passesFilters(rec) {
			return rec, nil
		}
	}

	for r.scanner.Scan() {
		line := r.scanner.Text()
		if strings.HasPrefix(line, "@") {
			if r.header == nil {
				r.header = NewSamHeader()
			}
			r.header.AddLine(line)
			continue
		}
		rec, err := parseSamLine(line)
		if err != nil {
			return nil, fmt.Errorf("parse SAM: %w", err)
		}
		if r.passesFilters(rec) {
			return rec, nil
		}
	}

	if err := r.scanner.Err(); err != nil {
		r.err = err
		return nil, err
	}

	r.err = io.EOF
	return nil, io.EOF
}

func (r *SamTextReader) passesFilters(rec *SamRecord) bool {
	if r.opts.flagReq != 0 && rec.Flag&r.opts.flagReq != r.opts.flagReq {
		return false
	}
	if r.opts.flagFilter != 0 && rec.Flag&r.opts.flagFilter != 0 {
		return false
	}
	if r.opts.minMapQ != 0 && rec.MapQ < r.opts.minMapQ {
		return false
	}
	for _, f := range r.opts.tagFilters {
		if !f.matchesRecord(rec) {
			return false
		}
	}
	return true
}

// Query is not supported for plain SAM files (no index).
func (r *SamTextReader) Query(ref string, start, end int) (iter.Seq2[*SamRecord, error], error) {
	return nil, fmt.Errorf("sam: Query not supported on plain SAM files (no index)")
}

// Close releases resources.
func (r *SamTextReader) Close() error {
	if r.gz != nil {
		r.gz.Close()
	}
	if r.src != nil {
		return r.src.Close()
	}
	return nil
}
