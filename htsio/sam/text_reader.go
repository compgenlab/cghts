package sam

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"iter"
	"os"
	"strconv"
	"strings"

	"github.com/compgen-io/cgkit/htsio"
)

func init() {
	// SAM text reader is the fallback for any format not matched by other
	// registered readers (BAM, CRAM, etc.).
	htsio.RegisterFallbackReader(htsio.ReaderRegistration{
		NewFromFile: func(filename string, opts *htsio.SamReaderOpts) (htsio.SamReader, error) {
			return NewTextReader(filename, opts)
		},
		NewFromStream: func(r io.ReadCloser, opts *htsio.SamReaderOpts) (htsio.SamReader, error) {
			return NewTextReaderFromReader(r, opts)
		},
	})
}

// TextReader reads plain SAM text files (.sam or .sam.gz).
// It implements the htsio.SamReader interface.
type TextReader struct {
	src      io.ReadCloser
	gz       *gzip.Reader
	scanner  *bufio.Scanner
	header   *htsio.SamHeader
	opts     *htsio.SamReaderOpts
	nextLine string
	err      error
	isGzip   bool
}

// NewTextReader creates a SAM text reader for the given file.
func NewTextReader(filename string, opts *htsio.SamReaderOpts) (*TextReader, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	return NewTextReaderFromReader(f, opts)
}

// NewTextReaderFromReader creates a SAM text reader from an io.ReadCloser.
// The reader may be plain text or gzip-compressed; gzip is auto-detected
// by peeking at the first two bytes.
func NewTextReaderFromReader(rc io.ReadCloser, opts *htsio.SamReaderOpts) (*TextReader, error) {
	if opts == nil {
		opts = htsio.NewSamReaderOpts()
	}

	r := &TextReader{
		src:  rc,
		opts: opts,
	}

	// Peek to detect gzip.
	br := bufio.NewReaderSize(rc, 64*1024)
	peek, _ := br.Peek(2)

	var reader io.Reader = br
	if len(peek) >= 2 && peek[0] == 0x1f && peek[1] == 0x8b {
		gz, err := gzip.NewReader(br)
		if err != nil {
			rc.Close()
			return nil, fmt.Errorf("sam: gzip reader: %w", err)
		}
		r.gz = gz
		r.isGzip = true
		reader = gz
	}

	r.scanner = bufio.NewScanner(reader)
	r.scanner.Buffer(make([]byte, 0, 1024*1024), 10*1024*1024)
	r.populateHeader()

	return r, nil
}

func (r *TextReader) populateHeader() {
	for r.scanner.Scan() {
		line := r.scanner.Text()
		if strings.HasPrefix(line, "@") {
			if r.header == nil {
				r.header = htsio.NewSamHeader()
			}
			r.header.AddLine(line)
			continue
		}
		r.nextLine = line
		return
	}
}

func (r *TextReader) Header() (*htsio.SamHeader, error) {
	return r.header, nil
}

func (r *TextReader) Records() iter.Seq2[*htsio.SamRecord, error] {
	return func(yield func(*htsio.SamRecord, error) bool) {
		if r.err != nil {
			yield(nil, r.err)
			return
		}

		if r.nextLine != "" {
			rec, err := parseSamLine(r.nextLine)
			r.nextLine = ""
			if err != nil {
				yield(nil, fmt.Errorf("parse SAM: %w", err))
				return
			}
			if r.opts.PassesFilters(rec) {
				if !yield(rec, nil) {
					return
				}
			}
		}

		for r.scanner.Scan() {
			line := r.scanner.Text()
			if strings.HasPrefix(line, "@") {
				if r.header == nil {
					r.header = htsio.NewSamHeader()
				}
				r.header.AddLine(line)
				continue
			}
			rec, err := parseSamLine(line)
			if err != nil {
				yield(nil, fmt.Errorf("parse SAM: %w", err))
				return
			}
			if r.opts.PassesFilters(rec) {
				if !yield(rec, nil) {
					return
				}
			}
		}

		if err := r.scanner.Err(); err != nil {
			r.err = err
			yield(nil, err)
			return
		}

		r.err = io.EOF
	}
}

func (r *TextReader) Query(ref string, start, end int) (iter.Seq2[*htsio.SamRecord, error], error) {
	return nil, fmt.Errorf("sam: Query not supported on plain SAM files (no index)")
}

func (r *TextReader) Close() error {
	if r.gz != nil {
		r.gz.Close()
	}
	if r.src != nil {
		return r.src.Close()
	}
	return nil
}

func parseSamLine(line string) (*htsio.SamRecord, error) {
	fields := strings.SplitN(line, "\t", 12)
	if len(fields) < 11 {
		return nil, fmt.Errorf("expected at least 11 fields, got %d", len(fields))
	}

	flag, err := strconv.Atoi(fields[1])
	if err != nil {
		return nil, fmt.Errorf("invalid flag %q: %w", fields[1], err)
	}

	pos, err := strconv.Atoi(fields[3])
	if err != nil {
		return nil, fmt.Errorf("invalid pos %q: %w", fields[3], err)
	}

	mapq, err := strconv.Atoi(fields[4])
	if err != nil {
		return nil, fmt.Errorf("invalid mapq %q: %w", fields[4], err)
	}

	pnext, err := strconv.Atoi(fields[7])
	if err != nil {
		return nil, fmt.Errorf("invalid pnext %q: %w", fields[7], err)
	}

	tlen, err := strconv.Atoi(fields[8])
	if err != nil {
		return nil, fmt.Errorf("invalid tlen %q: %w", fields[8], err)
	}

	rec := &htsio.SamRecord{
		ReadName:  fields[0],
		Flag:      flag,
		RefName:   fields[2],
		Pos:       pos,
		MapQ:      mapq,
		Cigar:     fields[5],
		RefNext:   fields[6],
		PosNext:   pnext,
		InsertLen: tlen,
		Seq:       fields[9],
		Qual:      fields[10],
		Tags:      make(map[string]htsio.SamTag),
	}

	if len(fields) > 11 {
		for _, raw := range strings.Split(fields[11], "\t") {
			parts := strings.SplitN(raw, ":", 3)
			// Require TAG:TYPE:VALUE with a non-empty single-character type.
			// Malformed fields are skipped (samtools is similarly tolerant);
			// the type-length check also guards the parts[1][0] access below.
			if len(parts) != 3 || len(parts[1]) == 0 {
				continue
			}
			tag := parts[0]
			rec.Tags[tag] = htsio.SamTag{
				Type:  parts[1][0],
				Value: parts[2],
			}
			rec.TagOrder = append(rec.TagOrder, tag)
		}
	}

	return rec, nil
}
