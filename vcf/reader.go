package vcf

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strings"
)

// VcfReader is a streaming, forward-only VCF parser. Records are read one at a
// time with [VcfReader.NextRecord]. A VcfReader is not safe for concurrent use.
type VcfReader struct {
	file       *os.File
	reader     *bufio.Reader
	header     *VcfHeader
	headerRead bool
	closed     bool
}

// NewVcfReader returns a streaming VCF parser over rd. The input is not
// inspected for gzip compression; wrap rd in a gzip reader yourself if needed.
func NewVcfReader(rd io.Reader) (*VcfReader, error) {
	if rd == nil {
		return nil, io.ErrUnexpectedEOF
	}
	return &VcfReader{reader: bufio.NewReaderSize(rd, 1024*1024)}, nil
}

// NewVcfFile opens the named VCF file for streaming reads. If the file begins
// with the gzip magic bytes it is transparently decompressed (this also handles
// BGZF, which is gzip-compatible). The caller should [VcfReader.Close] the
// reader when done. For random access on a tabix-indexed file, use
// [NewIndexedVcfReader] instead.
func NewVcfFile(filename string) (*VcfReader, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}

	r := bufio.NewReaderSize(f, 1024*1024)
	if magic, err := r.Peek(2); err == nil && magic[0] == 0x1f && magic[1] == 0x8b {
		gz, err := gzip.NewReader(r)
		if err != nil {
			f.Close()
			return nil, err
		}
		r = bufio.NewReaderSize(gz, 1024*1024)
	}

	return &VcfReader{file: f, reader: r}, nil
}

// Header returns the parsed VCF header, reading it from the input on first call.
func (r *VcfReader) Header() (*VcfHeader, error) {
	if !r.headerRead {
		if err := r.readHeader(); err != nil {
			return nil, err
		}
	}
	return r.header, nil
}

func (r *VcfReader) readHeader() error {
	var meta []string
	var chromLine string
	for {
		line, err := r.reader.ReadString('\n')
		trimmed := strings.TrimRight(line, "\r\n")
		if strings.HasPrefix(trimmed, "#CHROM") {
			chromLine = trimmed
			break
		}
		if strings.HasPrefix(trimmed, "##") {
			meta = append(meta, trimmed)
		}
		if err != nil {
			if err == io.EOF {
				if chromLine == "" && len(meta) == 0 {
					return fmt.Errorf("vcf: missing header")
				}
				break
			}
			return err
		}
	}

	h, err := parseHeaderLines(meta, chromLine)
	if err != nil {
		return err
	}
	r.header = h
	r.headerRead = true
	return nil
}

// NextRecord returns the next VCF record, skipping blank and comment lines. It
// reads the header first if it has not already been read. It returns [io.EOF]
// when the input is exhausted.
func (r *VcfReader) NextRecord() (*VcfRecord, error) {
	if r.closed {
		return nil, io.EOF
	}
	if !r.headerRead {
		if err := r.readHeader(); err != nil {
			return nil, err
		}
	}
	for {
		line, err := r.reader.ReadString('\n')
		if len(line) > 0 {
			trimmed := strings.TrimRight(line, "\r\n")
			if trimmed != "" && trimmed[0] != '#' {
				rec, perr := newRecord(trimmed, r.header)
				if perr != nil {
					return nil, perr
				}
				return rec, nil
			}
		}
		if err != nil {
			return nil, err
		}
	}
}

// Close closes the underlying file (if the reader was opened from a file) and
// marks the reader closed.
func (r *VcfReader) Close() {
	if r.file != nil {
		r.file.Close()
	}
	r.closed = true
}
