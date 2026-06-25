package vcf

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

// ChunkedVcfReader reads a numbered sequence of VCF chunk files (as produced by
// vcf-split: base.1.vcf.gz, base.2.vcf.gz, ...) as a single record stream. Only
// one underlying file is open at a time, so concatenating thousands of chunks
// stays within the file-descriptor limit. The numeric component of the supplied
// filename gives the starting index and the zero-padding width; reading walks
// upward until the next-numbered file does not exist.
type ChunkedVcfReader struct {
	dir     string
	prefix  string // filename part before the number ("" if the number is first)
	suffix  string // filename part after the number ("" if the number is last)
	padding int    // zero-pad width derived from the starting token
	index   int    // current chunk number

	cur    *VcfReader
	header *VcfHeader
}

// NewChunkedVcfReader opens the chunk sequence beginning at firstChunk. The
// filename must contain a purely-numeric, dot-delimited component (the chunk
// number), e.g. "sample.1.vcf.gz" or "sample.001.vcf.gz".
func NewChunkedVcfReader(firstChunk string) (*ChunkedVcfReader, error) {
	dir := filepath.Dir(firstChunk)
	base := filepath.Base(firstChunk)
	tokens := strings.Split(base, ".")

	numIdx := -1
	for i, tok := range tokens {
		if isAllDigits(tok) {
			numIdx = i
		}
	}
	if numIdx == -1 {
		return nil, fmt.Errorf("vcf: --chunks expects a numbered filename like base.1.vcf.gz, got %q", firstChunk)
	}
	start, _ := strconv.Atoi(tokens[numIdx])

	c := &ChunkedVcfReader{
		dir:     dir,
		prefix:  strings.Join(tokens[:numIdx], "."),
		suffix:  strings.Join(tokens[numIdx+1:], "."),
		padding: len(tokens[numIdx]),
		index:   start,
	}

	r, err := NewVcfFile(firstChunk)
	if err != nil {
		return nil, err
	}
	h, err := r.Header()
	if err != nil {
		r.Close()
		return nil, err
	}
	c.cur = r
	c.header = h
	return c, nil
}

// Header returns the header of the first chunk (chunks from one split share it).
func (c *ChunkedVcfReader) Header() (*VcfHeader, error) {
	return c.header, nil
}

// NextRecord returns the next record across the chunk sequence, transparently
// advancing to the next-numbered file when the current one is exhausted.
func (c *ChunkedVcfReader) NextRecord() (*VcfRecord, error) {
	for {
		if c.cur == nil {
			return nil, io.EOF
		}
		rec, err := c.cur.NextRecord()
		if err == nil {
			return rec, nil
		}
		if err != io.EOF {
			return nil, err
		}
		// Current chunk exhausted; try to open the next one.
		c.cur.Close()
		c.cur = nil
		c.index++
		name := c.chunkName(c.index)
		if _, statErr := os.Stat(name); statErr != nil {
			return nil, io.EOF
		}
		next, oerr := NewVcfFile(name)
		if oerr != nil {
			return nil, oerr
		}
		c.cur = next
	}
}

// Close closes the currently open chunk, if any.
func (c *ChunkedVcfReader) Close() error {
	if c.cur != nil {
		c.cur.Close()
		c.cur = nil
	}
	return nil
}

// chunkName builds the path for chunk number i, preserving the zero-padding.
func (c *ChunkedVcfReader) chunkName(i int) string {
	num := strconv.Itoa(i)
	for len(num) < c.padding {
		num = "0" + num
	}
	parts := make([]string, 0, 3)
	if c.prefix != "" {
		parts = append(parts, c.prefix)
	}
	parts = append(parts, num)
	if c.suffix != "" {
		parts = append(parts, c.suffix)
	}
	return filepath.Join(c.dir, strings.Join(parts, "."))
}

func isAllDigits(s string) bool {
	if s == "" {
		return false
	}
	for i := 0; i < len(s); i++ {
		if s[i] < '0' || s[i] > '9' {
			return false
		}
	}
	return true
}
