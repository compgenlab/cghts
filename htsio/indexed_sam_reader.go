package htsio

import (
	"encoding/binary"
	"fmt"
	"io"
	"os"

	"github.com/compgen-io/cgltk/htsio/bgzf"
)

// IndexedSamReader reads BAM files with random access via a BAI index.
// Use Query to retrieve records overlapping a genomic region.
type IndexedSamReader struct {
	ir     *bgzf.IndexedReader
	f      *os.File
	idx    *BinIndex
	refs   []bamRefInfo
	refMap map[string]int // ref name → index
	hdr    *SamHeader
}

// NewIndexedSamReader opens a BAM file and its corresponding .bai index.
// The BAI file is expected at filename.bai.
func NewIndexedSamReader(filename string) (*IndexedSamReader, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}

	// Read the header using a streaming reader first, then record where
	// the data starts for the indexed reader.
	streamReader := bgzf.NewReader(f)

	isr := &IndexedSamReader{f: f}

	// Read BAM header
	if err := isr.readHeader(streamReader); err != nil {
		f.Close()
		return nil, fmt.Errorf("bam: reading header: %w", err)
	}

	// Now create the indexed reader on the same file.
	isr.ir = bgzf.NewIndexedReader(f)

	// Load BAI index
	baiPath := filename + ".bai"
	idx, err := LoadBAI(baiPath)
	if err != nil {
		f.Close()
		return nil, fmt.Errorf("bam: loading BAI index: %w", err)
	}
	isr.idx = idx

	return isr, nil
}

// Header returns the parsed SAM header.
func (isr *IndexedSamReader) Header() (*SamHeader, error) {
	return isr.hdr, nil
}

// Close releases resources.
func (isr *IndexedSamReader) Close() error {
	if isr.f != nil {
		return isr.f.Close()
	}
	return nil
}

// Query returns a SamReader that yields records overlapping the 0-based
// half-open region [start, end) on the given reference. The returned
// reader must be consumed or discarded before calling Query again.
func (isr *IndexedSamReader) Query(ref string, start, end int) (SamReader, error) {
	refID, ok := isr.refMap[ref]
	if !ok {
		return nil, fmt.Errorf("bam: unknown reference %q", ref)
	}

	chunks := isr.idx.Query(refID, start, end)
	if len(chunks) == 0 {
		return &emptyReader{}, nil
	}

	return &regionReader{
		ir:     isr.ir,
		refs:   isr.refs,
		chunks: chunks,
		refID:  refID,
		start:  start,
		end:    end,
	}, nil
}

// readHeader reads the BAM header from a streaming bgzf reader.
func (isr *IndexedSamReader) readHeader(r *bgzf.Reader) error {
	var magic [4]byte
	if _, err := io.ReadFull(r, magic[:]); err != nil {
		return fmt.Errorf("reading magic: %w", err)
	}
	if magic != [4]byte{'B', 'A', 'M', 1} {
		return fmt.Errorf("invalid BAM magic: %x", magic)
	}

	var headerLen int32
	if err := binary.Read(r, binary.LittleEndian, &headerLen); err != nil {
		return fmt.Errorf("reading header length: %w", err)
	}

	headerText := make([]byte, headerLen)
	if _, err := io.ReadFull(r, headerText); err != nil {
		return fmt.Errorf("reading header text: %w", err)
	}

	isr.hdr = NewSamHeader()
	for _, line := range splitLines(string(headerText)) {
		if line != "" {
			isr.hdr.AddLine(line)
		}
	}

	var nRef int32
	if err := binary.Read(r, binary.LittleEndian, &nRef); err != nil {
		return fmt.Errorf("reading ref count: %w", err)
	}

	isr.refs = make([]bamRefInfo, nRef)
	isr.refMap = make(map[string]int, nRef)
	for i := int32(0); i < nRef; i++ {
		var nameLen int32
		if err := binary.Read(r, binary.LittleEndian, &nameLen); err != nil {
			return fmt.Errorf("reading ref name length: %w", err)
		}
		nameBuf := make([]byte, nameLen)
		if _, err := io.ReadFull(r, nameBuf); err != nil {
			return fmt.Errorf("reading ref name: %w", err)
		}
		name := string(nameBuf[:nameLen-1])

		var refLen int32
		if err := binary.Read(r, binary.LittleEndian, &refLen); err != nil {
			return fmt.Errorf("reading ref length: %w", err)
		}
		isr.refs[i] = bamRefInfo{name: name, length: refLen}
		isr.refMap[name] = int(i)
	}

	return nil
}

// splitLines splits a string by newlines, trimming \r.
func splitLines(s string) []string {
	var lines []string
	start := 0
	for i := 0; i < len(s); i++ {
		if s[i] == '\n' {
			line := s[start:i]
			if len(line) > 0 && line[len(line)-1] == '\r' {
				line = line[:len(line)-1]
			}
			lines = append(lines, line)
			start = i + 1
		}
	}
	if start < len(s) {
		line := s[start:]
		if len(line) > 0 && line[len(line)-1] == '\r' {
			line = line[:len(line)-1]
		}
		lines = append(lines, line)
	}
	return lines
}

// regionReader iterates BAM records within a set of chunks, filtering
// to those overlapping [start, end) on the target reference.
type regionReader struct {
	ir       *bgzf.IndexedReader
	refs     []bamRefInfo
	chunks   []Chunk
	chunkIdx int // current chunk index
	refID    int
	start    int
	end      int
	started  bool
	done     bool
}

func (rr *regionReader) Header() (*SamHeader, error) {
	return nil, nil
}

func (rr *regionReader) Close() error {
	rr.done = true
	return nil
}

func (rr *regionReader) Next() (*SamRecord, error) {
	if rr.done {
		return nil, io.EOF
	}

	for {
		if !rr.started || rr.pastChunkEnd() {
			if !rr.advanceChunk() {
				rr.done = true
				return nil, io.EOF
			}
		}

		rec, err := readBamRecord(rr.ir, rr.refs)
		if err != nil {
			if err == io.EOF {
				// Try next chunk.
				if !rr.advanceChunk() {
					rr.done = true
					return nil, io.EOF
				}
				continue
			}
			return nil, err
		}

		// Check if the record is on our target reference.
		recRefID := -1
		if rec.RefName != "*" {
			for i, r := range rr.refs {
				if r.name == rec.RefName {
					recRefID = i
					break
				}
			}
		}

		if recRefID != rr.refID {
			// Past our reference — done.
			rr.done = true
			return nil, io.EOF
		}

		// Record position (0-based).
		recStart := rec.Pos - 1
		recEnd := recStart + CigarRefLen(rec.Cigar)

		// Filter: record must overlap [start, end).
		if recEnd <= rr.start {
			continue // before query region
		}
		if recStart >= rr.end {
			// Past query region — done.
			rr.done = true
			return nil, io.EOF
		}

		return rec, nil
	}
}

func (rr *regionReader) pastChunkEnd() bool {
	if rr.chunkIdx >= len(rr.chunks) {
		return true
	}
	return rr.ir.VirtualTell() >= rr.chunks[rr.chunkIdx].End
}

func (rr *regionReader) advanceChunk() bool {
	if rr.started {
		rr.chunkIdx++
	}
	rr.started = true

	if rr.chunkIdx >= len(rr.chunks) {
		return false
	}

	if err := rr.ir.SeekToVirtualOffset(rr.chunks[rr.chunkIdx].Begin); err != nil {
		return false
	}
	return true
}

// emptyReader is a SamReader that immediately returns EOF.
type emptyReader struct{}

func (e *emptyReader) Next() (*SamRecord, error) { return nil, io.EOF }
func (e *emptyReader) Header() (*SamHeader, error) { return nil, nil }
func (e *emptyReader) Close() error                { return nil }
