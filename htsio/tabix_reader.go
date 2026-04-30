package htsio

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/compgen-io/cgltk/htsio/bgzf"
)

// TabixReader reads BGZF-compressed, tabix-indexed text files (BED, VCF, GFF,
// etc.) with random access by genomic region. The .tbi index is required — it
// provides column definitions, the coordinate system (0-based vs 1-based),
// comment character, and header skip count.
type TabixReader struct {
	ir  *bgzf.IndexedReader
	f   *os.File
	idx *BinIndex
}

// TabixRecord holds a single parsed line from a tabix query along with the
// extracted genomic coordinates.
type TabixRecord struct {
	Line  string
	Ref   string
	Start int // 0-based
	End   int // 0-based, exclusive
}

// TabixIterator yields lines from a tabix query region.
type TabixIterator struct {
	ir      *bgzf.IndexedReader
	idx     *BinIndex
	scanner *bufio.Scanner
	chunks  []Chunk
	ref     string
	start   int // query start, 0-based
	end     int // query end, 0-based exclusive
	started bool
	done    bool
}

// NewTabixReader opens a BGZF-compressed file and its .tbi index.
// The TBI index file must exist at filename.tbi.
func NewTabixReader(filename string) (*TabixReader, error) {
	tbiPath := filename + ".tbi"
	idx, err := LoadTBI(tbiPath)
	if err != nil {
		return nil, fmt.Errorf("tabix: loading index: %w", err)
	}

	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}

	ir := bgzf.NewIndexedReader(f)

	return &TabixReader{
		ir:  ir,
		f:   f,
		idx: idx,
	}, nil
}

// Close releases resources.
func (tr *TabixReader) Close() error {
	if tr.f != nil {
		return tr.f.Close()
	}
	return nil
}

// Index returns the parsed TBI index, which contains column definitions,
// coordinate system info, and reference sequence names.
func (tr *TabixReader) Index() *BinIndex {
	return tr.idx
}

// Query returns an iterator that yields TabixRecords overlapping the
// 0-based half-open region [start, end) on the given reference.
func (tr *TabixReader) Query(ref string, start, end int) (*TabixIterator, error) {
	refID := tr.idx.RefID(ref)
	if refID < 0 {
		return nil, fmt.Errorf("tabix: unknown reference %q", ref)
	}

	chunks := tr.idx.Query(refID, start, end)
	if len(chunks) == 0 {
		return &TabixIterator{done: true}, nil
	}

	return &TabixIterator{
		ir:     tr.ir,
		idx:    tr.idx,
		chunks: chunks,
		ref:    ref,
		start:  start,
		end:    end,
	}, nil
}

// Next returns the next TabixRecord overlapping the query region.
// Returns nil, io.EOF when done.
func (ti *TabixIterator) Next() (*TabixRecord, error) {
	if ti.done {
		return nil, io.EOF
	}

	// Seek to the first chunk on first call.
	if !ti.started {
		ti.started = true
		if err := ti.ir.SeekToVirtualOffset(ti.chunks[0].Begin); err != nil {
			ti.done = true
			return nil, io.EOF
		}
		ti.scanner = bufio.NewScanner(ti.ir)
		ti.scanner.Buffer(make([]byte, 0, 1024*1024), 10*1024*1024)
	}

	// Read lines sequentially. The chunks are merged and sorted, so we
	// scan until we find an overlapping record, hit the end of all chunks,
	// or pass the query region.
	for ti.scanner.Scan() {
		line := ti.scanner.Text()

		// Skip empty lines and comment lines.
		if line == "" {
			continue
		}
		if ti.idx.Meta != 0 && line[0] == byte(ti.idx.Meta) {
			continue
		}

		// Parse the line to extract coordinates.
		rec, err := ti.parseLine(line)
		if err != nil {
			continue // skip unparseable lines
		}

		// Check reference match.
		if rec.Ref != ti.ref {
			// Past our reference — done.
			ti.done = true
			return nil, io.EOF
		}

		// Filter: record must overlap [start, end).
		if rec.End <= ti.start {
			continue // before query region
		}
		if rec.Start >= ti.end {
			// Past query region.
			ti.done = true
			return nil, io.EOF
		}

		return rec, nil
	}

	if err := ti.scanner.Err(); err != nil {
		return nil, err
	}

	ti.done = true
	return nil, io.EOF
}

// parseLine extracts the reference name and coordinates from a tab-delimited
// line using the column definitions from the TBI index.
func (ti *TabixIterator) parseLine(line string) (*TabixRecord, error) {
	fields := strings.Split(line, "\t")

	colSeq := int(ti.idx.ColSeq) - 1 // 1-based → 0-based column index
	colBeg := int(ti.idx.ColBeg) - 1
	colEnd := int(ti.idx.ColEnd) - 1

	if colSeq < 0 || colSeq >= len(fields) {
		return nil, fmt.Errorf("seq column %d out of range", colSeq)
	}
	if colBeg < 0 || colBeg >= len(fields) {
		return nil, fmt.Errorf("beg column %d out of range", colBeg)
	}

	ref := fields[colSeq]

	begStr := fields[colBeg]
	beg, err := strconv.Atoi(begStr)
	if err != nil {
		return nil, fmt.Errorf("parsing start: %w", err)
	}

	// Convert to 0-based if the file uses 1-based coordinates.
	if !ti.idx.ZeroBased {
		beg--
	}

	// End coordinate.
	end := beg + 1 // default: point feature
	if ti.idx.ColEnd != 0 && colEnd >= 0 && colEnd < len(fields) {
		endStr := fields[colEnd]
		e, err := strconv.Atoi(endStr)
		if err == nil {
			if !ti.idx.ZeroBased {
				// 1-based inclusive end → 0-based exclusive: no change needed
				// (1-based [1,10] = 0-based [0,10))
				end = e
			} else {
				end = e
			}
		}
	}

	return &TabixRecord{
		Line:  line,
		Ref:   ref,
		Start: beg,
		End:   end,
	}, nil
}

