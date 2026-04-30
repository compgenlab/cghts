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

// tabixIndex is the interface shared by TBI (BinIndex) and CSI (CSIIndex)
// for tabix region queries.
type tabixIndex interface {
	Query(refID int, start, end int) []Chunk
	RefID(name string) int
}

// tabixMeta holds the column definitions and coordinate metadata from
// a tabix index (TBI or CSI).
type tabixMeta struct {
	Format    int32
	ColSeq    int32
	ColBeg    int32
	ColEnd    int32
	Meta      int32
	Skip      int32
	ZeroBased bool
}

// TabixReader reads BGZF-compressed, tabix-indexed text files (BED, VCF, GFF,
// etc.) with random access by genomic region. Supports both .tbi and .csi
// indexes, which provide column definitions, the coordinate system (0-based
// vs 1-based), comment character, and header skip count.
type TabixReader struct {
	ir   *bgzf.IndexedReader
	f    *os.File
	idx  tabixIndex
	meta tabixMeta
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
	meta    *tabixMeta
	scanner *bufio.Scanner
	chunks  []Chunk
	ref     string
	start   int // query start, 0-based
	end     int // query end, 0-based exclusive
	started bool
	done    bool
}

// NewTabixReader opens a BGZF-compressed file and its tabix index.
// It looks for a .tbi index first, then falls back to .csi.
func NewTabixReader(filename string) (*TabixReader, error) {
	var idx tabixIndex
	var meta tabixMeta

	tbiPath := filename + ".tbi"
	csiPath := filename + ".csi"

	if _, err := os.Stat(tbiPath); err == nil {
		tbi, err := LoadTBI(tbiPath)
		if err != nil {
			return nil, fmt.Errorf("tabix: loading TBI index: %w", err)
		}
		idx = tbi
		meta = tabixMeta{
			Format:    tbi.Format,
			ColSeq:    tbi.ColSeq,
			ColBeg:    tbi.ColBeg,
			ColEnd:    tbi.ColEnd,
			Meta:      tbi.Meta,
			Skip:      tbi.Skip,
			ZeroBased: tbi.ZeroBased,
		}
	} else if _, err := os.Stat(csiPath); err == nil {
		csi, err := LoadCSI(csiPath)
		if err != nil {
			return nil, fmt.Errorf("tabix: loading CSI index: %w", err)
		}
		idx = csi
		meta = tabixMeta{
			Format:    csi.Format,
			ColSeq:    csi.ColSeq,
			ColBeg:    csi.ColBeg,
			ColEnd:    csi.ColEnd,
			Meta:      csi.Meta,
			Skip:      csi.Skip,
			ZeroBased: csi.ZeroBased,
		}
	} else {
		return nil, fmt.Errorf("tabix: no index found (.tbi or .csi) for %s", filename)
	}

	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}

	ir := bgzf.NewIndexedReader(f)

	return &TabixReader{
		ir:   ir,
		f:    f,
		idx:  idx,
		meta: meta,
	}, nil
}

// Close releases resources.
func (tr *TabixReader) Close() error {
	if tr.f != nil {
		return tr.f.Close()
	}
	return nil
}

// Meta returns the tabix metadata (column definitions, coordinate system).
func (tr *TabixReader) Meta() tabixMeta {
	return tr.meta
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
		meta:   &tr.meta,
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
		if ti.meta.Meta != 0 && line[0] == byte(ti.meta.Meta) {
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

	colSeq := int(ti.meta.ColSeq) - 1 // 1-based → 0-based column index
	colBeg := int(ti.meta.ColBeg) - 1
	colEnd := int(ti.meta.ColEnd) - 1

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
	if !ti.meta.ZeroBased {
		beg--
	}

	// End coordinate.
	end := beg + 1 // default: point feature
	if ti.meta.ColEnd != 0 && colEnd >= 0 && colEnd < len(fields) {
		endStr := fields[colEnd]
		e, err := strconv.Atoi(endStr)
		if err == nil {
			if !ti.meta.ZeroBased {
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

