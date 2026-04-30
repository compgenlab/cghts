package htsio

import (
	"bufio"
	"fmt"
	"iter"
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

// Query returns an iterator over TabixRecords overlapping the 0-based
// half-open region [start, end) on the given reference.
func (tr *TabixReader) Query(ref string, start, end int) (iter.Seq2[*TabixRecord, error], error) {
	refID := tr.idx.RefID(ref)
	if refID < 0 {
		return nil, fmt.Errorf("tabix: unknown reference %q", ref)
	}

	chunks := tr.idx.Query(refID, start, end)
	if len(chunks) == 0 {
		return func(yield func(*TabixRecord, error) bool) {}, nil
	}

	return tr.iterChunks(chunks, ref, start, end), nil
}

// iterChunks returns an iterator that reads lines from the given chunks,
// parses coordinates, and yields records overlapping [start, end).
func (tr *TabixReader) iterChunks(chunks []Chunk, ref string, start, end int) iter.Seq2[*TabixRecord, error] {
	return func(yield func(*TabixRecord, error) bool) {
		if err := tr.ir.SeekToVirtualOffset(chunks[0].Begin); err != nil {
			yield(nil, fmt.Errorf("tabix: seeking to chunk: %w", err))
			return
		}

		scanner := bufio.NewScanner(tr.ir)
		scanner.Buffer(make([]byte, 0, 1024*1024), 10*1024*1024)

		for scanner.Scan() {
			line := scanner.Text()

			if line == "" {
				continue
			}
			if tr.meta.Meta != 0 && line[0] == byte(tr.meta.Meta) {
				continue
			}

			rec, err := parseTabulatedLine(line, &tr.meta)
			if err != nil {
				continue
			}

			if rec.Ref != ref {
				return // past our reference
			}
			if rec.End <= start {
				continue
			}
			if rec.Start >= end {
				return
			}

			if !yield(rec, nil) {
				return
			}
		}

		if err := scanner.Err(); err != nil {
			yield(nil, err)
		}
	}
}

// parseTabulatedLine extracts the reference name and coordinates from a
// tab-delimited line using the column definitions from the tabix metadata.
func parseTabulatedLine(line string, meta *tabixMeta) (*TabixRecord, error) {
	fields := strings.Split(line, "\t")

	colSeq := int(meta.ColSeq) - 1 // 1-based → 0-based column index
	colBeg := int(meta.ColBeg) - 1
	colEnd := int(meta.ColEnd) - 1

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
	if !meta.ZeroBased {
		beg--
	}

	// End coordinate.
	end := beg + 1 // default: point feature
	if meta.ColEnd != 0 && colEnd >= 0 && colEnd < len(fields) {
		endStr := fields[colEnd]
		e, err := strconv.Atoi(endStr)
		if err == nil {
			end = e
		}
	}

	return &TabixRecord{
		Line:  line,
		Ref:   ref,
		Start: beg,
		End:   end,
	}, nil
}

