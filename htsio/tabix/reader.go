package tabix

import (
	"bufio"
	"fmt"
	"iter"
	"os"
	"strconv"
	"strings"

	"github.com/compgenlab/hts/htsio/bgzf"
)

// tabixIndex is the interface shared by TBI (BinIndex) and CSI (CSIIndex)
// for tabix region queries.
type tabixIndex interface {
	Query(refID int, start, end int) []Chunk
	RefID(name string) int
	RefNames() []string
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

// Reader reads BGZF-compressed, tabix-indexed text files (BED, VCF, GFF,
// etc.) with random access by genomic region.
type Reader struct {
	ir   *bgzf.IndexedReader
	f    *os.File
	idx  tabixIndex
	meta tabixMeta

	colNames     []string
	colNamesRead bool

	recCache *recordCache // parsed-record LRU; nil = disabled
}

// Record holds a single parsed line from a tabix query along with the
// extracted genomic coordinates.
type Record struct {
	Line  string
	Ref   string
	Start int // 0-based
	End   int // 0-based, exclusive
}

// NewReader opens a BGZF-compressed file and its tabix index.
// It looks for a .tbi index first, then falls back to .csi.
//
// The returned reader caches parsed records per 16 kb window
// ([DefaultRecordCacheWindows] windows) so a position-sorted stream of queries does
// not re-scan/re-parse the same chunk; use [NewReaderSize] to size or disable that
// cache. The reader is not safe for concurrent use.
func NewReader(filename string) (*Reader, error) {
	return NewReaderSize(filename, DefaultRecordCacheWindows)
}

// NewReaderSize is [NewReader] with an explicit parsed-record cache size, in 16 kb
// windows. A cacheWindows <= 0 disables the record cache (every query re-scans the
// chunk, the pre-cache behavior).
func NewReaderSize(filename string, cacheWindows int) (*Reader, error) {
	tr, err := openReader(filename)
	if err != nil {
		return nil, err
	}
	if cacheWindows > 0 {
		tr.recCache = newRecordCache(cacheWindows)
	}
	return tr, nil
}

func openReader(filename string) (*Reader, error) {
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

	return &Reader{
		ir:   ir,
		f:    f,
		idx:  idx,
		meta: meta,
	}, nil
}

// DisableCache turns off the parsed-record cache (freeing any cached windows), so
// every subsequent query re-scans the chunk. Used by higher layers that cache parsed
// records themselves (e.g. vcf.IndexedVcfReader) and by tests comparing against the
// uncached path.
func (tr *Reader) DisableCache() { tr.recCache = nil }

// Close releases resources.
func (tr *Reader) Close() error {
	if tr.f != nil {
		return tr.f.Close()
	}
	return nil
}

// Meta returns the tabix metadata (column definitions, coordinate system).
func (tr *Reader) Meta() tabixMeta {
	return tr.meta
}

// HasRef reports whether the index contains the given reference name.
func (tr *Reader) HasRef(ref string) bool {
	return tr.idx.RefID(ref) >= 0
}

// RefNames returns the reference sequence names present in the index, in
// reference order. It is the contig list used to build a contig-name converter.
func (tr *Reader) RefNames() []string {
	return tr.idx.RefNames()
}

// ColumnNames returns the column names from the file's header line, with a
// leading meta character (e.g. '#') stripped. The header is identified the
// tabix way: the last of the index's skipped lines. A file with no skipped
// lines (Skip == 0) has no column header and this returns nil. The result is
// cached.
func (tr *Reader) ColumnNames() ([]string, error) {
	if tr.colNamesRead {
		return tr.colNames, nil
	}
	tr.colNamesRead = true

	if tr.meta.Skip <= 0 {
		return nil, nil // no skipped line => no header
	}
	if err := tr.ir.SeekToVirtualOffset(0); err != nil {
		return nil, fmt.Errorf("tabix: reading header: %w", err)
	}
	sc := bufio.NewScanner(tr.ir)
	sc.Buffer(make([]byte, 0, 1024*1024), 10*1024*1024)

	header := ""
	for n := 0; n < int(tr.meta.Skip); n++ {
		if !sc.Scan() {
			break
		}
		header = sc.Text()
	}
	if err := sc.Err(); err != nil {
		return nil, fmt.Errorf("tabix: reading header: %w", err)
	}
	if header == "" {
		return nil, nil
	}
	if tr.meta.Meta != 0 && header[0] == byte(tr.meta.Meta) {
		header = header[1:]
	}
	tr.colNames = strings.Split(header, "\t")
	return tr.colNames, nil
}

// ColumnByName returns the 1-based column number for the named column, matching
// it against the header (see [Reader.ColumnNames]). It returns an error when the
// file has no header or the name is not found.
func (tr *Reader) ColumnByName(name string) (int, error) {
	names, err := tr.ColumnNames()
	if err != nil {
		return 0, err
	}
	if len(names) == 0 {
		return 0, fmt.Errorf("tabix: file has no column header to resolve %q", name)
	}
	for i, n := range names {
		if n == name {
			return i + 1, nil
		}
	}
	return 0, fmt.Errorf("tabix: column %q not found in header", name)
}

// Query returns an iterator over Records overlapping the 0-based
// half-open region [start, end) on the given reference.
func (tr *Reader) Query(ref string, start, end int) (iter.Seq2[*Record, error], error) {
	refID := tr.idx.RefID(ref)
	if refID < 0 {
		return nil, fmt.Errorf("tabix: unknown reference %q", ref)
	}

	// Cache path: for a query contained in a single 16 kb window (the common
	// point/near-point case) serve from — or hydrate — the record cache. A window's
	// hydrated records cover every sub-query within it, so no coverage guard is
	// needed. Wider queries (start<0, empty, or spanning windows) fall through to the
	// direct chunk scan and are not cached.
	if tr.recCache != nil && start >= 0 && start < end && (start>>windowShift) == ((end-1)>>windowShift) {
		key := windowKey{refID: refID, win: start >> windowShift}
		recs, ok := tr.recCache.get(key)
		if !ok {
			winStart := key.win << windowShift
			var err error
			recs, err = tr.hydrateWindow(refID, ref, winStart, winStart+(1<<windowShift))
			if err != nil {
				return nil, err
			}
			tr.recCache.put(key, recs)
		}
		return sliceSeq(recs, start, end), nil
	}

	chunks := tr.idx.Query(refID, start, end)
	if len(chunks) == 0 {
		return func(yield func(*Record, error) bool) {}, nil
	}

	return tr.iterChunks(chunks, ref, start, end), nil
}

// hydrateWindow scans the full 16 kb window [winStart, winEnd) once and returns
// every record overlapping it, in file order. It reuses the same seek + scanner +
// parse loop as iterChunks, so records straddling a bgzf block boundary reassemble
// transparently; only the stop bound (the window, not a caller's sub-region) differs.
func (tr *Reader) hydrateWindow(refID int, ref string, winStart, winEnd int) ([]*Record, error) {
	chunks := tr.idx.Query(refID, winStart, winEnd)
	if len(chunks) == 0 {
		return nil, nil
	}
	if err := tr.ir.SeekToVirtualOffset(chunks[0].Begin); err != nil {
		return nil, fmt.Errorf("tabix: seeking to chunk: %w", err)
	}

	scanner := bufio.NewScanner(tr.ir)
	scanner.Buffer(make([]byte, 0, 1024*1024), 10*1024*1024)

	var out []*Record
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
			break
		}
		if rec.End <= winStart {
			continue
		}
		if rec.Start >= winEnd {
			break
		}
		out = append(out, rec)
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return out, nil
}

func (tr *Reader) iterChunks(chunks []Chunk, ref string, start, end int) iter.Seq2[*Record, error] {
	return func(yield func(*Record, error) bool) {
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
				return
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

func parseTabulatedLine(line string, meta *tabixMeta) (*Record, error) {
	fields := strings.Split(line, "\t")

	colSeq := int(meta.ColSeq) - 1
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

	if !meta.ZeroBased {
		beg--
	}

	end := beg + 1
	if meta.ColEnd != 0 && colEnd >= 0 && colEnd < len(fields) {
		endStr := fields[colEnd]
		e, err := strconv.Atoi(endStr)
		if err == nil {
			end = e
		}
	}

	return &Record{
		Line:  line,
		Ref:   ref,
		Start: beg,
		End:   end,
	}, nil
}
