package tabix

import (
	"bufio"
	"fmt"
	"io"
	"iter"
	"os"
	"strconv"
	"strings"

	"github.com/compgenlab/cghts/htsio/bgzf"
)

// tabixIndex is the interface shared by TBI (BinIndex) and CSI (CSIIndex)
// for tabix region queries.
type tabixIndex interface {
	Query(refID int, start, end int) []Chunk
	RefID(name string) int
	RefNames() []string
	RecordCount(refID int) (uint64, bool)
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
	ir *bgzf.IndexedReader
	// closer releases whatever backs the reader: the file this package opened,
	// or a source the caller handed over with [WithCloser]. It is nil when the
	// caller keeps ownership.
	closer io.Closer
	idx    tabixIndex
	meta   tabixMeta

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
		meta = metaOfTBI(tbi)
	} else if _, err := os.Stat(csiPath); err == nil {
		csi, err := LoadCSI(csiPath)
		if err != nil {
			return nil, fmt.Errorf("tabix: loading CSI index: %w", err)
		}
		idx = csi
		meta = metaOfCSI(csi)
	} else {
		return nil, fmt.Errorf("tabix: no index found (.tbi or .csi) for %s", filename)
	}

	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}

	ir := bgzf.NewIndexedReader(f)

	return &Reader{
		ir:     ir,
		closer: f,
		idx:    idx,
		meta:   meta,
	}, nil
}

// ReaderOption configures a [Reader] built by [NewReaderFromSource].
type ReaderOption func(*readerOpts)

type readerOpts struct {
	cacheWindows int
	closer       io.Closer
}

// WithRecordCacheWindows sets the parsed-record cache size, in 16 kb windows.
// A value <= 0 disables the cache.
func WithRecordCacheWindows(n int) ReaderOption {
	return func(o *readerOpts) { o.cacheWindows = n }
}

// WithCloser hands ownership of a resource to the reader, so [Reader.Close]
// releases it.
//
// Ownership is explicit rather than inferred from whether the data stream
// happens to implement io.Closer: a remote source is usually wrapped in an
// io.SectionReader by the time it gets here, which is not a Closer even though
// the thing underneath it very much needs closing.
func WithCloser(c io.Closer) ReaderOption {
	return func(o *readerOpts) { o.closer = c }
}

// NewReaderFromSource opens a tabix-indexed stream that is not a local file.
//
// data must be seekable — BGZF decompression seeks to block boundaries, and the
// whole point of an index is to avoid reading the rest. index is the raw,
// still-compressed .tbi or .csi stream; which one it is is detected from its
// magic, so a caller that resolved the sidecar by trying both suffixes does not
// have to tell us which it found.
//
// For a remote source, wrap it with [iosource.ReadSeeker] and pass the source
// itself to [WithCloser].
func NewReaderFromSource(data io.ReadSeeker, index io.Reader, opts ...ReaderOption) (*Reader, error) {
	o := readerOpts{cacheWindows: DefaultRecordCacheWindows}
	for _, opt := range opts {
		opt(&o)
	}
	idx, meta, err := loadIndexFrom(index)
	if err != nil {
		return nil, err
	}
	tr := &Reader{
		ir:     bgzf.NewIndexedReader(data),
		closer: o.closer,
		idx:    idx,
		meta:   meta,
	}
	if o.cacheWindows > 0 {
		tr.recCache = newRecordCache(o.cacheWindows)
	}
	return tr, nil
}

// loadIndexFrom parses either index format, choosing by magic.
func loadIndexFrom(rd io.Reader) (tabixIndex, tabixMeta, error) {
	r := bgzf.NewReader(rd)

	var magic [4]byte
	if _, err := io.ReadFull(r, magic[:]); err != nil {
		return nil, tabixMeta{}, fmt.Errorf("tabix: reading index magic: %w", err)
	}
	switch magic {
	case [4]byte{'T', 'B', 'I', 1}:
		tbi, err := parseTBIBody(r)
		if err != nil {
			return nil, tabixMeta{}, fmt.Errorf("tabix: loading TBI index: %w", err)
		}
		return tbi, metaOfTBI(tbi), nil
	case [4]byte{'C', 'S', 'I', 1}:
		csi, err := parseCSIBody(r)
		if err != nil {
			return nil, tabixMeta{}, fmt.Errorf("tabix: loading CSI index: %w", err)
		}
		return csi, metaOfCSI(csi), nil
	default:
		return nil, tabixMeta{}, fmt.Errorf("tabix: not a TBI or CSI index (magic %x)", magic)
	}
}

func metaOfTBI(tbi *BinIndex) tabixMeta {
	return tabixMeta{
		Format: tbi.Format, ColSeq: tbi.ColSeq, ColBeg: tbi.ColBeg, ColEnd: tbi.ColEnd,
		Meta: tbi.Meta, Skip: tbi.Skip, ZeroBased: tbi.ZeroBased,
	}
}

func metaOfCSI(csi *CSIIndex) tabixMeta {
	return tabixMeta{
		Format: csi.Format, ColSeq: csi.ColSeq, ColBeg: csi.ColBeg, ColEnd: csi.ColEnd,
		Meta: csi.Meta, Skip: csi.Skip, ZeroBased: csi.ZeroBased,
	}
}

// DisableCache turns off the parsed-record cache (freeing any cached windows), so
// every subsequent query re-scans the chunk. Used by higher layers that cache parsed
// records themselves (e.g. vcf.IndexedVcfReader) and by tests comparing against the
// uncached path.
func (tr *Reader) DisableCache() { tr.recCache = nil }

// Close releases resources. It is a no-op when the caller retained ownership of
// the underlying source (see [WithCloser]).
func (tr *Reader) Close() error {
	if tr.closer != nil {
		return tr.closer.Close()
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

// RecordCount returns how many records the index says ref holds, without
// reading them; ok is false when the index does not record counts (see
// BinIndex.RecordCount) or does not know ref. Names match verbatim, as HasRef.
func (tr *Reader) RecordCount(ref string) (n uint64, ok bool) {
	return tr.idx.RecordCount(tr.idx.RefID(ref))
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

	if tr.meta.Skip <= 0 {
		tr.colNamesRead = true
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
		tr.colNamesRead = true
		return nil, nil
	}
	if tr.meta.Meta != 0 && header[0] == byte(tr.meta.Meta) {
		header = header[1:]
	}
	tr.colNames = strings.Split(header, "\t")
	// Marked read only now that it has been. Setting this up front cached a
	// failure as "no header": the first call surfaced the real I/O error and
	// every call after it returned (nil, nil), so a caller resolving several
	// columns in a loop saw a transport failure once and then "file has no
	// column header" for the rest.
	tr.colNamesRead = true
	return tr.colNames, nil
}

// HeaderLines returns the leading meta lines — those beginning with the index's
// meta character — from the start of the file.
//
// It reads through the same byte source as queries do, which is what lets a
// format whose header is meta-marked rather than counted (VCF, where every line
// up to and including #CHROM starts with '#') recover its header from a remote
// source instead of reopening the file by name.
//
// Returns nil when the index declares no meta character.
func (tr *Reader) HeaderLines() ([]string, error) {
	if tr.meta.Meta == 0 {
		return nil, nil
	}
	if err := tr.ir.SeekToVirtualOffset(0); err != nil {
		return nil, fmt.Errorf("tabix: reading header: %w", err)
	}
	sc := bufio.NewScanner(tr.ir)
	sc.Buffer(make([]byte, 0, 1024*1024), 10*1024*1024)

	var out []string
	for sc.Scan() {
		line := sc.Text()
		if line == "" || line[0] != byte(tr.meta.Meta) {
			break // first data line: the header is done
		}
		out = append(out, line)
	}
	if err := sc.Err(); err != nil {
		return nil, fmt.Errorf("tabix: reading header: %w", err)
	}
	return out, nil
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
	var out []*Record
	err := tr.scanChunks(chunks, func(line string) bool {
		if line == "" || (tr.meta.Meta != 0 && line[0] == byte(tr.meta.Meta)) {
			return true
		}
		rec, err := parseTabulatedLine(line, &tr.meta)
		if err != nil {
			return true
		}
		if rec.Ref != ref || rec.Start >= winEnd {
			return false
		}
		if rec.End > winStart {
			out = append(out, rec)
		}
		return true
	})
	if err != nil {
		return nil, err
	}
	return out, nil
}

func (tr *Reader) iterChunks(chunks []Chunk, ref string, start, end int) iter.Seq2[*Record, error] {
	return func(yield func(*Record, error) bool) {
		stopped := false
		err := tr.scanChunks(chunks, func(line string) bool {
			if line == "" || (tr.meta.Meta != 0 && line[0] == byte(tr.meta.Meta)) {
				return true
			}
			rec, err := parseTabulatedLine(line, &tr.meta)
			if err != nil {
				return true
			}
			if rec.Ref != ref || rec.Start >= end {
				return false
			}
			if rec.End <= start {
				return true
			}
			if !yield(rec, nil) {
				stopped = true
				return false
			}
			return true
		})
		if err != nil && !stopped {
			yield(nil, err)
		}
	}
}

// scanChunks hands fn each line in the index's chunks, in file order, until
// fn returns false or the chunks run out.
//
// It reads each chunk from its start to its end and jumps the gaps between
// them, which is how a tabix index is meant to be read. It used to seek to the
// first chunk and read every line from there until the region was passed. A
// region's first chunk is often a top-level bin's: in a GTF the long genes
// live there, so it starts far upstream of the region, and every query read
// -- and decompressed and parsed -- megabytes it then threw away.
//
// It seeks only when behind a chunk's start, never back to one: a reader at
// the end of a block is at the same file position as the next block's start,
// so the last line read for one chunk can be the first of the next, and
// seeking back to it would hand fn that line twice.
func (tr *Reader) scanChunks(chunks []Chunk, fn func(line string) bool) error {
	var buf []byte
	for i, ch := range chunks {
		if i == 0 || tr.ir.VirtualTell() < ch.Begin {
			if err := tr.ir.SeekToVirtualOffset(ch.Begin); err != nil {
				return fmt.Errorf("tabix: seeking to chunk: %w", err)
			}
		}
		for tr.ir.VirtualTell() < ch.End {
			line, err := readLine(tr.ir, buf[:0])
			buf = line
			if len(line) > 0 {
				if !fn(string(line)) {
					return nil
				}
			}
			if err == io.EOF {
				return nil
			}
			if err != nil {
				return err
			}
		}
	}
	return nil
}

// readLine appends one line from r to buf, without its newline (or a trailing
// carriage return). It returns io.EOF only once nothing is left; a final line
// with no newline is returned with a nil error first.
func readLine(r io.ByteReader, buf []byte) ([]byte, error) {
	for {
		c, err := r.ReadByte()
		if err != nil {
			if err == io.EOF && len(buf) > 0 {
				return trimCR(buf), nil
			}
			return trimCR(buf), err
		}
		if c == '\n' {
			return trimCR(buf), nil
		}
		buf = append(buf, c)
	}
}

func trimCR(b []byte) []byte {
	if n := len(b); n > 0 && b[n-1] == '\r' {
		return b[:n-1]
	}
	return b
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
	switch {
	case meta.Format&presetMask == PresetVCF:
		// Must mirror the writer: this end feeds the overlap test, so a record
		// understated here is dropped even when the index pointed straight at
		// it.
		end = vcfSpanEnd(fields, beg)
	case meta.ColEnd != 0 && colEnd >= 0 && colEnd < len(fields):
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
