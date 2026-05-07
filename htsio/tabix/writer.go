package tabix

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/compgen-io/cgkit/htsio/bgzf"
)

const (
	defaultMaxMemory = 768 * 1024 * 1024
)

// WriterOpts configures a Writer.
type WriterOpts struct {
	colSeq    int32
	colBeg    int32
	colEnd    int32
	meta      int32
	skip      int32
	zeroBased bool
	autoIndex bool
}

// NewWriterOpts creates a WriterOpts with default values.
func NewWriterOpts() *WriterOpts {
	return &WriterOpts{
		colSeq: 1,
		colBeg: 2,
		colEnd: 3,
	}
}

func (o *WriterOpts) Columns(seq, beg, end int) *WriterOpts {
	o.colSeq = int32(seq)
	o.colBeg = int32(beg)
	o.colEnd = int32(end)
	return o
}

func (o *WriterOpts) Meta(ch byte) *WriterOpts {
	o.meta = int32(ch)
	return o
}

func (o *WriterOpts) Skip(n int) *WriterOpts {
	o.skip = int32(n)
	return o
}

func (o *WriterOpts) ZeroBased() *WriterOpts {
	o.zeroBased = true
	return o
}

func (o *WriterOpts) AutoIndex() *WriterOpts {
	o.autoIndex = true
	return o
}

func (o *WriterOpts) BED() *WriterOpts {
	o.colSeq = 1
	o.colBeg = 2
	o.colEnd = 3
	o.meta = 0
	o.zeroBased = true
	return o
}

func (o *WriterOpts) VCF() *WriterOpts {
	o.colSeq = 1
	o.colBeg = 2
	o.colEnd = 0
	o.meta = '#'
	o.zeroBased = false
	return o
}

func (o *WriterOpts) GFF() *WriterOpts {
	o.colSeq = 1
	o.colBeg = 4
	o.colEnd = 5
	o.meta = '#'
	o.zeroBased = false
	return o
}

type tabixLine struct {
	line  string
	ref   string
	start int
}

// Writer writes sorted, BGZF-compressed tabular files with optional
// .tbi index generation.
type Writer struct {
	filename string
	opts     *WriterOpts

	headerLines []string

	buf     []tabixLine
	bufSize int
	maxMem  int

	tmpPrefix string
	tmpFiles  []string
	tmpCount  int

	refOrder []string
	refIdx   map[string]int

	mu     sync.Mutex
	closed bool
	err    error
}

// NewWriter creates a sorted BGZF tabular writer.
func NewWriter(filename string, opts *WriterOpts) *Writer {
	return &Writer{
		filename: filename,
		opts:     opts,
		maxMem:   defaultMaxMemory,
		refIdx:   make(map[string]int),
	}
}

func (tw *Writer) WriteHeader(line string) {
	tw.headerLines = append(tw.headerLines, line)
}

func (tw *Writer) Write(line string) error {
	tw.mu.Lock()
	defer tw.mu.Unlock()

	if tw.closed {
		return fmt.Errorf("tabix: writer is closed")
	}
	if tw.err != nil {
		return tw.err
	}

	parsed, err := tw.parseLine(line)
	if err != nil {
		return err
	}

	if _, ok := tw.refIdx[parsed.ref]; !ok {
		tw.refIdx[parsed.ref] = len(tw.refOrder)
		tw.refOrder = append(tw.refOrder, parsed.ref)
	}

	tw.buf = append(tw.buf, parsed)
	tw.bufSize += len(line) + 64

	if tw.bufSize >= tw.maxMem {
		if err := tw.flushBuffer(); err != nil {
			tw.err = err
			return err
		}
	}

	return nil
}

func (tw *Writer) Close() error {
	tw.mu.Lock()
	defer tw.mu.Unlock()

	if tw.closed {
		return nil
	}
	tw.closed = true

	if tw.err != nil {
		tw.cleanup()
		return tw.err
	}

	if len(tw.tmpFiles) == 0 {
		tw.sortBuffer()
		return tw.writeFinal(tw.buf)
	}

	if len(tw.buf) > 0 {
		if err := tw.flushBuffer(); err != nil {
			tw.cleanup()
			return err
		}
	}

	if err := tw.mergeFiles(); err != nil {
		tw.cleanup()
		return err
	}

	tw.cleanup()
	return nil
}

func (tw *Writer) parseLine(line string) (tabixLine, error) {
	fields := strings.Split(line, "\t")
	colSeq := int(tw.opts.colSeq) - 1
	colBeg := int(tw.opts.colBeg) - 1

	if colSeq < 0 || colSeq >= len(fields) {
		return tabixLine{}, fmt.Errorf("tabix: seq column %d out of range for line with %d fields", colSeq+1, len(fields))
	}
	if colBeg < 0 || colBeg >= len(fields) {
		return tabixLine{}, fmt.Errorf("tabix: beg column %d out of range for line with %d fields", colBeg+1, len(fields))
	}

	ref := fields[colSeq]
	beg, err := strconv.Atoi(fields[colBeg])
	if err != nil {
		return tabixLine{}, fmt.Errorf("tabix: parsing start: %w", err)
	}

	if !tw.opts.zeroBased {
		beg--
	}

	return tabixLine{line: line, ref: ref, start: beg}, nil
}

func (tw *Writer) sortBuffer() {
	sort.Slice(tw.buf, func(i, j int) bool {
		return tw.less(tw.buf[i], tw.buf[j])
	})
}

func (tw *Writer) less(a, b tabixLine) bool {
	ai, aok := tw.refIdx[a.ref]
	bi, bok := tw.refIdx[b.ref]
	if !aok {
		ai = len(tw.refOrder)
	}
	if !bok {
		bi = len(tw.refOrder)
	}
	if ai != bi {
		return ai < bi
	}
	return a.start < b.start
}

func (tw *Writer) flushBuffer() error {
	if len(tw.buf) == 0 {
		return nil
	}

	tw.sortBuffer()

	tw.tmpCount++
	prefix := tw.tmpPrefix
	if prefix == "" {
		prefix = tw.filename + ".tmp"
	}
	tmpPath := fmt.Sprintf("%s.%05d.gz", prefix, tw.tmpCount)

	if err := tw.writeBGZFSimple(tmpPath, tw.buf); err != nil {
		return fmt.Errorf("writing temp file %s: %w", tmpPath, err)
	}

	tw.tmpFiles = append(tw.tmpFiles, tmpPath)
	tw.buf = tw.buf[:0]
	tw.bufSize = 0
	return nil
}

func (tw *Writer) writeBGZFSimple(filename string, lines []tabixLine) error {
	w, err := bgzf.NewBGZipFile(filename)
	if err != nil {
		return err
	}

	for _, l := range lines {
		if _, err := io.WriteString(w, l.line+"\n"); err != nil {
			w.Close()
			return err
		}
	}

	return w.Close()
}

func (tw *Writer) writeFinal(lines []tabixLine) error {
	w, err := bgzf.NewBGZipFile(tw.filename)
	if err != nil {
		return err
	}

	for _, h := range tw.headerLines {
		if _, err := io.WriteString(w, h+"\n"); err != nil {
			w.Close()
			return err
		}
	}

	var ib *tbiIndexBuilder
	if tw.opts.autoIndex {
		ib = newTBIIndexBuilder(tw.opts, tw.refOrder)
	}

	for _, l := range lines {
		if ib != nil {
			ib.addRecord(l, w.VirtualTell())
		}
		if _, err := io.WriteString(w, l.line+"\n"); err != nil {
			w.Close()
			return err
		}
	}

	if err := w.Close(); err != nil {
		return err
	}

	if ib != nil {
		return ib.writeTBI(tw.filename + ".tbi")
	}
	return nil
}

func (tw *Writer) mergeFiles() error {
	type lineReader struct {
		r       *bgzf.Reader
		f       *os.File
		scanner *bufio.Scanner
	}

	readers := make([]*lineReader, len(tw.tmpFiles))
	for i, path := range tw.tmpFiles {
		f, err := os.Open(path)
		if err != nil {
			for j := 0; j < i; j++ {
				readers[j].f.Close()
			}
			return fmt.Errorf("opening temp file %s: %w", path, err)
		}
		r := bgzf.NewReader(f)
		s := bufio.NewScanner(r)
		s.Buffer(make([]byte, 0, 1024*1024), 10*1024*1024)
		readers[i] = &lineReader{r: r, f: f, scanner: s}
	}
	defer func() {
		for _, r := range readers {
			r.f.Close()
		}
	}()

	outFile, err := bgzf.NewBGZipFile(tw.filename)
	if err != nil {
		return err
	}

	for _, h := range tw.headerLines {
		if _, err := io.WriteString(outFile, h+"\n"); err != nil {
			outFile.Close()
			return err
		}
	}

	var ib *tbiIndexBuilder
	if tw.opts.autoIndex {
		ib = newTBIIndexBuilder(tw.opts, tw.refOrder)
	}

	h := &tabixMergeHeap{less: tw.less}

	for i, r := range readers {
		if r.scanner.Scan() {
			parsed, err := tw.parseLine(r.scanner.Text())
			if err != nil {
				outFile.Close()
				return err
			}
			h.push(&tabixMergeItem{line: parsed, readerIdx: i})
		}
	}

	for h.Len() > 0 {
		item := h.pop()
		if ib != nil {
			ib.addRecord(item.line, outFile.VirtualTell())
		}
		if _, err := io.WriteString(outFile, item.line.line+"\n"); err != nil {
			outFile.Close()
			return err
		}

		r := readers[item.readerIdx]
		if r.scanner.Scan() {
			parsed, err := tw.parseLine(r.scanner.Text())
			if err != nil {
				outFile.Close()
				return err
			}
			h.push(&tabixMergeItem{line: parsed, readerIdx: item.readerIdx})
		}
	}

	if err := outFile.Close(); err != nil {
		return err
	}

	if ib != nil {
		return ib.writeTBI(tw.filename + ".tbi")
	}
	return nil
}

type tbiRefBuilder struct {
	bins      map[uint32][]Chunk
	linearIdx map[int]bgzf.VirtualOffset
	lastVO    bgzf.VirtualOffset
}

type tbiIndexBuilder struct {
	opts     *WriterOpts
	refOrder []string
	refIdx   map[string]int
	refs     map[string]*tbiRefBuilder
}

func newTBIIndexBuilder(opts *WriterOpts, refOrder []string) *tbiIndexBuilder {
	refIdx := make(map[string]int, len(refOrder))
	for i, name := range refOrder {
		refIdx[name] = i
	}
	return &tbiIndexBuilder{
		opts:     opts,
		refOrder: refOrder,
		refIdx:   refIdx,
		refs:     make(map[string]*tbiRefBuilder),
	}
}

func (ib *tbiIndexBuilder) addRecord(l tabixLine, vo bgzf.VirtualOffset) {
	rb, ok := ib.refs[l.ref]
	if !ok {
		rb = &tbiRefBuilder{
			bins:      make(map[uint32][]Chunk),
			linearIdx: make(map[int]bgzf.VirtualOffset),
		}
		ib.refs[l.ref] = rb
	}

	end := l.start + 1
	fields := strings.Split(l.line, "\t")
	colEnd := int(ib.opts.colEnd) - 1
	if ib.opts.colEnd != 0 && colEnd >= 0 && colEnd < len(fields) {
		if e, err := strconv.Atoi(fields[colEnd]); err == nil {
			end = e
		}
	}

	bin := Reg2Bin(l.start, end)

	chunks := rb.bins[uint32(bin)]
	if len(chunks) > 0 {
		last := &chunks[len(chunks)-1]
		if vo.BlockOffset() == last.End.BlockOffset() || vo == last.End {
			last.End = vo
		} else {
			chunks = append(chunks, Chunk{Begin: vo, End: vo})
		}
	} else {
		chunks = append(chunks, Chunk{Begin: vo, End: vo})
	}
	rb.bins[uint32(bin)] = chunks
	rb.lastVO = vo

	window := l.start >> 14
	if existing, ok := rb.linearIdx[window]; !ok || vo < existing {
		rb.linearIdx[window] = vo
	}
}

func (ib *tbiIndexBuilder) writeTBI(path string) error {
	tbiFile, err := os.Create(path)
	if err != nil {
		return err
	}

	w := bgzf.NewWriter(tbiFile)

	w.Write([]byte("TBI\x01"))
	binary.Write(w, binary.LittleEndian, int32(len(ib.refOrder)))

	fmtVal := int32(0)
	if ib.opts.zeroBased {
		fmtVal |= 0x10000
	}
	binary.Write(w, binary.LittleEndian, fmtVal)
	binary.Write(w, binary.LittleEndian, ib.opts.colSeq)
	binary.Write(w, binary.LittleEndian, ib.opts.colBeg)
	binary.Write(w, binary.LittleEndian, ib.opts.colEnd)
	binary.Write(w, binary.LittleEndian, ib.opts.meta)
	binary.Write(w, binary.LittleEndian, ib.opts.skip)

	var namesData []byte
	for _, name := range ib.refOrder {
		namesData = append(namesData, []byte(name)...)
		namesData = append(namesData, 0)
	}
	binary.Write(w, binary.LittleEndian, int32(len(namesData)))
	w.Write(namesData)

	for _, refName := range ib.refOrder {
		rb, ok := ib.refs[refName]
		if !ok {
			binary.Write(w, binary.LittleEndian, int32(0))
			binary.Write(w, binary.LittleEndian, int32(0))
			continue
		}

		mergedBins := make(map[uint32][]Chunk)
		for bin, chunks := range rb.bins {
			if len(chunks) == 0 {
				continue
			}
			sort.Slice(chunks, func(i, j int) bool {
				return chunks[i].Begin < chunks[j].Begin
			})
			merged := []Chunk{chunks[0]}
			for i := 1; i < len(chunks); i++ {
				last := &merged[len(merged)-1]
				if chunks[i].Begin <= last.End {
					if chunks[i].End > last.End {
						last.End = chunks[i].End
					}
				} else {
					merged = append(merged, chunks[i])
				}
			}
			if rb.lastVO > merged[len(merged)-1].End {
				merged[len(merged)-1].End = rb.lastVO
			}
			mergedBins[bin] = merged
		}

		binary.Write(w, binary.LittleEndian, int32(len(mergedBins)))
		for bin, chunks := range mergedBins {
			binary.Write(w, binary.LittleEndian, bin)
			binary.Write(w, binary.LittleEndian, int32(len(chunks)))
			for _, c := range chunks {
				binary.Write(w, binary.LittleEndian, uint64(c.Begin))
				binary.Write(w, binary.LittleEndian, uint64(c.End))
			}
		}

		maxWindow := 0
		for wnd := range rb.linearIdx {
			if wnd > maxWindow {
				maxWindow = wnd
			}
		}
		nIntervals := int32(maxWindow + 1)
		binary.Write(w, binary.LittleEndian, nIntervals)
		for i := int32(0); i < nIntervals; i++ {
			if vo, ok := rb.linearIdx[int(i)]; ok {
				binary.Write(w, binary.LittleEndian, uint64(vo))
			} else {
				binary.Write(w, binary.LittleEndian, uint64(0))
			}
		}
	}

	if err := w.Close(); err != nil {
		tbiFile.Close()
		return err
	}
	return tbiFile.Close()
}

func (tw *Writer) cleanup() {
	for _, path := range tw.tmpFiles {
		os.Remove(path)
	}
	tw.tmpFiles = nil
}

type tabixMergeItem struct {
	line      tabixLine
	readerIdx int
}

type tabixMergeHeap struct {
	items []*tabixMergeItem
	less  func(a, b tabixLine) bool
}

func (h *tabixMergeHeap) Len() int { return len(h.items) }

func (h *tabixMergeHeap) push(item *tabixMergeItem) {
	h.items = append(h.items, item)
	i := len(h.items) - 1
	for i > 0 {
		parent := (i - 1) / 2
		if !h.less(h.items[i].line, h.items[parent].line) {
			break
		}
		h.items[i], h.items[parent] = h.items[parent], h.items[i]
		i = parent
	}
}

func (h *tabixMergeHeap) pop() *tabixMergeItem {
	n := len(h.items)
	result := h.items[0]
	h.items[0] = h.items[n-1]
	h.items[n-1] = nil
	h.items = h.items[:n-1]
	i := 0
	for {
		left := 2*i + 1
		if left >= len(h.items) {
			break
		}
		smallest := left
		right := left + 1
		if right < len(h.items) && h.less(h.items[right].line, h.items[left].line) {
			smallest = right
		}
		if !h.less(h.items[smallest].line, h.items[i].line) {
			break
		}
		h.items[i], h.items[smallest] = h.items[smallest], h.items[i]
		i = smallest
	}
	return result
}
