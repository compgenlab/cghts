package bam

import (
	"container/heap"
	"fmt"
	"iter"
	"os"
	"sort"
	"sync"

	"github.com/compgen-io/cgltk/htsio"
)

const (
	defaultMaxMemory = 768 * 1024 * 1024
)

// sortedWriter buffers SamRecords in memory, sorts them, and writes
// sorted temp BAM files. On Close, temp files are merge-sorted into the
// final output. Supports coordinate sort and name sort.
type sortedWriter struct {
	filename  string
	header    *htsio.SamHeader
	refs      []bamRefInfo
	refIdx    map[string]int32
	sortCoord bool

	tmpPrefix string
	tmpFiles  []string
	tmpCount  int

	buf     []*htsio.SamRecord
	bufSize int
	maxMem  int

	mu     sync.Mutex
	closed bool
	err    error
}

// NewSortedWriter creates a sorted BAM writer. If sortCoord is true, records
// are coordinate-sorted; otherwise name-sorted.
func NewSortedWriter(filename string, header *htsio.SamHeader, sortCoord bool, tmpPrefix ...string) (*sortedWriter, error) {
	prefix := ""
	if len(tmpPrefix) > 0 {
		prefix = tmpPrefix[0]
	}
	return newSortedWriter(filename, header, sortCoord, prefix)
}

func newSortedWriter(filename string, header *htsio.SamHeader, sortCoord bool, tmpPrefix string) (*sortedWriter, error) {
	sw := &sortedWriter{
		filename:  filename,
		header:    header,
		sortCoord: sortCoord,
		tmpPrefix: tmpPrefix,
		maxMem:    defaultMaxMemory,
	}

	if sw.tmpPrefix == "" {
		sw.tmpPrefix = filename + ".tmp"
	}

	if header != nil {
		hrefs := header.References()
		sw.refs = make([]bamRefInfo, len(hrefs))
		sw.refIdx = make(map[string]int32, len(hrefs))
		for i, hr := range hrefs {
			sw.refs[i] = bamRefInfo{name: hr.Name, length: int32(hr.Length)}
			sw.refIdx[hr.Name] = int32(i)
		}
	} else {
		sw.refIdx = make(map[string]int32)
	}

	return sw, nil
}

func (sw *sortedWriter) Write(rec *htsio.SamRecord) error {
	sw.mu.Lock()
	defer sw.mu.Unlock()

	if sw.closed {
		return fmt.Errorf("bam: writer is closed")
	}
	if sw.err != nil {
		return sw.err
	}

	sw.buf = append(sw.buf, rec)
	sw.bufSize += estimateRecordSize(rec)

	if sw.bufSize >= sw.maxMem {
		if err := sw.flushBuffer(); err != nil {
			sw.err = err
			return err
		}
	}

	return nil
}

func (sw *sortedWriter) Close() error {
	sw.mu.Lock()
	defer sw.mu.Unlock()

	if sw.closed {
		return nil
	}
	sw.closed = true

	if sw.err != nil {
		sw.cleanup()
		return sw.err
	}

	if len(sw.tmpFiles) == 0 {
		sw.sortBuffer()
		if err := sw.writeBAM(sw.filename, sw.buf); err != nil {
			return err
		}
		return nil
	}

	if len(sw.buf) > 0 {
		if err := sw.flushBuffer(); err != nil {
			sw.cleanup()
			return err
		}
	}

	if err := sw.mergeFiles(); err != nil {
		sw.cleanup()
		return err
	}

	sw.cleanup()
	return nil
}

func (sw *sortedWriter) flushBuffer() error {
	if len(sw.buf) == 0 {
		return nil
	}

	sw.sortBuffer()

	sw.tmpCount++
	tmpPath := fmt.Sprintf("%s.%05d.bam", sw.tmpPrefix, sw.tmpCount)

	if err := sw.writeBAM(tmpPath, sw.buf); err != nil {
		return fmt.Errorf("writing temp file %s: %w", tmpPath, err)
	}

	sw.tmpFiles = append(sw.tmpFiles, tmpPath)
	sw.buf = sw.buf[:0]
	sw.bufSize = 0
	return nil
}

func (sw *sortedWriter) sortBuffer() {
	if sw.sortCoord {
		sort.Slice(sw.buf, func(i, j int) bool {
			return sw.lessCoord(sw.buf[i], sw.buf[j])
		})
	} else {
		sort.Slice(sw.buf, func(i, j int) bool {
			return sw.lessName(sw.buf[i], sw.buf[j])
		})
	}
}

func (sw *sortedWriter) lessCoord(a, b *htsio.SamRecord) bool {
	aRef := sw.refID(a)
	bRef := sw.refID(b)

	aUnmapped := aRef == -1
	bUnmapped := bRef == -1
	if aUnmapped != bUnmapped {
		return !aUnmapped
	}
	if aRef != bRef {
		return aRef < bRef
	}
	return a.Pos < b.Pos
}

func (sw *sortedWriter) lessName(a, b *htsio.SamRecord) bool {
	if a.ReadName != b.ReadName {
		return a.ReadName < b.ReadName
	}
	return a.Flag < b.Flag
}

func (sw *sortedWriter) refID(rec *htsio.SamRecord) int32 {
	if rec.RefName == "*" {
		return -1
	}
	if id, ok := sw.refIdx[rec.RefName]; ok {
		return id
	}
	return -1
}

func (sw *sortedWriter) writeBAM(filename string, records []*htsio.SamRecord) error {
	w, err := NewWriter(filename, sw.header)
	if err != nil {
		return err
	}
	for _, rec := range records {
		if err := w.Write(rec); err != nil {
			w.Close()
			return err
		}
	}
	return w.Close()
}

func (sw *sortedWriter) mergeFiles() error {
	readers := make([]*Reader, len(sw.tmpFiles))
	for i, path := range sw.tmpFiles {
		f, err := os.Open(path)
		if err != nil {
			for j := 0; j < i; j++ {
				readers[j].Close()
			}
			return fmt.Errorf("opening temp file %s: %w", path, err)
		}
		r, err := NewReader(f, path, nil)
		if err != nil {
			f.Close()
			for j := 0; j < i; j++ {
				readers[j].Close()
			}
			return fmt.Errorf("reading temp file %s: %w", path, err)
		}
		readers[i] = r
	}
	defer func() {
		for _, r := range readers {
			r.Close()
		}
	}()

	outWriter, err := NewWriter(sw.filename, sw.header)
	if err != nil {
		return err
	}

	type pullReader struct {
		next func() (*htsio.SamRecord, error, bool)
		stop func()
	}
	pullers := make([]pullReader, len(readers))
	for i, r := range readers {
		next, stop := iter.Pull2(r.Records())
		pullers[i] = pullReader{next: next, stop: stop}
	}
	defer func() {
		for _, p := range pullers {
			p.stop()
		}
	}()

	h := &mergeHeap{less: sw.lessCoord}
	if !sw.sortCoord {
		h.less = sw.lessName
	}

	for i := range pullers {
		rec, err, ok := pullers[i].next()
		if !ok {
			continue
		}
		if err != nil {
			outWriter.Close()
			return fmt.Errorf("reading from temp file: %w", err)
		}
		heap.Push(h, &mergeEntry{rec: rec, readerIdx: i})
	}

	for h.Len() > 0 {
		entry := heap.Pop(h).(*mergeEntry)
		if err := outWriter.Write(entry.rec); err != nil {
			outWriter.Close()
			return err
		}

		rec, err, ok := pullers[entry.readerIdx].next()
		if !ok {
			continue
		}
		if err != nil {
			outWriter.Close()
			return fmt.Errorf("reading from temp file: %w", err)
		}
		heap.Push(h, &mergeEntry{rec: rec, readerIdx: entry.readerIdx})
	}

	return outWriter.Close()
}

func (sw *sortedWriter) cleanup() {
	for _, path := range sw.tmpFiles {
		os.Remove(path)
	}
	sw.tmpFiles = nil
}

func estimateRecordSize(rec *htsio.SamRecord) int {
	size := 200
	size += len(rec.ReadName)
	size += len(rec.Cigar)
	size += len(rec.Seq)
	size += len(rec.Qual)
	size += len(rec.RefName)
	size += len(rec.RefNext)
	for k, v := range rec.Tags {
		size += len(k) + len(v.Value) + 16
	}
	return size
}

type mergeEntry struct {
	rec       *htsio.SamRecord
	readerIdx int
}

type mergeHeap struct {
	entries []*mergeEntry
	less    func(a, b *htsio.SamRecord) bool
}

func (h *mergeHeap) Len() int { return len(h.entries) }

func (h *mergeHeap) Less(i, j int) bool {
	return h.less(h.entries[i].rec, h.entries[j].rec)
}

func (h *mergeHeap) Swap(i, j int) {
	h.entries[i], h.entries[j] = h.entries[j], h.entries[i]
}

func (h *mergeHeap) Push(x interface{}) {
	h.entries = append(h.entries, x.(*mergeEntry))
}

func (h *mergeHeap) Pop() interface{} {
	old := h.entries
	n := len(old)
	entry := old[n-1]
	old[n-1] = nil
	h.entries = old[:n-1]
	return entry
}
