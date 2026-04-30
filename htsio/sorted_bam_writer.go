package htsio

import (
	"container/heap"
	"fmt"
	"io"
	"os"
	"sort"
	"sync"
)

const (
	// defaultMaxMemory is the maximum memory (in bytes) to buffer before
	// flushing to a temp file. ~768MB matches samtools sort default.
	defaultMaxMemory = 768 * 1024 * 1024
)

// sortedBamWriter buffers SamRecords in memory, sorts them, and writes
// sorted temp BAM files. On Close, temp files are merge-sorted into the
// final output. Supports coordinate sort and name sort.
type sortedBamWriter struct {
	filename  string
	header    *SamHeader
	refs      []bamRefInfo
	refIdx    map[string]int32
	sortCoord bool // true for coordinate sort, false for name sort

	tmpPrefix string   // prefix for temp file names
	tmpFiles  []string // paths of sorted temp BAM files
	tmpCount  int      // number of temp files created

	buf     []*SamRecord // in-memory buffer
	bufSize int          // approximate memory usage of buf in bytes
	maxMem  int          // max memory before flush

	mu     sync.Mutex
	closed bool
	err    error
}

func newSortedBamWriter(filename string, header *SamHeader, sortCoord bool, tmpPrefix string) (*sortedBamWriter, error) {
	sw := &sortedBamWriter{
		filename:  filename,
		header:    header,
		sortCoord: sortCoord,
		tmpPrefix: tmpPrefix,
		maxMem:    defaultMaxMemory,
	}

	if sw.tmpPrefix == "" {
		sw.tmpPrefix = filename + ".tmp"
	}

	// Build reference list from header.
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

// Write buffers a record. When the buffer exceeds maxMem, it is sorted
// and flushed to a temp BAM file.
func (sw *sortedBamWriter) Write(rec *SamRecord) error {
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

// Close sorts and flushes any remaining buffered records, then merge-sorts
// all temp files into the final output.
func (sw *sortedBamWriter) Close() error {
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

	// If everything fits in one buffer (no temp files), write directly.
	if len(sw.tmpFiles) == 0 {
		sw.sortBuffer()
		if err := sw.writeBAM(sw.filename, sw.buf); err != nil {
			return err
		}
		return nil
	}

	// Flush remaining buffer to a temp file.
	if len(sw.buf) > 0 {
		if err := sw.flushBuffer(); err != nil {
			sw.cleanup()
			return err
		}
	}

	// Merge-sort temp files into final output.
	if err := sw.mergeFiles(); err != nil {
		sw.cleanup()
		return err
	}

	sw.cleanup()
	return nil
}

// flushBuffer sorts the current buffer and writes it to a temp BAM file.
func (sw *sortedBamWriter) flushBuffer() error {
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

// sortBuffer sorts the buffer using the configured sort order.
func (sw *sortedBamWriter) sortBuffer() {
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

// lessCoord implements samtools coordinate sort order:
// unmapped (refID -1) at end, then by refID, then by pos.
func (sw *sortedBamWriter) lessCoord(a, b *SamRecord) bool {
	aRef := sw.refID(a)
	bRef := sw.refID(b)

	// Unmapped reads go at the end.
	aUnmapped := aRef == -1
	bUnmapped := bRef == -1
	if aUnmapped != bUnmapped {
		return !aUnmapped // mapped before unmapped
	}
	if aRef != bRef {
		return aRef < bRef
	}
	return a.Pos < b.Pos
}

// lessName implements samtools name sort order: plain lexicographic
// comparison of the read name.
func (sw *sortedBamWriter) lessName(a, b *SamRecord) bool {
	if a.ReadName != b.ReadName {
		return a.ReadName < b.ReadName
	}
	return a.Flag < b.Flag
}

func (sw *sortedBamWriter) refID(rec *SamRecord) int32 {
	if rec.RefName == "*" {
		return -1
	}
	if id, ok := sw.refIdx[rec.RefName]; ok {
		return id
	}
	return -1
}

// writeBAM writes a slice of records to a BAM file.
func (sw *sortedBamWriter) writeBAM(filename string, records []*SamRecord) error {
	w, err := newBamWriter(filename, sw.header)
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

// mergeFiles performs a k-way merge of sorted temp BAM files into the
// final output file.
func (sw *sortedBamWriter) mergeFiles() error {
	// Open all temp files as BAM readers.
	readers := make([]*BamReader, len(sw.tmpFiles))
	for i, path := range sw.tmpFiles {
		f, err := os.Open(path)
		if err != nil {
			// Close already-opened readers.
			for j := 0; j < i; j++ {
				readers[j].Close()
			}
			return fmt.Errorf("opening temp file %s: %w", path, err)
		}
		r, err := NewBamReader(f)
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

	// Open the final output writer.
	outWriter, err := newBamWriter(sw.filename, sw.header)
	if err != nil {
		return err
	}

	// Initialize the merge heap.
	h := &mergeHeap{less: sw.lessCoord}
	if !sw.sortCoord {
		h.less = sw.lessName
	}

	for i, r := range readers {
		rec, err := r.Next()
		if err == io.EOF {
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

		// Refill from the same reader.
		rec, err := readers[entry.readerIdx].Next()
		if err == io.EOF {
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

// cleanup removes all temp files.
func (sw *sortedBamWriter) cleanup() {
	for _, path := range sw.tmpFiles {
		os.Remove(path)
	}
	sw.tmpFiles = nil
}

// estimateRecordSize returns a rough estimate of memory used by a SamRecord.
func estimateRecordSize(rec *SamRecord) int {
	size := 200 // base struct overhead
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

// mergeEntry holds a record and the index of the reader it came from.
type mergeEntry struct {
	rec       *SamRecord
	readerIdx int
}

// mergeHeap implements heap.Interface for k-way merge.
type mergeHeap struct {
	entries []*mergeEntry
	less    func(a, b *SamRecord) bool
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
