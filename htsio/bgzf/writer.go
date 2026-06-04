package bgzf

import (
	"compress/flate"
	"encoding/binary"
	"hash/crc32"
	"io"
	"os"
	"runtime"
)

// Writer writes BGZF-compressed data. It implements io.WriteCloser.
//
// Data is buffered until the uncompressed block reaches MaxUncompressedSize,
// at which point it is flushed as a complete BGZF block. Call Close to flush
// the final partial block and write the EOF marker.
//
// When created with NewParallelWriter, compression is performed by a pool of
// goroutines and written in order by a background drain goroutine. The
// Write/Flush/Close API is unchanged.
type Writer struct {
	w     io.Writer
	f     *os.File // non-nil if opened by NewBGZipFile
	buf   []byte   // uncompressed data pending flush
	level int      // flate compression level

	closed        bool  // prevents double-close
	err           error // sticky error
	compressedOff int64 // compressed bytes written so far (block boundaries)

	// Parallel compression (zero values → single-threaded mode).
	nWorkers int
	workCh   chan compressJob       // uncompressed blocks → workers
	drainCh  chan chan compressResult // result channels → drain goroutine (preserves order)
	drainErr chan error              // drain goroutine final error
}

// compressJob is a unit of work sent to a compression worker.
type compressJob struct {
	data   []byte
	level  int
	result chan<- compressResult
}

// compressResult is the output of compressing one input buffer.
// blocks contains one or more complete BGZF blocks concatenated together
// (split may occur for incompressible data).
type compressResult struct {
	blocks []byte
	err    error
}

// NewWriter creates a BGZF writer using the default compression level
// (single-threaded).
func NewWriter(w io.Writer) *Writer {
	return &Writer{
		w:     w,
		buf:   make([]byte, 0, MaxUncompressedSize),
		level: flate.DefaultCompression,
	}
}

// NewParallelWriter creates a BGZF writer that compresses blocks in parallel
// using the given number of goroutines. If threads <= 1, behaves identically
// to NewWriter. A threads value of 0 uses runtime.NumCPU().
func NewParallelWriter(w io.Writer, threads int) *Writer {
	return newParallelWriter(w, flate.DefaultCompression, threads)
}

// NewParallelWriterLevel creates a parallel BGZF writer with a custom
// compression level.
func NewParallelWriterLevel(w io.Writer, level int, threads int) (*Writer, error) {
	if level < flate.HuffmanOnly || level > flate.BestCompression {
		if level != flate.DefaultCompression {
			return nil, flate.CorruptInputError(int64(level))
		}
	}
	return newParallelWriter(w, level, threads), nil
}

func newParallelWriter(w io.Writer, level int, threads int) *Writer {
	if threads == 0 {
		threads = runtime.NumCPU()
	}
	bw := &Writer{
		w:     w,
		buf:   make([]byte, 0, MaxUncompressedSize),
		level: level,
	}
	if threads > 1 {
		bw.nWorkers = threads
		bw.workCh = make(chan compressJob, threads)
		bw.drainCh = make(chan chan compressResult, threads*2)
		bw.drainErr = make(chan error, 1)

		// Start compression workers.
		for range threads {
			go compressWorker(bw.workCh)
		}

		// Start drain goroutine — writes compressed blocks to output in order.
		go bw.drainLoop()
	}
	return bw
}

// compressWorker reads jobs from workCh and sends results back on each job's
// result channel.
func compressWorker(workCh <-chan compressJob) {
	for job := range workCh {
		blocks, err := compressBlocks(job.data, job.level)
		job.result <- compressResult{blocks: blocks, err: err}
	}
}

// drainLoop reads result channels from drainCh in order, waits for each
// result, and writes the compressed blocks to the output writer.
func (w *Writer) drainLoop() {
	var err error
	for resultCh := range w.drainCh {
		if err != nil {
			// Already errored — drain remaining channels to unblock workers.
			<-resultCh
			continue
		}
		res := <-resultCh
		if res.err != nil {
			err = res.err
			continue
		}
		if _, werr := w.w.Write(res.blocks); werr != nil {
			err = werr
			continue
		}
		w.compressedOff += int64(len(res.blocks))
	}
	w.drainErr <- err
}

// NewBGZipFile creates a Writer that writes to the named file. The file is
// created (or truncated) and will be closed when the writer is closed.
func NewBGZipFile(filename string) (*Writer, error) {
	f, err := os.Create(filename)
	if err != nil {
		return nil, err
	}
	w := NewWriter(f)
	w.f = f
	return w, nil
}

// NewWriterLevel creates a BGZF writer with the specified flate compression
// level (single-threaded).
func NewWriterLevel(w io.Writer, level int) (*Writer, error) {
	if level < flate.HuffmanOnly || level > flate.BestCompression {
		if level != flate.DefaultCompression {
			return nil, flate.CorruptInputError(int64(level))
		}
	}
	return &Writer{
		w:     w,
		buf:   make([]byte, 0, MaxUncompressedSize),
		level: level,
	}, nil
}

// Write implements io.Writer. Data is buffered and flushed as complete BGZF
// blocks when the buffer reaches MaxUncompressedSize.
func (w *Writer) Write(p []byte) (int, error) {
	if w.err != nil {
		return 0, w.err
	}

	total := 0
	for len(p) > 0 {
		space := MaxUncompressedSize - len(w.buf)
		if space == 0 {
			if err := w.flush(); err != nil {
				w.err = err
				return total, err
			}
			space = MaxUncompressedSize
		}

		n := len(p)
		if n > space {
			n = space
		}
		w.buf = append(w.buf, p[:n]...)
		p = p[n:]
		total += n
	}
	return total, nil
}

// VirtualTell returns the virtual offset where the next byte will be written.
// This is the compressed block offset of the current (unflushed) block
// combined with the uncompressed position within it.
//
// Note: in parallel mode, compressedOff reflects only blocks that have been
// drained to the output so far. Call Flush to ensure all pending blocks are
// written before relying on VirtualTell.
func (w *Writer) VirtualTell() VirtualOffset {
	return NewVirtualOffset(w.compressedOff, uint16(len(w.buf)))
}

// Flush writes the current buffer as a BGZF block, even if it is not full.
func (w *Writer) Flush() error {
	if w.err != nil {
		return w.err
	}
	if len(w.buf) == 0 {
		return nil
	}
	return w.flush()
}

// Close flushes any remaining data, writes the BGZF EOF block, and closes
// the underlying file if the writer was created with NewBGZipFile. Close is
// idempotent.
func (w *Writer) Close() error {
	if w.closed {
		return nil
	}
	w.closed = true
	if w.err != nil {
		w.shutdownWorkers()
		return w.err
	}
	if len(w.buf) > 0 {
		if err := w.flush(); err != nil {
			w.err = err
			w.shutdownWorkers()
			return err
		}
	}

	// Shut down parallel pipeline and collect drain error.
	if err := w.shutdownWorkers(); err != nil {
		w.err = err
		return err
	}

	// Write the standard EOF block.
	if _, err := w.w.Write(bgzfEOFBlock); err != nil {
		w.err = err
		return err
	}
	if w.f != nil {
		return w.f.Close()
	}
	return nil
}

// shutdownWorkers closes the worker pool and drain goroutine (if parallel),
// returning any error from the drain goroutine.
func (w *Writer) shutdownWorkers() error {
	if w.workCh == nil {
		return nil
	}
	close(w.workCh)
	close(w.drainCh)
	err := <-w.drainErr
	w.workCh = nil
	w.drainCh = nil
	return err
}

// maxCompressedPayload is the space available for DEFLATE output within a
// single BGZF block after accounting for the header and trailer.
const maxCompressedPayload = MaxBlockSize - bgzfHeaderSize - gzipTrailerSize

// flush compresses w.buf into one or more BGZF blocks and writes them.
func (w *Writer) flush() error {
	data := w.buf
	if len(data) == 0 {
		return nil
	}

	if w.nWorkers > 1 {
		return w.flushAsync(data)
	}
	return w.flushSync(data)
}

// flushSync compresses and writes a block synchronously (single-threaded path).
func (w *Writer) flushSync(data []byte) error {
	blocks, err := compressBlocks(data, w.level)
	if err != nil {
		return err
	}
	if _, err := w.w.Write(blocks); err != nil {
		return err
	}
	w.compressedOff += int64(len(blocks))
	w.buf = w.buf[:0]
	return nil
}

// flushAsync sends the current buffer to the worker pool for compression.
func (w *Writer) flushAsync(data []byte) error {
	// Copy data — the buffer will be reused.
	cp := make([]byte, len(data))
	copy(cp, data)
	w.buf = w.buf[:0]

	resultCh := make(chan compressResult, 1)
	w.workCh <- compressJob{
		data:   cp,
		level:  w.level,
		result: resultCh,
	}
	w.drainCh <- resultCh
	return nil
}

// compressBlocks compresses data into one or more complete BGZF blocks.
// If the compressed output exceeds the block size limit (incompressible data),
// the input is split in half and each half is compressed recursively.
// Returns the concatenated block bytes.
func compressBlocks(data []byte, level int) ([]byte, error) {
	compressed, err := deflateData(data, level)
	if err != nil {
		return nil, err
	}

	if len(compressed) > maxCompressedPayload {
		mid := len(data) / 2
		first, err := compressBlocks(data[:mid], level)
		if err != nil {
			return nil, err
		}
		second, err := compressBlocks(data[mid:], level)
		if err != nil {
			return nil, err
		}
		return append(first, second...), nil
	}

	return buildBGZFBlock(data, compressed), nil
}

// buildBGZFBlock assembles a complete BGZF block from uncompressed and
// compressed data.
func buildBGZFBlock(uncompressed, compressed []byte) []byte {
	crc := crc32.ChecksumIEEE(uncompressed)
	isize := uint32(len(uncompressed))
	blockSize := bgzfHeaderSize + len(compressed) + gzipTrailerSize
	bsize := uint16(blockSize - 1)

	block := make([]byte, blockSize)

	// Header.
	block[0] = 0x1f // ID1
	block[1] = 0x8b // ID2
	block[2] = 0x08 // CM = deflate
	block[3] = 0x04 // FLG = FEXTRA
	// MTIME, XFL, OS all zero
	block[9] = 0xff // OS = unknown
	binary.LittleEndian.PutUint16(block[10:12], 6)    // XLEN
	block[12] = 'B'                                    // SI1
	block[13] = 'C'                                    // SI2
	binary.LittleEndian.PutUint16(block[14:16], 2)     // SLEN
	binary.LittleEndian.PutUint16(block[16:18], bsize) // BSIZE

	// Compressed data.
	copy(block[bgzfHeaderSize:], compressed)

	// Trailer: CRC32 + ISIZE.
	trailerOff := bgzfHeaderSize + len(compressed)
	binary.LittleEndian.PutUint32(block[trailerOff:], crc)
	binary.LittleEndian.PutUint32(block[trailerOff+4:], isize)

	return block
}

// deflateData compresses data using the given flate compression level.
func deflateData(data []byte, level int) ([]byte, error) {
	cbw := &sliceWriter{buf: make([]byte, 0, len(data)+256)}
	fw, err := flate.NewWriter(cbw, level)
	if err != nil {
		return nil, err
	}
	if _, err := fw.Write(data); err != nil {
		return nil, err
	}
	if err := fw.Close(); err != nil {
		return nil, err
	}
	return cbw.buf, nil
}

// sliceWriter is a simple byte slice writer used to capture compressed output.
type sliceWriter struct {
	buf []byte
}

func (w *sliceWriter) Write(p []byte) (int, error) {
	w.buf = append(w.buf, p...)
	return len(p), nil
}
