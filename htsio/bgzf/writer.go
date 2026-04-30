package bgzf

import (
	"compress/flate"
	"encoding/binary"
	"hash/crc32"
	"io"
	"os"
)

// Writer writes BGZF-compressed data. It implements io.WriteCloser.
//
// Data is buffered until the uncompressed block reaches MaxUncompressedSize,
// at which point it is flushed as a complete BGZF block. Call Close to flush
// the final partial block and write the EOF marker.
type Writer struct {
	w              io.Writer
	f              *os.File // non-nil if opened by NewBGZipFile
	buf            []byte   // uncompressed data pending flush
	level          int      // flate compression level
	closed         bool     // prevents double-close
	err            error    // sticky error
	compressedOff  int64    // compressed bytes written so far (block boundaries)
}

// NewWriter creates a BGZF writer using the default compression level.
func NewWriter(w io.Writer) *Writer {
	return &Writer{
		w:     w,
		buf:   make([]byte, 0, MaxUncompressedSize),
		level: flate.DefaultCompression,
	}
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

// NewWriterLevel creates a BGZF writer with the specified flate compression level.
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
		return w.err
	}
	if len(w.buf) > 0 {
		if err := w.flush(); err != nil {
			w.err = err
			return err
		}
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

// maxCompressedPayload is the space available for DEFLATE output within a
// single BGZF block after accounting for the header and trailer.
const maxCompressedPayload = MaxBlockSize - bgzfHeaderSize - gzipTrailerSize

// flush compresses w.buf into one or more BGZF blocks and writes them.
// If the compressed output exceeds the block size limit (possible with
// incompressible data), the input is split in half and each half is
// flushed recursively.
func (w *Writer) flush() error {
	data := w.buf
	if len(data) == 0 {
		return nil
	}

	compressed, err := w.deflate(data)
	if err != nil {
		return err
	}

	if len(compressed) > maxCompressedPayload {
		// Compressed output doesn't fit in one block. Split the input
		// in half and flush each part separately.
		mid := len(data) / 2
		w.buf = data[:mid]
		if err := w.flush(); err != nil {
			return err
		}
		w.buf = data[mid:]
		if err := w.flush(); err != nil {
			return err
		}
		w.buf = w.buf[:0]
		return nil
	}

	crc := crc32.ChecksumIEEE(data)
	isize := uint32(len(data))
	blockSize := bgzfHeaderSize + len(compressed) + gzipTrailerSize
	bsize := uint16(blockSize - 1)

	// Write the header.
	var hdr [bgzfHeaderSize]byte
	hdr[0] = 0x1f // ID1
	hdr[1] = 0x8b // ID2
	hdr[2] = 0x08 // CM = deflate
	hdr[3] = 0x04 // FLG = FEXTRA
	// MTIME, XFL, OS all zero
	hdr[9] = 0xff // OS = unknown
	binary.LittleEndian.PutUint16(hdr[10:12], 6)    // XLEN
	hdr[12] = 'B'                                    // SI1
	hdr[13] = 'C'                                    // SI2
	binary.LittleEndian.PutUint16(hdr[14:16], 2)     // SLEN
	binary.LittleEndian.PutUint16(hdr[16:18], bsize) // BSIZE

	if _, err := w.w.Write(hdr[:]); err != nil {
		return err
	}
	if _, err := w.w.Write(compressed); err != nil {
		return err
	}

	// Write trailer: CRC32 + ISIZE.
	var trailer [gzipTrailerSize]byte
	binary.LittleEndian.PutUint32(trailer[0:4], crc)
	binary.LittleEndian.PutUint32(trailer[4:8], isize)
	if _, err := w.w.Write(trailer[:]); err != nil {
		return err
	}

	w.compressedOff += int64(blockSize)
	w.buf = w.buf[:0]
	return nil
}

// deflate compresses data using the writer's compression level and returns
// the raw DEFLATE output.
func (w *Writer) deflate(data []byte) ([]byte, error) {
	cbw := &sliceWriter{buf: make([]byte, 0, len(data)+256)}
	fw, err := flate.NewWriter(cbw, w.level)
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
