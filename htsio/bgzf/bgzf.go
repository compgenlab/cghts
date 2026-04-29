// Package bgzf implements reading and writing of BGZF (Blocked GNU Zip Format)
// files as defined in the SAM/BAM specification. BGZF is a multi-member gzip
// format where each member (block) contains at most 64 KiB of uncompressed data
// and carries an extra field (BSIZE) recording the total block size minus one.
//
// Virtual offsets encode both the compressed block offset and the uncompressed
// offset within the block as (blockOffset << 16 | withinBlockOffset).
package bgzf

import (
	"bufio"
	"compress/flate"
	"encoding/binary"
	"fmt"
	"hash/crc32"
	"io"
)

const (
	// MaxBlockSize is the maximum compressed block size (including header/trailer).
	MaxBlockSize = 65536

	// MaxUncompressedSize is the maximum uncompressed payload per block.
	MaxUncompressedSize = 65536

	// bgzfHeaderSize is the fixed size of a BGZF block header (standard gzip
	// header fields + the BGZF extra subfield).
	bgzfHeaderSize = 18

	// gzipTrailerSize is the size of the gzip trailer (CRC32 + ISIZE).
	gzipTrailerSize = 8
)

// bgzfEOFBlock is the 28-byte empty BGZF block that marks end-of-file.
var bgzfEOFBlock = []byte{
	0x1f, 0x8b, 0x08, 0x04, // ID1, ID2, CM, FLG (FEXTRA)
	0x00, 0x00, 0x00, 0x00, // MTIME
	0x00, 0xff, // XFL, OS
	0x06, 0x00, // XLEN = 6
	0x42, 0x43, // SI1, SI2 ('BC')
	0x02, 0x00, // SLEN = 2
	0x1b, 0x00, // BSIZE = 27 (block size - 1)
	0x03, 0x00, // compressed empty DEFLATE block
	0x00, 0x00, 0x00, 0x00, // CRC32
	0x00, 0x00, 0x00, 0x00, // ISIZE = 0
}

// VirtualOffset encodes a position in a BGZF file as the compressed block
// offset in the upper 48 bits and the uncompressed offset within the block
// in the lower 16 bits.
type VirtualOffset uint64

// NewVirtualOffset creates a VirtualOffset from a block offset and an
// uncompressed offset within the block.
func NewVirtualOffset(blockOffset int64, withinBlock uint16) VirtualOffset {
	return VirtualOffset(uint64(blockOffset)<<16 | uint64(withinBlock))
}

// BlockOffset returns the compressed file offset of the BGZF block.
func (v VirtualOffset) BlockOffset() int64 {
	return int64(uint64(v) >> 16)
}

// WithinBlock returns the uncompressed offset within the BGZF block.
func (v VirtualOffset) WithinBlock() uint16 {
	return uint16(v)
}

// blockHeader holds the parsed fields from a BGZF block header.
type blockHeader struct {
	bsize uint16 // total block size minus 1
}

// readBlockHeader reads and validates a BGZF block header from r.
// Returns io.EOF if there are no more bytes to read.
func readBlockHeader(r io.Reader) (*blockHeader, error) {
	var buf [bgzfHeaderSize]byte
	_, err := io.ReadFull(r, buf[:])
	if err != nil {
		if err == io.ErrUnexpectedEOF {
			return nil, fmt.Errorf("bgzf: truncated block header")
		}
		return nil, err // io.EOF if no bytes at all
	}

	// Validate gzip magic
	if buf[0] != 0x1f || buf[1] != 0x8b {
		return nil, fmt.Errorf("bgzf: invalid gzip magic: %#x %#x", buf[0], buf[1])
	}
	// CM must be 8 (deflate)
	if buf[2] != 8 {
		return nil, fmt.Errorf("bgzf: unsupported compression method: %d", buf[2])
	}
	// FLG must have FEXTRA set (bit 2)
	if buf[3]&0x04 == 0 {
		return nil, fmt.Errorf("bgzf: FEXTRA flag not set")
	}
	// XLEN = 6
	xlen := binary.LittleEndian.Uint16(buf[10:12])
	if xlen != 6 {
		return nil, fmt.Errorf("bgzf: unexpected XLEN: %d (expected 6)", xlen)
	}
	// BC subfield
	if buf[12] != 'B' || buf[13] != 'C' {
		return nil, fmt.Errorf("bgzf: missing BC extra subfield: %c%c", buf[12], buf[13])
	}
	// SLEN = 2
	slen := binary.LittleEndian.Uint16(buf[14:16])
	if slen != 2 {
		return nil, fmt.Errorf("bgzf: unexpected SLEN: %d (expected 2)", slen)
	}

	bsize := binary.LittleEndian.Uint16(buf[16:18])
	return &blockHeader{bsize: bsize}, nil
}

// Reader reads BGZF-compressed data. It implements io.Reader and io.ByteReader.
//
// Data is decompressed one block at a time. The current virtual offset can be
// queried with VirtualOffset() at any point during reading.
type Reader struct {
	r   *bufio.Reader
	buf []byte // decompressed data for the current block
	pos int    // read position within buf

	blockOffset int64 // compressed file offset of the current block
	fileOffset  int64 // compressed file offset (total bytes consumed from r)

	err error // sticky error
}

// NewReader creates a BGZF reader that decompresses data from r.
func NewReader(r io.Reader) *Reader {
	return &Reader{
		r: bufio.NewReaderSize(r, MaxBlockSize*2),
	}
}

// VirtualTell returns the virtual offset of the next byte that will be read.
func (r *Reader) VirtualTell() VirtualOffset {
	return NewVirtualOffset(r.blockOffset, uint16(r.pos))
}

// Read implements io.Reader.
func (r *Reader) Read(p []byte) (int, error) {
	if r.err != nil {
		return 0, r.err
	}

	total := 0
	for len(p) > 0 {
		// If current block is exhausted, load the next one.
		if r.pos >= len(r.buf) {
			if err := r.nextBlock(); err != nil {
				r.err = err
				if total > 0 {
					return total, nil
				}
				return 0, err
			}
		}

		n := copy(p, r.buf[r.pos:])
		r.pos += n
		p = p[n:]
		total += n
	}
	return total, nil
}

// ReadByte implements io.ByteReader.
func (r *Reader) ReadByte() (byte, error) {
	if r.err != nil {
		return 0, r.err
	}
	if r.pos >= len(r.buf) {
		if err := r.nextBlock(); err != nil {
			r.err = err
			return 0, err
		}
	}
	b := r.buf[r.pos]
	r.pos++
	return b, nil
}

// nextBlock reads and decompresses the next BGZF block.
func (r *Reader) nextBlock() error {
	r.blockOffset = r.fileOffset

	hdr, err := readBlockHeader(r.r)
	if err != nil {
		return err
	}

	// Total block size = bsize + 1. The header is already consumed.
	blockSize := int(hdr.bsize) + 1
	remaining := blockSize - bgzfHeaderSize
	if remaining < gzipTrailerSize {
		return fmt.Errorf("bgzf: block too small: bsize=%d", hdr.bsize)
	}

	// Read the rest of the block (compressed data + trailer).
	blockData := make([]byte, remaining)
	if _, err := io.ReadFull(r.r, blockData); err != nil {
		if err == io.ErrUnexpectedEOF {
			return fmt.Errorf("bgzf: truncated block data")
		}
		return fmt.Errorf("bgzf: reading block: %w", err)
	}

	r.fileOffset += int64(blockSize)

	compressedData := blockData[:len(blockData)-gzipTrailerSize]
	trailer := blockData[len(blockData)-gzipTrailerSize:]

	expectedCRC := binary.LittleEndian.Uint32(trailer[0:4])
	expectedSize := binary.LittleEndian.Uint32(trailer[4:8])

	// Empty block (EOF marker) — size 0 is valid.
	if expectedSize == 0 {
		r.buf = r.buf[:0]
		r.pos = 0
		// Peek to see if there's more data. An EOF block at the end
		// of the file should return io.EOF.
		if _, err := r.r.Peek(1); err != nil {
			return io.EOF
		}
		return nil
	}

	if expectedSize > MaxUncompressedSize {
		return fmt.Errorf("bgzf: uncompressed size %d exceeds maximum %d", expectedSize, MaxUncompressedSize)
	}

	// Decompress using raw DEFLATE (no gzip wrapper).
	fr := flate.NewReader(nil)
	resetter, ok := fr.(flate.Resetter)
	if !ok {
		return fmt.Errorf("bgzf: flate.Reader does not implement Resetter")
	}
	if err := resetter.Reset(io.NopCloser(newBytesReader(compressedData)), nil); err != nil {
		return fmt.Errorf("bgzf: flate reset: %w", err)
	}

	// Reuse r.buf if it has enough capacity.
	if cap(r.buf) >= int(expectedSize) {
		r.buf = r.buf[:expectedSize]
	} else {
		r.buf = make([]byte, expectedSize)
	}

	n, err := io.ReadFull(fr, r.buf)
	if err != nil {
		return fmt.Errorf("bgzf: decompressing block: %w", err)
	}
	r.buf = r.buf[:n]
	fr.Close()

	// Verify CRC32.
	actualCRC := crc32.ChecksumIEEE(r.buf)
	if actualCRC != expectedCRC {
		return fmt.Errorf("bgzf: CRC32 mismatch: expected %#x, got %#x", expectedCRC, actualCRC)
	}

	r.pos = 0
	return nil
}

// bytesReader is a minimal io.Reader over a byte slice, avoiding
// the overhead of bytes.NewReader which also implements io.Seeker.
type bytesReader struct {
	data []byte
	pos  int
}

func newBytesReader(data []byte) *bytesReader {
	return &bytesReader{data: data}
}

func (b *bytesReader) Read(p []byte) (int, error) {
	if b.pos >= len(b.data) {
		return 0, io.EOF
	}
	n := copy(p, b.data[b.pos:])
	b.pos += n
	return n, nil
}
