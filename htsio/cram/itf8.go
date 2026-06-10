package cram

import (
	"fmt"
	"io"
)

// readITF8 reads an ITF8-encoded int32 from a byte stream.
// ITF8 uses a variable number of bytes (1-5) determined by prefix bits:
//
//	0xxxxxxx                             → 1 byte,  7 bits
//	10xxxxxx xxxxxxxx                    → 2 bytes, 14 bits
//	110xxxxx xxxxxxxx xxxxxxxx           → 3 bytes, 21 bits
//	1110xxxx xxxxxxxx xxxxxxxx xxxxxxxx  → 4 bytes, 28 bits
//	11110000 xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx → 5 bytes, 32 bits
func readITF8(r io.Reader) (int32, error) {
	var b [1]byte
	if _, err := io.ReadFull(r, b[:]); err != nil {
		return 0, err
	}
	first := b[0]

	if first&0x80 == 0 {
		// 0xxxxxxx
		return int32(first), nil
	}
	if first&0x40 == 0 {
		// 10xxxxxx
		if _, err := io.ReadFull(r, b[:]); err != nil {
			return 0, err
		}
		return int32(first&0x3f)<<8 | int32(b[0]), nil
	}
	if first&0x20 == 0 {
		// 110xxxxx
		var buf [2]byte
		if _, err := io.ReadFull(r, buf[:]); err != nil {
			return 0, err
		}
		return int32(first&0x1f)<<16 | int32(buf[0])<<8 | int32(buf[1]), nil
	}
	if first&0x10 == 0 {
		// 1110xxxx
		var buf [3]byte
		if _, err := io.ReadFull(r, buf[:]); err != nil {
			return 0, err
		}
		return int32(first&0x0f)<<24 | int32(buf[0])<<16 | int32(buf[1])<<8 | int32(buf[2]), nil
	}

	// 1111xxxx — 5 bytes: 4 bits from first byte + 8+8+8+4 from next 4 bytes = 32 bits
	var buf [4]byte
	if _, err := io.ReadFull(r, buf[:]); err != nil {
		return 0, err
	}
	return int32(first&0x0f)<<28 | int32(buf[0])<<20 | int32(buf[1])<<12 | int32(buf[2])<<4 | int32(buf[3]&0x0f), nil
}

// readLTF8 reads an LTF8-encoded int64 from a byte stream.
// LTF8 extends ITF8 to 64 bits using up to 9 bytes.
func readLTF8(r io.Reader) (int64, error) {
	var b [1]byte
	if _, err := io.ReadFull(r, b[:]); err != nil {
		return 0, err
	}
	first := b[0]

	if first&0x80 == 0 {
		return int64(first), nil
	}
	if first&0x40 == 0 {
		if _, err := io.ReadFull(r, b[:]); err != nil {
			return 0, err
		}
		return int64(first&0x3f)<<8 | int64(b[0]), nil
	}
	if first&0x20 == 0 {
		var buf [2]byte
		if _, err := io.ReadFull(r, buf[:]); err != nil {
			return 0, err
		}
		return int64(first&0x1f)<<16 | int64(buf[0])<<8 | int64(buf[1]), nil
	}
	if first&0x10 == 0 {
		var buf [3]byte
		if _, err := io.ReadFull(r, buf[:]); err != nil {
			return 0, err
		}
		return int64(first&0x0f)<<24 | int64(buf[0])<<16 | int64(buf[1])<<8 | int64(buf[2]), nil
	}
	if first&0x08 == 0 {
		var buf [4]byte
		if _, err := io.ReadFull(r, buf[:]); err != nil {
			return 0, err
		}
		return int64(first&0x07)<<32 | int64(buf[0])<<24 | int64(buf[1])<<16 | int64(buf[2])<<8 | int64(buf[3]), nil
	}
	if first&0x04 == 0 {
		var buf [5]byte
		if _, err := io.ReadFull(r, buf[:]); err != nil {
			return 0, err
		}
		return int64(first&0x03)<<40 | int64(buf[0])<<32 | int64(buf[1])<<24 | int64(buf[2])<<16 | int64(buf[3])<<8 | int64(buf[4]), nil
	}
	if first&0x02 == 0 {
		var buf [6]byte
		if _, err := io.ReadFull(r, buf[:]); err != nil {
			return 0, err
		}
		return int64(first&0x01)<<48 | int64(buf[0])<<40 | int64(buf[1])<<32 | int64(buf[2])<<24 | int64(buf[3])<<16 | int64(buf[4])<<8 | int64(buf[5]), nil
	}
	if first&0x01 == 0 {
		var buf [7]byte
		if _, err := io.ReadFull(r, buf[:]); err != nil {
			return 0, err
		}
		return int64(buf[0])<<48 | int64(buf[1])<<40 | int64(buf[2])<<32 | int64(buf[3])<<24 | int64(buf[4])<<16 | int64(buf[5])<<8 | int64(buf[6]), nil
	}

	// 0xFF prefix — full 64-bit value
	var buf [8]byte
	if _, err := io.ReadFull(r, buf[:]); err != nil {
		return 0, err
	}
	return int64(buf[0])<<56 | int64(buf[1])<<48 | int64(buf[2])<<40 | int64(buf[3])<<32 |
		int64(buf[4])<<24 | int64(buf[5])<<16 | int64(buf[6])<<8 | int64(buf[7]), nil
}

// readSizedBytes reads exactly size bytes from r, where size was read from
// untrusted input. It rejects a negative size and reads through a LimitReader
// so a bogus (huge) size cannot force a multi-gigabyte allocation up front: the
// buffer grows with the bytes actually present, and a short stream is reported
// as a truncation error rather than returning a partially-zeroed buffer.
func readSizedBytes(r io.Reader, size int32) ([]byte, error) {
	if size < 0 {
		return nil, fmt.Errorf("negative block size: %d", size)
	}
	data, err := io.ReadAll(io.LimitReader(r, int64(size)))
	if err != nil {
		return nil, err
	}
	if len(data) != int(size) {
		return nil, fmt.Errorf("truncated block: got %d of %d bytes", len(data), size)
	}
	return data, nil
}

// readITF8Array reads an ITF8 length followed by that many ITF8 values.
func readITF8Array(r io.Reader) ([]int32, error) {
	n, err := readITF8(r)
	if err != nil {
		return nil, fmt.Errorf("reading array length: %w", err)
	}
	if n < 0 {
		return nil, fmt.Errorf("negative array length: %d", n)
	}
	// Grow via append rather than make([]int32, n): n is attacker-controlled,
	// so a huge value would otherwise allocate gigabytes up front. Appending
	// caps memory at the number of elements actually present in the stream — a
	// lying length simply fails the read below at EOF.
	vals := make([]int32, 0, min(int(n), 1024))
	for i := int32(0); i < n; i++ {
		v, err := readITF8(r)
		if err != nil {
			return nil, fmt.Errorf("reading array element %d: %w", i, err)
		}
		vals = append(vals, v)
	}
	return vals, nil
}

// bitReader reads individual bits from a byte slice, MSB first.
type bitReader struct {
	data []byte
	pos  int // bit position
}

func newBitReader(data []byte) *bitReader {
	return &bitReader{data: data}
}

// readBits reads n bits (up to 32) and returns them as an int32.
func (br *bitReader) readBits(n int) (int32, error) {
	if n <= 0 || n > 32 {
		return 0, fmt.Errorf("invalid bit count: %d", n)
	}

	var val int32
	for i := 0; i < n; i++ {
		byteIdx := br.pos / 8
		bitIdx := 7 - (br.pos % 8) // MSB first
		if byteIdx >= len(br.data) {
			return 0, io.ErrUnexpectedEOF
		}
		bit := (br.data[byteIdx] >> uint(bitIdx)) & 1
		val = (val << 1) | int32(bit)
		br.pos++
	}
	return val, nil
}

// readBit reads a single bit.
func (br *bitReader) readBit() (int, error) {
	byteIdx := br.pos / 8
	bitIdx := 7 - (br.pos % 8)
	if byteIdx >= len(br.data) {
		return 0, io.ErrUnexpectedEOF
	}
	bit := int((br.data[byteIdx] >> uint(bitIdx)) & 1)
	br.pos++
	return bit, nil
}

// readByte reads 8 bits as a byte.
func (br *bitReader) readByte() (byte, error) {
	v, err := br.readBits(8)
	if err != nil {
		return 0, err
	}
	return byte(v), nil
}

// writeITF8 writes an ITF8-encoded int32 to a byte stream.
func writeITF8(w io.Writer, val int32) error {
	v := uint32(val)
	switch {
	case v < 0x80:
		_, err := w.Write([]byte{byte(v)})
		return err
	case v < 0x4000:
		_, err := w.Write([]byte{byte(v>>8) | 0x80, byte(v)})
		return err
	case v < 0x200000:
		_, err := w.Write([]byte{byte(v>>16) | 0xC0, byte(v >> 8), byte(v)})
		return err
	case v < 0x10000000:
		_, err := w.Write([]byte{byte(v>>24) | 0xE0, byte(v >> 16), byte(v >> 8), byte(v)})
		return err
	default:
		_, err := w.Write([]byte{byte(v>>28) | 0xF0, byte(v >> 20), byte(v >> 12), byte(v >> 4), byte(v & 0x0f)})
		return err
	}
}

// writeLTF8 writes an LTF8-encoded int64 to a byte stream.
func writeLTF8(w io.Writer, val int64) error {
	v := uint64(val)
	switch {
	case v < 0x80:
		_, err := w.Write([]byte{byte(v)})
		return err
	case v < 0x4000:
		_, err := w.Write([]byte{byte(v>>8) | 0x80, byte(v)})
		return err
	case v < 0x200000:
		_, err := w.Write([]byte{byte(v>>16) | 0xC0, byte(v >> 8), byte(v)})
		return err
	case v < 0x10000000:
		_, err := w.Write([]byte{byte(v>>24) | 0xE0, byte(v >> 16), byte(v >> 8), byte(v)})
		return err
	case v < 0x800000000:
		_, err := w.Write([]byte{byte(v>>32) | 0xF0, byte(v >> 24), byte(v >> 16), byte(v >> 8), byte(v)})
		return err
	case v < 0x40000000000:
		_, err := w.Write([]byte{byte(v>>40) | 0xF8, byte(v >> 32), byte(v >> 24), byte(v >> 16), byte(v >> 8), byte(v)})
		return err
	case v < 0x2000000000000:
		_, err := w.Write([]byte{byte(v>>48) | 0xFC, byte(v >> 40), byte(v >> 32), byte(v >> 24), byte(v >> 16), byte(v >> 8), byte(v)})
		return err
	case v < 0x100000000000000:
		_, err := w.Write([]byte{0xFE, byte(v >> 48), byte(v >> 40), byte(v >> 32), byte(v >> 24), byte(v >> 16), byte(v >> 8), byte(v)})
		return err
	default:
		_, err := w.Write([]byte{0xFF, byte(v >> 56), byte(v >> 48), byte(v >> 40), byte(v >> 32), byte(v >> 24), byte(v >> 16), byte(v >> 8), byte(v)})
		return err
	}
}

// writeITF8Array writes an ITF8 length followed by that many ITF8 values.
func writeITF8Array(w io.Writer, vals []int32) error {
	if err := writeITF8(w, int32(len(vals))); err != nil {
		return err
	}
	for _, v := range vals {
		if err := writeITF8(w, v); err != nil {
			return err
		}
	}
	return nil
}

// itf8Size returns the number of bytes needed to encode val as ITF8.
func itf8Size(val int32) int {
	v := uint32(val)
	switch {
	case v < 0x80:
		return 1
	case v < 0x4000:
		return 2
	case v < 0x200000:
		return 3
	case v < 0x10000000:
		return 4
	default:
		return 5
	}
}

// bitWriter writes individual bits to a byte slice, MSB first.
type bitWriter struct {
	data []byte
	pos  int // bit position
}

func newBitWriter() *bitWriter {
	return &bitWriter{data: make([]byte, 0, 1024)}
}

// writeBits writes n bits (up to 32) from val, MSB first.
func (bw *bitWriter) writeBits(val int32, n int) {
	for i := n - 1; i >= 0; i-- {
		byteIdx := bw.pos / 8
		bitIdx := 7 - (bw.pos % 8)
		for byteIdx >= len(bw.data) {
			bw.data = append(bw.data, 0)
		}
		if val&(1<<uint(i)) != 0 {
			bw.data[byteIdx] |= 1 << uint(bitIdx)
		}
		bw.pos++
	}
}

// writeBit writes a single bit.
func (bw *bitWriter) writeBit(bit int) {
	byteIdx := bw.pos / 8
	bitIdx := 7 - (bw.pos % 8)
	for byteIdx >= len(bw.data) {
		bw.data = append(bw.data, 0)
	}
	if bit != 0 {
		bw.data[byteIdx] |= 1 << uint(bitIdx)
	}
	bw.pos++
}

// bytes returns the accumulated bytes (padded with zero bits if needed).
func (bw *bitWriter) bytes() []byte {
	// Ensure the last partial byte is included.
	needed := (bw.pos + 7) / 8
	for len(bw.data) < needed {
		bw.data = append(bw.data, 0)
	}
	return bw.data[:needed]
}

// readITF8 reads an ITF8 value from the bit stream.
func (br *bitReader) readITF8() (int32, error) {
	first, err := br.readByte()
	if err != nil {
		return 0, err
	}

	if first&0x80 == 0 {
		return int32(first), nil
	}
	if first&0x40 == 0 {
		b, err := br.readByte()
		if err != nil {
			return 0, err
		}
		return int32(first&0x3f)<<8 | int32(b), nil
	}
	if first&0x20 == 0 {
		b1, err := br.readByte()
		if err != nil {
			return 0, err
		}
		b2, err := br.readByte()
		if err != nil {
			return 0, err
		}
		return int32(first&0x1f)<<16 | int32(b1)<<8 | int32(b2), nil
	}
	if first&0x10 == 0 {
		b1, err := br.readByte()
		if err != nil {
			return 0, err
		}
		b2, err := br.readByte()
		if err != nil {
			return 0, err
		}
		b3, err := br.readByte()
		if err != nil {
			return 0, err
		}
		return int32(first&0x0f)<<24 | int32(b1)<<16 | int32(b2)<<8 | int32(b3), nil
	}

	b1, err := br.readByte()
	if err != nil {
		return 0, err
	}
	b2, err := br.readByte()
	if err != nil {
		return 0, err
	}
	b3, err := br.readByte()
	if err != nil {
		return 0, err
	}
	b4, err := br.readByte()
	if err != nil {
		return 0, err
	}
	return int32(first&0x0f)<<28 | int32(b1)<<20 | int32(b2)<<12 | int32(b3)<<4 | int32(b4&0x0f), nil
}
