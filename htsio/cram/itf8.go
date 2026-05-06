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

// readITF8Array reads an ITF8 length followed by that many ITF8 values.
func readITF8Array(r io.Reader) ([]int32, error) {
	n, err := readITF8(r)
	if err != nil {
		return nil, fmt.Errorf("reading array length: %w", err)
	}
	if n < 0 {
		return nil, fmt.Errorf("negative array length: %d", n)
	}
	vals := make([]int32, n)
	for i := int32(0); i < n; i++ {
		vals[i], err = readITF8(r)
		if err != nil {
			return nil, fmt.Errorf("reading array element %d: %w", i, err)
		}
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
