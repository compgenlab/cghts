package cram

import (
	"bytes"
	"fmt"
	"io"
)

// Codec IDs from the CRAM spec.
const (
	codecNull         = 0
	codecExternal     = 1
	codecGolomb       = 2 // not used in v3
	codecHuffman      = 3
	codecByteArrayLen = 4
	codecByteArrayStop = 5
	codecBeta         = 6
	codecSubexp       = 7
	codecGolombRice   = 8 // not used in v3
	codecGamma        = 9
)

// codec is the interface for all CRAM encoding codecs.
// Each codec reads values from either the core bit stream or external blocks.
type codec interface{}

// intCodec decodes integer values.
type intCodec interface {
	codec
	decodeInt(core *bitReader, external map[int32][]byte, extPos map[int32]*int) (int32, error)
}

// byteCodec decodes single byte values.
type byteCodec interface {
	codec
	decodeByte(core *bitReader, external map[int32][]byte, extPos map[int32]*int) (byte, error)
}

// byteArrayCodec decodes byte array values.
type byteArrayCodec interface {
	codec
	decodeByteArray(core *bitReader, external map[int32][]byte, extPos map[int32]*int) ([]byte, error)
}

// readEncoding reads an encoding descriptor and returns the appropriate codec.
func readEncoding(r io.Reader) (codec, error) {
	codecID, err := readITF8(r)
	if err != nil {
		return nil, fmt.Errorf("reading codec ID: %w", err)
	}

	paramLen, err := readITF8(r)
	if err != nil {
		return nil, fmt.Errorf("reading param length: %w", err)
	}

	paramData := make([]byte, paramLen)
	if paramLen > 0 {
		if _, err := io.ReadFull(r, paramData); err != nil {
			return nil, fmt.Errorf("reading params: %w", err)
		}
	}

	pr := bytes.NewReader(paramData)

	switch codecID {
	case codecNull:
		return &nullCodec{}, nil
	case codecExternal:
		return readExternalCodec(pr)
	case codecHuffman:
		return readHuffmanCodec(pr)
	case codecByteArrayLen:
		return readByteArrayLenCodec(pr)
	case codecByteArrayStop:
		return readByteArrayStopCodec(pr)
	case codecBeta:
		return readBetaCodec(pr)
	case codecSubexp:
		return readSubexpCodec(pr)
	case codecGamma:
		return readGammaCodec(pr)
	default:
		return nil, fmt.Errorf("unsupported codec: %d", codecID)
	}
}

// nullCodec: no data.
type nullCodec struct{}

// externalCodec reads values from an external block.
type externalCodec struct {
	blockID int32
}

func readExternalCodec(r io.Reader) (*externalCodec, error) {
	blockID, err := readITF8(r)
	if err != nil {
		return nil, err
	}
	return &externalCodec{blockID: blockID}, nil
}

func (c *externalCodec) decodeInt(core *bitReader, external map[int32][]byte, extPos map[int32]*int) (int32, error) {
	data, ok := external[c.blockID]
	if !ok {
		return 0, fmt.Errorf("external block %d not found", c.blockID)
	}
	pos := extPos[c.blockID]
	br := bytes.NewReader(data[*pos:])
	val, err := readITF8(br)
	if err != nil {
		return 0, err
	}
	*pos = len(data) - br.Len()
	return val, nil
}

func (c *externalCodec) decodeByte(core *bitReader, external map[int32][]byte, extPos map[int32]*int) (byte, error) {
	data, ok := external[c.blockID]
	if !ok {
		return 0, fmt.Errorf("external block %d not found", c.blockID)
	}
	pos := extPos[c.blockID]
	if *pos >= len(data) {
		return 0, io.ErrUnexpectedEOF
	}
	b := data[*pos]
	*pos++
	return b, nil
}

func (c *externalCodec) decodeByteArray(core *bitReader, external map[int32][]byte, extPos map[int32]*int) ([]byte, error) {
	// For byte arrays, external codec reads an ITF8 length followed by that many bytes.
	data, ok := external[c.blockID]
	if !ok {
		return nil, fmt.Errorf("external block %d not found", c.blockID)
	}
	pos := extPos[c.blockID]
	br := bytes.NewReader(data[*pos:])
	length, err := readITF8(br)
	if err != nil {
		return nil, err
	}
	*pos = len(data) - br.Len()
	if *pos+int(length) > len(data) {
		return nil, io.ErrUnexpectedEOF
	}
	result := make([]byte, length)
	copy(result, data[*pos:*pos+int(length)])
	*pos += int(length)
	return result, nil
}

// huffmanCodec implements canonical Huffman coding.
type huffmanCodec struct {
	alphabet   []int32
	bitLengths []int32
	// Precomputed lookup: for each code length, sorted symbols and their codes.
	codes []huffCode
}

type huffCode struct {
	symbol    int32
	code      int32
	bitLength int32
}

func readHuffmanCodec(r io.Reader) (*huffmanCodec, error) {
	alphabet, err := readITF8Array(r)
	if err != nil {
		return nil, fmt.Errorf("reading huffman alphabet: %w", err)
	}
	bitLengths, err := readITF8Array(r)
	if err != nil {
		return nil, fmt.Errorf("reading huffman bit lengths: %w", err)
	}
	if len(alphabet) != len(bitLengths) {
		return nil, fmt.Errorf("huffman: alphabet size %d != bit lengths size %d", len(alphabet), len(bitLengths))
	}

	c := &huffmanCodec{
		alphabet:   alphabet,
		bitLengths: bitLengths,
	}
	c.buildCodes()
	return c, nil
}

func (c *huffmanCodec) buildCodes() {
	// Build canonical Huffman codes.
	// Sort by bit length, then by alphabet order.
	type symLen struct {
		symbol int32
		length int32
	}
	syms := make([]symLen, len(c.alphabet))
	for i := range c.alphabet {
		syms[i] = symLen{c.alphabet[i], c.bitLengths[i]}
	}

	// Sort: by length ascending, then symbol ascending.
	for i := 1; i < len(syms); i++ {
		for j := i; j > 0 && (syms[j].length < syms[j-1].length || (syms[j].length == syms[j-1].length && syms[j].symbol < syms[j-1].symbol)); j-- {
			syms[j], syms[j-1] = syms[j-1], syms[j]
		}
	}

	c.codes = make([]huffCode, len(syms))
	code := int32(0)
	for i, sl := range syms {
		if i > 0 && sl.length > syms[i-1].length {
			code <<= uint(sl.length - syms[i-1].length)
		}
		c.codes[i] = huffCode{
			symbol:    sl.symbol,
			code:      code,
			bitLength: sl.length,
		}
		code++
	}
}

func (c *huffmanCodec) decodeInt(core *bitReader, external map[int32][]byte, extPos map[int32]*int) (int32, error) {
	// Single-symbol Huffman: bit length 0 means just return the symbol.
	if len(c.codes) == 1 && c.codes[0].bitLength == 0 {
		return c.codes[0].symbol, nil
	}

	val := int32(0)
	bits := int32(0)
	for _, hc := range c.codes {
		for bits < hc.bitLength {
			bit, err := core.readBit()
			if err != nil {
				return 0, err
			}
			val = (val << 1) | int32(bit)
			bits++
		}
		if bits == hc.bitLength && val == hc.code {
			return hc.symbol, nil
		}
	}
	return 0, fmt.Errorf("huffman: no matching code")
}

func (c *huffmanCodec) decodeByte(core *bitReader, external map[int32][]byte, extPos map[int32]*int) (byte, error) {
	v, err := c.decodeInt(core, external, extPos)
	if err != nil {
		return 0, err
	}
	return byte(v), nil
}

// betaCodec: fixed-width binary encoding.
type betaCodec struct {
	offset int32
	length int32 // number of bits
}

func readBetaCodec(r io.Reader) (*betaCodec, error) {
	offset, err := readITF8(r)
	if err != nil {
		return nil, err
	}
	length, err := readITF8(r)
	if err != nil {
		return nil, err
	}
	return &betaCodec{offset: offset, length: length}, nil
}

func (c *betaCodec) decodeInt(core *bitReader, external map[int32][]byte, extPos map[int32]*int) (int32, error) {
	val, err := core.readBits(int(c.length))
	if err != nil {
		return 0, err
	}
	return val - c.offset, nil
}

func (c *betaCodec) decodeByte(core *bitReader, external map[int32][]byte, extPos map[int32]*int) (byte, error) {
	v, err := c.decodeInt(core, external, extPos)
	if err != nil {
		return 0, err
	}
	return byte(v), nil
}

// subexpCodec: subexponential encoding.
type subexpCodec struct {
	offset int32
	k      int32
}

func readSubexpCodec(r io.Reader) (*subexpCodec, error) {
	offset, err := readITF8(r)
	if err != nil {
		return nil, err
	}
	k, err := readITF8(r)
	if err != nil {
		return nil, err
	}
	return &subexpCodec{offset: offset, k: k}, nil
}

func (c *subexpCodec) decodeInt(core *bitReader, external map[int32][]byte, extPos map[int32]*int) (int32, error) {
	// Read unary prefix.
	uLen := 0
	for {
		bit, err := core.readBit()
		if err != nil {
			return 0, err
		}
		if bit == 0 {
			break
		}
		uLen++
	}

	var val int32
	if uLen == 0 {
		// First group: k bits
		if c.k > 0 {
			v, err := core.readBits(int(c.k))
			if err != nil {
				return 0, err
			}
			val = v
		}
	} else {
		// Subsequent groups: k + uLen bits
		nbits := int(c.k) + uLen
		v, err := core.readBits(nbits)
		if err != nil {
			return 0, err
		}
		val = (int32(1) << uint(nbits)) - (int32(1) << uint(c.k)) + v
	}

	return val - c.offset, nil
}

// gammaCodec: Elias gamma encoding.
type gammaCodec struct {
	offset int32
}

func readGammaCodec(r io.Reader) (*gammaCodec, error) {
	offset, err := readITF8(r)
	if err != nil {
		return nil, err
	}
	return &gammaCodec{offset: offset}, nil
}

func (c *gammaCodec) decodeInt(core *bitReader, external map[int32][]byte, extPos map[int32]*int) (int32, error) {
	// Count leading zeros.
	nZeros := 0
	for {
		bit, err := core.readBit()
		if err != nil {
			return 0, err
		}
		if bit == 1 {
			break
		}
		nZeros++
	}

	// Read nZeros more bits.
	val := int32(1)
	for i := 0; i < nZeros; i++ {
		bit, err := core.readBit()
		if err != nil {
			return 0, err
		}
		val = (val << 1) | int32(bit)
	}

	return val - c.offset, nil
}

// byteArrayLenCodec: byte array with explicit length.
type byteArrayLenCodec struct {
	lenCodec codec
	valCodec codec
}

func readByteArrayLenCodec(r io.Reader) (*byteArrayLenCodec, error) {
	lenCodec, err := readEncoding(r)
	if err != nil {
		return nil, fmt.Errorf("reading length codec: %w", err)
	}
	valCodec, err := readEncoding(r)
	if err != nil {
		return nil, fmt.Errorf("reading value codec: %w", err)
	}
	return &byteArrayLenCodec{lenCodec: lenCodec, valCodec: valCodec}, nil
}

func (c *byteArrayLenCodec) decodeByteArray(core *bitReader, external map[int32][]byte, extPos map[int32]*int) ([]byte, error) {
	// Decode length.
	lc, ok := c.lenCodec.(intCodec)
	if !ok {
		return nil, fmt.Errorf("byteArrayLen: length codec is not intCodec")
	}
	length, err := lc.decodeInt(core, external, extPos)
	if err != nil {
		return nil, fmt.Errorf("byteArrayLen: decoding length: %w", err)
	}

	// Decode value bytes.
	result := make([]byte, length)
	if vc, ok := c.valCodec.(*externalCodec); ok {
		// Optimization: read directly from external block.
		data, exists := external[vc.blockID]
		if !exists {
			return nil, fmt.Errorf("external block %d not found", vc.blockID)
		}
		pos := extPos[vc.blockID]
		if *pos+int(length) > len(data) {
			return nil, io.ErrUnexpectedEOF
		}
		copy(result, data[*pos:*pos+int(length)])
		*pos += int(length)
	} else if vc, ok := c.valCodec.(byteCodec); ok {
		for i := int32(0); i < length; i++ {
			result[i], err = vc.decodeByte(core, external, extPos)
			if err != nil {
				return nil, err
			}
		}
	} else {
		return nil, fmt.Errorf("byteArrayLen: value codec cannot decode bytes")
	}

	return result, nil
}

// byteArrayStopCodec: byte array terminated by a stop byte.
type byteArrayStopCodec struct {
	stopByte byte
	blockID  int32
}

func readByteArrayStopCodec(r io.Reader) (*byteArrayStopCodec, error) {
	var stop [1]byte
	if _, err := io.ReadFull(r, stop[:]); err != nil {
		return nil, err
	}
	blockID, err := readITF8(r)
	if err != nil {
		return nil, err
	}
	return &byteArrayStopCodec{stopByte: stop[0], blockID: blockID}, nil
}

func (c *byteArrayStopCodec) decodeByteArray(core *bitReader, external map[int32][]byte, extPos map[int32]*int) ([]byte, error) {
	data, ok := external[c.blockID]
	if !ok {
		return nil, fmt.Errorf("external block %d not found", c.blockID)
	}
	pos := extPos[c.blockID]
	start := *pos
	for *pos < len(data) {
		if data[*pos] == c.stopByte {
			result := make([]byte, *pos-start)
			copy(result, data[start:*pos])
			*pos++ // skip stop byte
			return result, nil
		}
		*pos++
	}
	return nil, fmt.Errorf("byteArrayStop: stop byte 0x%02x not found", c.stopByte)
}
