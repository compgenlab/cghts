package codec

import (
	"bytes"
	"compress/bzip2"
	"fmt"
	"io"
)

// Adaptive arithmetic coder (CRAM v3.1 method 6).
// Supports order-0, order-1, with optional PACK, RLE, STRIPE, CAT, and EXT transforms.
// Based on the htscodecs arith_dynamic.c reference implementation.

// Arith flags — same bit positions as rANS Nx16 except X_EXT replaces X32.
const (
	arithOrderMask   = 0x03
	arithOrderExt    = 0x04 // external compression (bzip2)
	arithOrderStripe = 0x08
	arithOrderNoSize = 0x10
	arithOrderCat    = 0x20
	arithOrderRLE    = 0x40
	arithOrderPack   = 0x80
)

const arithMaxRun = 4

func DecodeArithDynamic(data []byte) ([]byte, error) {
	return decodeArithDynamicWithSize(data, 0)
}

func decodeArithDynamicWithSize(data []byte, expectedSize uint32) ([]byte, error) {
	if len(data) == 0 {
		return nil, fmt.Errorf("arith: empty data")
	}

	// Check for STRIPE first.
	if data[0]&arithOrderStripe != 0 {
		return decodeArithStripe(data)
	}

	flags := data[0]
	data = data[1:]

	doPack := flags&arithOrderPack != 0
	doRLE := flags&arithOrderRLE != 0
	doCat := flags&arithOrderCat != 0
	noSize := flags&arithOrderNoSize != 0
	doExt := flags&arithOrderExt != 0
	order := int(flags & arithOrderMask)

	// Read uncompressed size.
	var outSize uint32
	if !noSize {
		n, val := varGetU32(data)
		if n == 0 {
			return nil, fmt.Errorf("arith: truncated size")
		}
		outSize = val
		data = data[n:]
	} else {
		outSize = expectedSize
	}

	// Read PACK metadata (same format as rANS Nx16).
	var packMap [16]byte
	var nPackedSym int
	var unpackedSize uint32
	if doPack {
		n, nsym, mp := readPackMeta(data)
		if n == 0 {
			return nil, fmt.Errorf("arith: truncated pack metadata")
		}
		packMap = mp
		nPackedSym = nsym
		unpackedSize = outSize
		data = data[n:]

		n2, packedSz := varGetU32(data)
		if n2 == 0 {
			return nil, fmt.Errorf("arith: truncated packed size")
		}
		outSize = packedSz
		data = data[n2:]
	}

	// Core decode.
	var tmp1 []byte
	if doCat {
		if int(outSize) > len(data) {
			return nil, fmt.Errorf("arith: CAT size %d exceeds data %d", outSize, len(data))
		}
		tmp1 = make([]byte, outSize)
		copy(tmp1, data[:outSize])
	} else if doExt {
		// External compression (bzip2).
		var err error
		tmp1, err = io.ReadAll(bzip2.NewReader(bytes.NewReader(data)))
		if err != nil {
			return nil, fmt.Errorf("arith: ext decompress: %w", err)
		}
	} else if len(data) > 0 {
		var err error
		if doRLE {
			if order == 1 {
				tmp1, err = arithUncompressO1RLE(data, outSize)
			} else {
				tmp1, err = arithUncompressO0RLE(data, outSize)
			}
		} else {
			if order == 1 {
				tmp1, err = arithUncompressO1(data, outSize)
			} else {
				tmp1, err = arithUncompressO0(data, outSize)
			}
		}
		if err != nil {
			return nil, err
		}
	} else {
		tmp1 = []byte{}
	}

	// Unpack.
	if doPack {
		var err error
		tmp1, err = unpackNx16(tmp1, unpackedSize, nPackedSym, packMap)
		if err != nil {
			return nil, fmt.Errorf("arith: unpack: %w", err)
		}
	}

	return tmp1, nil
}

func decodeArithStripe(data []byte) ([]byte, error) {
	if len(data) < 2 {
		return nil, fmt.Errorf("arith_stripe: truncated")
	}

	pos := 1 // skip flags byte

	n, ulen := varGetU32(data[pos:])
	if n == 0 {
		return nil, fmt.Errorf("arith_stripe: truncated ulen")
	}
	pos += n

	if pos >= len(data) {
		return nil, fmt.Errorf("arith_stripe: truncated N")
	}
	N := int(data[pos])
	pos++
	if N < 1 {
		return nil, fmt.Errorf("arith_stripe: N must be >= 1")
	}

	clenN := make([]int, N)
	ulenN := make([]int, N)
	idxN := make([]int, N)
	for i := 0; i < N; i++ {
		ulenN[i] = int(ulen)/N + boolToInt(int(ulen)%N > i)
		if i > 0 {
			idxN[i] = idxN[i-1] + ulenN[i-1]
		}
		n, cl := varGetU32(data[pos:])
		if n == 0 {
			return nil, fmt.Errorf("arith_stripe: truncated clen %d", i)
		}
		clenN[i] = int(cl)
		pos += n
	}

	outN := make([]byte, ulen)
	for i := 0; i < N; i++ {
		if pos+clenN[i] > len(data) {
			return nil, fmt.Errorf("arith_stripe: stripe %d overflows input", i)
		}
		decoded, err := decodeArithDynamicWithSize(data[pos:pos+clenN[i]], uint32(ulenN[i]))
		if err != nil {
			return nil, fmt.Errorf("arith_stripe: stripe %d: %w", i, err)
		}
		copy(outN[idxN[i]:], decoded)
		pos += clenN[i]
	}

	// Unstripe: interleave.
	out := make([]byte, ulen)
	for i := 0; i < int(ulen); i++ {
		stripe := i % N
		idx := idxN[stripe] + i/N
		out[i] = outN[idx]
	}

	return out, nil
}

func arithUncompressO0(data []byte, outSize uint32) ([]byte, error) {
	if len(data) < 1 {
		return nil, fmt.Errorf("arith_o0: data too short")
	}

	m := int(data[0])
	if m == 0 {
		m = 256
	}

	model := newSimpleModel(256, m)
	rc := newRangeDecoder(data[1:])

	out := make([]byte, outSize)
	for i := uint32(0); i < outSize; i++ {
		out[i] = byte(model.decodeSymbol(rc))
	}

	if !rc.finish() {
		return nil, fmt.Errorf("arith_o0: decode error")
	}
	return out, nil
}

func arithUncompressO1(data []byte, outSize uint32) ([]byte, error) {
	if len(data) < 1 {
		return nil, fmt.Errorf("arith_o1: data too short")
	}

	m := int(data[0])
	if m == 0 {
		m = 256
	}

	models := make([]*simpleModel, 256)
	for i := 0; i < 256; i++ {
		models[i] = newSimpleModel(256, m)
	}

	rc := newRangeDecoder(data[1:])

	out := make([]byte, outSize)
	var last byte
	for i := uint32(0); i < outSize; i++ {
		out[i] = byte(models[last].decodeSymbol(rc))
		last = out[i]
	}

	if !rc.finish() {
		return nil, fmt.Errorf("arith_o1: decode error")
	}
	return out, nil
}

func arithUncompressO0RLE(data []byte, outSize uint32) ([]byte, error) {
	if len(data) < 1 {
		return nil, fmt.Errorf("arith_o0_rle: data too short")
	}

	m := int(data[0])
	if m == 0 {
		m = 256
	}

	byteModel := newSimpleModel(256, m)

	// Run models: 258 models (one per symbol + 2 extra contexts), each with maxSym = arithMaxRun.
	const nRunSym = 258
	runModels := make([]*simpleModel, nRunSym)
	for i := 0; i < nRunSym; i++ {
		runModels[i] = newSimpleModel(nRunSym, arithMaxRun)
	}

	rc := newRangeDecoder(data[1:])

	out := make([]byte, outSize)
	for i := uint32(0); i < outSize; i++ {
		last := byte(byteModel.decodeSymbol(rc))
		out[i] = last

		run := 0
		rctx := int(last)
		for {
			r := int(runModels[rctx].decodeSymbol(rc))
			if rctx == int(last) {
				rctx = 256
			} else if rctx < nRunSym-1 {
				rctx++
			}
			run += r
			if r != arithMaxRun-1 || run >= int(outSize) {
				break
			}
		}
		for run > 0 && i+1 < outSize {
			i++
			out[i] = last
			run--
		}
	}

	if !rc.finish() {
		return nil, fmt.Errorf("arith_o0_rle: decode error")
	}
	return out, nil
}

func arithUncompressO1RLE(data []byte, outSize uint32) ([]byte, error) {
	if len(data) < 1 {
		return nil, fmt.Errorf("arith_o1_rle: data too short")
	}

	m := int(data[0])
	if m == 0 {
		m = 256
	}

	byteModels := make([]*simpleModel, 256)
	for i := 0; i < 256; i++ {
		byteModels[i] = newSimpleModel(256, m)
	}

	const nRunSym = 258
	runModels := make([]*simpleModel, nRunSym)
	for i := 0; i < nRunSym; i++ {
		runModels[i] = newSimpleModel(nRunSym, arithMaxRun)
	}

	rc := newRangeDecoder(data[1:])

	out := make([]byte, outSize)
	var ctx byte
	for i := uint32(0); i < outSize; i++ {
		last := byte(byteModels[ctx].decodeSymbol(rc))
		out[i] = last

		run := 0
		rctx := int(last)
		for {
			r := int(runModels[rctx].decodeSymbol(rc))
			if rctx == int(last) {
				rctx = 256
			} else if rctx < nRunSym-1 {
				rctx++
			}
			run += r
			if r != arithMaxRun-1 || run >= int(outSize) {
				break
			}
		}
		for run > 0 && i+1 < outSize {
			i++
			out[i] = last
		}
		ctx = last
	}

	if !rc.finish() {
		return nil, fmt.Errorf("arith_o1_rle: decode error")
	}
	return out, nil
}

