package codec

import (
	"encoding/binary"
	"fmt"
)

// rANS Nx16 decoder (CRAM v3.1).
// Supports order-0, order-1, with optional PACK, RLE, STRIPE, and CAT transforms.
// Based on the htscodecs reference implementation.

// Nx16 flags.
const (
	ransOrderMask   = 0x01
	ransOrderX32    = 0x04 // 32-way interleaving (unused for decode, same loop)
	ransOrderStripe = 0x08
	ransOrderNoSize = 0x10
	ransOrderCat    = 0x20
	ransOrderRLE    = 0x40
	ransOrderPack   = 0x80
)

// Nx16 constants.
const (
	ransNx16L       = 1 << 15 // renormalization threshold
	ransNx16TFShift = 12
	ransNx16TotFreq = 1 << ransNx16TFShift
)

func DecodeRansNx16(data []byte) ([]byte, error) {
	return decodeRansNx16WithSize(data, 0)
}

// decodeRansNx16WithSize decodes rANS Nx16 data. If expectedSize > 0, it is used
// when the NOSIZE flag is set (e.g., for STRIPE sub-streams).
func decodeRansNx16WithSize(data []byte, expectedSize uint32) ([]byte, error) {
	if len(data) == 0 {
		return nil, fmt.Errorf("rans_nx16: empty data")
	}

	// Check for STRIPE first (before consuming the flags byte).
	if data[0]&ransOrderStripe != 0 {
		return decodeRansNx16Stripe(data)
	}

	flags := data[0]
	data = data[1:]

	doPack := flags&ransOrderPack != 0
	doRLE := flags&ransOrderRLE != 0
	doCat := flags&ransOrderCat != 0
	noSize := flags&ransOrderNoSize != 0
	order := int(flags & ransOrderMask)

	// Read uncompressed size.
	var outSize uint32
	if !noSize {
		n, val := varGetU32(data)
		if n == 0 {
			return nil, fmt.Errorf("rans_nx16: truncated size")
		}
		outSize = val
		data = data[n:]
	} else {
		outSize = expectedSize
	}

	// Read PACK metadata.
	var packMap [16]byte
	var nPackedSym int
	var unpackedSize uint32
	if doPack {
		n, nsym, mp := readPackMeta(data)
		if n == 0 {
			return nil, fmt.Errorf("rans_nx16: truncated pack metadata")
		}
		packMap = mp
		nPackedSym = nsym
		unpackedSize = outSize
		data = data[n:]

		// Read the packed size (size after packing, before rANS).
		n2, packedSz := varGetU32(data)
		if n2 == 0 {
			return nil, fmt.Errorf("rans_nx16: truncated packed size")
		}
		outSize = packedSz
		data = data[n2:]
	}

	// Read RLE metadata.
	var rleMeta []byte
	var rleSyms []byte
	var rleNSyms int
	var rleLen uint32
	if doRLE {
		n1, uMetaSize := varGetU32(data)
		if n1 == 0 {
			return nil, fmt.Errorf("rans_nx16: truncated rle meta size")
		}
		n2, rl := varGetU32(data[n1:])
		if n2 == 0 {
			return nil, fmt.Errorf("rans_nx16: truncated rle len")
		}
		rleLen = rl

		if uMetaSize&1 != 0 {
			// Uncompressed meta.
			cMetaSize := uMetaSize / 2
			off := n1 + n2
			if uint32(len(data)-off) < cMetaSize {
				return nil, fmt.Errorf("rans_nx16: truncated rle meta")
			}
			rleMeta = data[off : off+int(cMetaSize)]
			data = data[off+int(cMetaSize):]
		} else {
			// Compressed meta (order-0 rANS).
			n3, cMetaSize := varGetU32(data[n1+n2:])
			if n3 == 0 {
				return nil, fmt.Errorf("rans_nx16: truncated rle cmeta size")
			}
			off := n1 + n2 + n3
			uMetaSize /= 2
			meta, err := decodeRansNx16Order0(data[off:off+int(cMetaSize)], uMetaSize, 4)
			if err != nil {
				return nil, fmt.Errorf("rans_nx16: decompressing rle meta: %w", err)
			}
			rleMeta = meta
			data = data[off+int(cMetaSize):]
		}

		// Parse RLE symbols from meta.
		if len(rleMeta) < 1 {
			return nil, fmt.Errorf("rans_nx16: empty rle meta")
		}
		rleNSyms = int(rleMeta[0])
		if rleNSyms == 0 {
			rleNSyms = 256
		}
		if len(rleMeta) < 1+rleNSyms {
			return nil, fmt.Errorf("rans_nx16: truncated rle syms")
		}
		rleSyms = rleMeta[1 : 1+rleNSyms]
		rleMeta = rleMeta[1+rleNSyms:]
		outSize = rleLen
	}

	// Determine interleaving width: 32 if X32 flag, else 4.
	nx := 4
	if flags&ransOrderX32 != 0 {
		nx = 32
	}

	// Core decode: in -> tmp1
	var tmp1 []byte
	if doCat {
		if int(outSize) > len(data) {
			return nil, fmt.Errorf("rans_nx16: CAT size %d exceeds data %d", outSize, len(data))
		}
		tmp1 = make([]byte, outSize)
		copy(tmp1, data[:outSize])
	} else if len(data) > 0 {
		var err error
		if order == 0 {
			tmp1, err = decodeRansNx16Order0(data, outSize, nx)
		} else {
			tmp1, err = decodeRansNx16Order1(data, outSize, nx)
		}
		if err != nil {
			return nil, err
		}
	} else {
		tmp1 = []byte{}
	}

	// Un-RLE: tmp1 -> tmp2
	tmp2 := tmp1
	if doRLE {
		var err error
		tmp2, err = rleDecodeNx16(tmp1, rleMeta, rleSyms)
		if err != nil {
			return nil, fmt.Errorf("rans_nx16: rle decode: %w", err)
		}
	}

	// Unpack: tmp2 -> tmp3
	tmp3 := tmp2
	if doPack {
		var err error
		tmp3, err = unpackNx16(tmp2, unpackedSize, nPackedSym, packMap)
		if err != nil {
			return nil, fmt.Errorf("rans_nx16: unpack: %w", err)
		}
	}

	return tmp3, nil
}

// decodeRansNx16Stripe handles the STRIPE transform.
func decodeRansNx16Stripe(data []byte) ([]byte, error) {
	if len(data) < 2 {
		return nil, fmt.Errorf("rans_nx16_stripe: truncated")
	}

	pos := 1 // skip flags byte

	// Read uncompressed length.
	n, ulen := varGetU32(data[pos:])
	if n == 0 {
		return nil, fmt.Errorf("rans_nx16_stripe: truncated ulen")
	}
	pos += n

	if pos >= len(data) {
		return nil, fmt.Errorf("rans_nx16_stripe: truncated N")
	}
	N := int(data[pos])
	pos++
	if N < 1 {
		return nil, fmt.Errorf("rans_nx16_stripe: N must be >= 1")
	}

	// Read compressed lengths for each stripe.
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
			return nil, fmt.Errorf("rans_nx16_stripe: truncated clen %d", i)
		}
		clenN[i] = int(cl)
		pos += n
	}

	// Decompress each stripe.
	outN := make([]byte, ulen)
	for i := 0; i < N; i++ {
		if pos+clenN[i] > len(data) {
			return nil, fmt.Errorf("rans_nx16_stripe: stripe %d overflows input", i)
		}
		decoded, err := decodeRansNx16WithSize(data[pos:pos+clenN[i]], uint32(ulenN[i]))
		if err != nil {
			return nil, fmt.Errorf("rans_nx16_stripe: stripe %d: %w", i, err)
		}
		if len(decoded) != ulenN[i] {
			return nil, fmt.Errorf("rans_nx16_stripe: stripe %d size mismatch: got %d, want %d", i, len(decoded), ulenN[i])
		}
		copy(outN[idxN[i]:], decoded)
		pos += clenN[i]
	}

	// Unstripe: interleave the N streams.
	out := make([]byte, ulen)
	for i := 0; i < int(ulen); i++ {
		stripe := i % N
		idx := idxN[stripe] + i/N
		out[i] = outN[idx]
	}

	return out, nil
}

func boolToInt(b bool) int {
	if b {
		return 1
	}
	return 0
}

// decodeRansNx16Order0 decodes order-0 rANS with 16-bit renormalization.
// nx is the interleaving width (4 or 32).
func decodeRansNx16Order0(data []byte, outSize uint32, nx int) ([]byte, error) {
	if len(data) < 4 {
		return nil, fmt.Errorf("rans_nx16_o0: data too short")
	}

	// Read frequency table.
	F, fsum, pos, err := nx16DecodeFreq(data)
	if err != nil {
		return nil, err
	}
	data = data[pos:]

	// Normalize frequencies.
	nx16NormaliseFreqShift(&F, fsum, ransNx16TotFreq)

	// Build lookup tables.
	var lutSym [ransNx16TotFreq]byte
	var lutFreq [ransNx16TotFreq]uint16
	var lutBase [ransNx16TotFreq]uint16
	x := uint32(0)
	for j := 0; j < 256; j++ {
		if F[j] > 0 {
			for y := uint32(0); y < F[j]; y++ {
				lutSym[x+y] = byte(j)
				lutFreq[x+y] = uint16(F[j])
				lutBase[x+y] = uint16(y)
			}
			x += F[j]
		}
	}
	if x != ransNx16TotFreq {
		return nil, fmt.Errorf("rans_nx16_o0: frequency sum %d != %d", x, ransNx16TotFreq)
	}

	// Initialize NX rANS states.
	stateBytes := nx * 4
	if len(data) < stateBytes {
		return nil, fmt.Errorf("rans_nx16_o0: not enough data for %d states", nx)
	}
	state := make([]uint32, nx)
	for i := 0; i < nx; i++ {
		state[i] = binary.LittleEndian.Uint32(data[i*4:])
	}
	data = data[stateBytes:]

	// Decode with round-robin interleaving.
	out := make([]byte, outSize)
	dataPos := 0

	for i := 0; i < int(outSize); i++ {
		si := i % nx
		m := state[si] & (ransNx16TotFreq - 1)
		out[i] = lutSym[m]
		state[si] = uint32(lutFreq[m])*(state[si]>>ransNx16TFShift) + uint32(lutBase[m])

		// 16-bit renormalization.
		if state[si] < ransNx16L && dataPos+1 < len(data) {
			state[si] = (state[si] << 16) | uint32(binary.LittleEndian.Uint16(data[dataPos:]))
			dataPos += 2
		}
	}

	return out, nil
}

// decodeRansNx16Order1 decodes order-1 rANS with 16-bit renormalization.
// nx is the interleaving width (4 or 32).
func decodeRansNx16Order1(data []byte, outSize uint32, nx int) ([]byte, error) {
	if len(data) < 4 {
		return nil, fmt.Errorf("rans_nx16_o1: data too short")
	}

	// First byte: upper 4 bits = TF_SHIFT, bit 0 = compressed freq table.
	firstByte := data[0]
	tfShift := uint(firstByte >> 4)
	if tfShift < 1 || tfShift > 12 {
		tfShift = 12
	}
	totFreq := uint32(1) << tfShift
	compressedFreqTable := firstByte&1 != 0

	var freqData []byte
	pos := 1 // skip first byte
	if compressedFreqTable {
		n1, uSize := varGetU32(data[pos:])
		pos += n1
		n2, cSize := varGetU32(data[pos:])
		pos += n2
		if pos+int(cSize) > len(data) {
			return nil, fmt.Errorf("rans_nx16_o1: truncated compressed freq table")
		}
		var err error
		freqData, err = decodeRansNx16Order0(data[pos:pos+int(cSize)], uSize, 4)
		if err != nil {
			return nil, fmt.Errorf("rans_nx16_o1: decompressing freq table: %w", err)
		}
		pos += int(cSize)
		data = data[pos:]
	} else {
		freqData = data[pos:]
	}

	// Read order-0 alphabet (which symbols are present).
	var F0 [256]uint32
	apos := nx16DecodeAlphabet(freqData, &F0)
	if apos == 0 {
		return nil, fmt.Errorf("rans_nx16_o1: failed to decode alphabet")
	}
	freqData = freqData[apos:]

	// Read per-context frequency tables.
	type lutEntry struct {
		sym  byte
		freq uint16
		base uint16
	}
	// Use slices since totFreq is determined at runtime.
	syms := make([][256]ransDecSymbol, 256)
	luts := make([][]lutEntry, 256)
	for i := range luts {
		luts[i] = make([]lutEntry, totFreq)
	}

	for ctx := 0; ctx < 256; ctx++ {
		if F0[ctx] == 0 {
			continue
		}

		// Read frequencies for this context using decode_freq_d format.
		var F [256]uint32
		var total uint32
		n := nx16DecodeFreqD(freqData, &F0, &F, &total)
		if n == 0 {
			return nil, fmt.Errorf("rans_nx16_o1: failed to decode freq for ctx %d", ctx)
		}
		freqData = freqData[n:]

		// Normalize.
		nx16NormaliseFreqShift(&F, total, totFreq)

		// Build LUT for this context.
		x := uint32(0)
		for j := 0; j < 256; j++ {
			if F[j] > 0 {
				syms[ctx][j] = ransDecSymbol{cumFreq: x, freq: F[j]}
				for y := uint32(0); y < F[j] && x+y < totFreq; y++ {
					luts[ctx][x+y] = lutEntry{
						sym:  byte(j),
						freq: uint16(F[j]),
						base: uint16(y),
					}
				}
				x += F[j]
			}
		}
		if x != totFreq {
			return nil, fmt.Errorf("rans_nx16_o1: ctx %d freq sum %d != %d", ctx, x, totFreq)
		}
	}

	if !compressedFreqTable {
		// If the freq table was inline, advance data past it.
		consumed := len(data) - 1 - len(freqData) // -1 for the first byte already consumed
		data = data[1+consumed:]
	}

	// Initialize NX rANS states.
	stateBytes := nx * 4
	if len(data) < stateBytes {
		return nil, fmt.Errorf("rans_nx16_o1: not enough data for %d states", nx)
	}
	state := make([]uint32, nx)
	for i := 0; i < nx; i++ {
		state[i] = binary.LittleEndian.Uint32(data[i*4:])
	}
	data = data[stateBytes:]

	// Decode NX-way split: each state decodes 1/NX of the output.
	out := make([]byte, outSize)
	dataPos := 0
	iszN := int(outSize) / nx
	lastSym := make([]byte, nx)
	iN := make([]int, nx)
	for z := 0; z < nx; z++ {
		iN[z] = z * iszN
	}

	mask := totFreq - 1
	for iN[0] < iszN {
		for z := 0; z < nx; z++ {
			ctx := lastSym[z]
			m := state[z] & mask
			entry := luts[ctx][m]
			out[iN[z]] = entry.sym
			lastSym[z] = entry.sym
			iN[z]++

			state[z] = uint32(entry.freq)*(state[z]>>tfShift) + uint32(entry.base)

			// 16-bit renormalization.
			if state[z] < ransNx16L && dataPos+1 < len(data) {
				state[z] = (state[z] << 16) | uint32(binary.LittleEndian.Uint16(data[dataPos:]))
				dataPos += 2
			}
		}
	}

	// Remainder: last state handles any leftover bytes.
	last := nx - 1
	for iN[last] < int(outSize) {
		ctx := lastSym[last]
		m := state[last] & mask
		entry := luts[ctx][m]
		out[iN[last]] = entry.sym
		lastSym[last] = entry.sym
		iN[last]++

		state[last] = uint32(entry.freq)*(state[last]>>tfShift) + uint32(entry.base)
		if state[last] < ransNx16L && dataPos+1 < len(data) {
			state[last] = (state[last] << 16) | uint32(binary.LittleEndian.Uint16(data[dataPos:]))
			dataPos += 2
		}
	}

	return out, nil
}

// varGetU32 reads a variable-length uint32 (7-bit chunks, MSB continuation).
// Returns (bytes consumed, value).
func varGetU32(data []byte) (int, uint32) {
	if len(data) == 0 {
		return 0, 0
	}
	val := uint32(0)
	n := 0
	for n < len(data) && n < 5 {
		c := data[n]
		n++
		val = (val << 7) | uint32(c&0x7f)
		if c&0x80 == 0 {
			break
		}
	}
	return n, val
}

// nx16DecodeAlphabet reads the RLE-compressed alphabet.
// Sets F[j] = 1 for each symbol j present.
func nx16DecodeAlphabet(data []byte, F *[256]uint32) int {
	if len(data) == 0 {
		return 0
	}

	pos := 0
	rle := 0
	j := int(data[pos])
	pos++

	for {
		if j > 255 {
			return 0
		}
		F[j] = 1

		if rle > 0 {
			rle--
			j++
		} else if pos < len(data) && j+1 == int(data[pos]) {
			j = int(data[pos])
			pos++
			if pos >= len(data) {
				break
			}
			rle = int(data[pos])
			pos++
		} else {
			if pos >= len(data) {
				break
			}
			j = int(data[pos])
			pos++
			if j == 0 {
				break
			}
		}
	}

	return pos
}

// nx16DecodeFreq reads an order-0 frequency table.
func nx16DecodeFreq(data []byte) ([256]uint32, uint32, int, error) {
	var F [256]uint32

	apos := nx16DecodeAlphabet(data, &F)
	if apos == 0 {
		return F, 0, 0, fmt.Errorf("rans_nx16: failed to decode alphabet")
	}

	pos := apos
	tot := uint32(0)
	for j := 0; j < 256; j++ {
		if F[j] != 0 {
			n, val := varGetU32(data[pos:])
			if n == 0 {
				return F, 0, 0, fmt.Errorf("rans_nx16: truncated freq at sym %d", j)
			}
			F[j] = val
			tot += val
			pos += n
		}
	}

	return F, tot, pos, nil
}

// nx16DecodeFreqD reads an order-1 per-context frequency table (delta/compact format).
func nx16DecodeFreqD(data []byte, F0 *[256]uint32, F *[256]uint32, total *uint32) int {
	pos := 0
	dz := 0
	T := uint32(0)

	for j := 0; j < 256; j++ {
		if F0[j] == 0 {
			continue
		}

		var f uint32
		if dz > 0 {
			f = 0
			dz--
		} else {
			if pos >= len(data) {
				return 0
			}
			n, val := varGetU32(data[pos:])
			if n == 0 {
				return 0
			}
			f = val
			pos += n
			if f == 0 {
				if pos >= len(data) {
					return 0
				}
				dz = int(data[pos])
				pos++
			}
		}
		F[j] = f
		T += f
	}

	*total = T
	return pos
}

// nx16NormaliseFreqShift normalises frequencies by shifting up to fill maxTot.
// This only works when size is a power-of-2 factor of maxTot.
func nx16NormaliseFreqShift(F *[256]uint32, size, maxTot uint32) {
	if size == 0 || size == maxTot {
		return
	}
	if size > maxTot {
		return // already larger, can't shift
	}
	shift := uint(0)
	s := size
	for s < maxTot {
		s *= 2
		shift++
	}
	if s != maxTot {
		return // not an exact power-of-2 ratio, leave as is
	}
	for i := 0; i < 256; i++ {
		F[i] <<= shift
	}
}

// readPackMeta reads the PACK metadata: nsyms and symbol map.
// Returns (bytes consumed, npacked_sym, map).
func readPackMeta(data []byte) (int, int, [16]byte) {
	var mp [16]byte
	if len(data) < 1 {
		return 0, 0, mp
	}

	pos := 0
	nsym := int(data[pos])
	pos++
	if nsym == 0 {
		nsym = 256
	}

	if nsym <= 1 {
		// 0 bits per symbol — constant value.
		if pos < len(data) {
			mp[0] = data[pos]
			pos++
		}
		return pos, 0, mp
	}

	// Read the symbol map.
	for i := 0; i < nsym && pos < len(data); i++ {
		if i < 16 {
			mp[i] = data[pos]
		}
		pos++
	}

	// Determine packed_sym value for unpack.
	var nPacked int
	switch {
	case nsym <= 2:
		nPacked = 8 // 1 bit per sym → 8 syms per byte
	case nsym <= 4:
		nPacked = 4 // 2 bits per sym → 4 syms per byte
	case nsym <= 16:
		nPacked = 2 // 4 bits per sym → 2 syms per byte
	default:
		nPacked = 1 // no packing
	}

	return pos, nPacked, mp
}

// unpackNx16 reverses the PACK transform.
func unpackNx16(data []byte, outSize uint32, nPackedSym int, mp [16]byte) ([]byte, error) {
	out := make([]byte, outSize)

	if nPackedSym == 0 {
		// All same symbol.
		for i := range out {
			out[i] = mp[0]
		}
		return out, nil
	}

	if nPackedSym == 1 {
		// No packing, just copy.
		copy(out, data)
		return out, nil
	}

	j := 0
	i := 0

	switch nPackedSym {
	case 8: // 1 bit per symbol, 8 per byte
		for j < len(data) && i < int(outSize) {
			c := data[j]
			j++
			for b := 0; b < 8 && i < int(outSize); b++ {
				out[i] = mp[c&1]
				c >>= 1
				i++
			}
		}
	case 4: // 2 bits per symbol, 4 per byte
		for j < len(data) && i < int(outSize) {
			c := data[j]
			j++
			for b := 0; b < 4 && i < int(outSize); b++ {
				out[i] = mp[c&3]
				c >>= 2
				i++
			}
		}
	case 2: // 4 bits per symbol, 2 per byte
		for j < len(data) && i < int(outSize) {
			c := data[j]
			j++
			out[i] = mp[c&0xf]
			i++
			if i < int(outSize) {
				out[i] = mp[(c>>4)&0xf]
				i++
			}
		}
	default:
		return nil, fmt.Errorf("unpack: unsupported nPackedSym=%d", nPackedSym)
	}

	return out, nil
}

// rleDecodeNx16 reverses the RLE transform.
// lit is the literal data, run is the run-length data, rleSyms lists
// which symbols can be RLE-encoded.
func rleDecodeNx16(lit []byte, run []byte, rleSyms []byte) ([]byte, error) {
	// Build lookup of RLE-eligible symbols.
	var isRLE [256]bool
	for _, s := range rleSyms {
		isRLE[s] = true
	}

	// Estimate output size (we'll grow if needed).
	out := make([]byte, 0, len(lit)*2)
	runPos := 0

	for _, b := range lit {
		if isRLE[b] {
			rlen := uint32(0)
			if runPos < len(run) {
				n, val := varGetU32(run[runPos:])
				rlen = val
				runPos += n
			}
			for k := uint32(0); k <= rlen; k++ {
				out = append(out, b)
			}
		} else {
			out = append(out, b)
		}
	}

	return out, nil
}
