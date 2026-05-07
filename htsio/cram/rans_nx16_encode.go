package cram

import (
	"encoding/binary"
	"sort"
)

// rANS Nx16 encoder for CRAM v3.1.
// Supports order-0 with optional PACK transform.

// varPutU32 encodes a uint32 as a 7-bit varint, MSB first.
func varPutU32(buf []byte, val uint32) int {
	n := 0
	if val >= (1 << 28) {
		buf[n] = byte((val>>28)&0x7f) | 0x80
		n++
	}
	if val >= (1 << 21) {
		buf[n] = byte((val>>21)&0x7f) | 0x80
		n++
	}
	if val >= (1 << 14) {
		buf[n] = byte((val>>14)&0x7f) | 0x80
		n++
	}
	if val >= (1 << 7) {
		buf[n] = byte((val>>7)&0x7f) | 0x80
		n++
	}
	buf[n] = byte(val & 0x7f)
	n++
	return n
}

// varPutU32Slice appends a 7-bit varint to a slice.
func varPutU32Slice(out []byte, val uint32) []byte {
	var buf [5]byte
	n := varPutU32(buf[:], val)
	return append(out, buf[:n]...)
}

// encodeRansNx16 encodes data using rANS Nx16 codec.
// Returns compressed data including flags byte prefix.
func encodeRansNx16(data []byte) []byte {
	if len(data) == 0 {
		return []byte{0, 0}
	}

	// Count distinct symbols.
	var seen [256]bool
	nsyms := 0
	for _, b := range data {
		if !seen[b] {
			seen[b] = true
			nsyms++
		}
	}

	// Try PACK transform if few symbols.
	if nsyms <= 16 {
		packed, packMeta := packNx16Encode(data, nsyms)
		// Encode packed data with order-0.
		encoded := encodeRansNx16Order0Core(packed)
		flags := byte(ransOrderPack)

		var out []byte
		out = append(out, flags)
		out = varPutU32Slice(out, uint32(len(data))) // uncompressed size
		out = append(out, packMeta...)
		out = varPutU32Slice(out, uint32(len(packed))) // packed size
		out = append(out, encoded...)

		// Also try without PACK and pick the smaller.
		nopack := encodeRansNx16Order0(data)
		if len(nopack) < len(out) {
			return nopack
		}
		return out
	}

	return encodeRansNx16Order0(data)
}

// encodeRansNx16Order0 encodes data with order-0 rANS Nx16 (no transforms).
func encodeRansNx16Order0(data []byte) []byte {
	flags := byte(0) // order-0, no transforms
	encoded := encodeRansNx16Order0Core(data)

	var out []byte
	out = append(out, flags)
	out = varPutU32Slice(out, uint32(len(data)))
	out = append(out, encoded...)
	return out
}

// encodeRansNx16Order0Core performs the core rANS Nx16 order-0 encoding.
// Returns freq table + states + encoded data (no flags/size prefix).
func encodeRansNx16Order0Core(data []byte) []byte {
	if len(data) == 0 {
		return nil
	}

	// Count frequencies.
	var freqs [256]uint32
	for _, b := range data {
		freqs[b]++
	}

	// Normalize to ransNx16TotFreq (4096).
	normFreqs := nx16NormalizeFreqs(freqs[:], uint32(len(data)))

	// Write frequency table.
	freqTable := nx16EncodeFreq(normFreqs[:])

	// Build cumulative frequencies for encoding.
	var cumFreq [256]uint32
	cum := uint32(0)
	for i := 0; i < 256; i++ {
		cumFreq[i] = cum
		cum += normFreqs[i]
	}

	// Encode with 4 interleaved states, 16-bit renormalization.
	nx := 4
	n := len(data)

	var state [4]uint32
	for i := 0; i < nx; i++ {
		state[i] = ransNx16L
	}

	// Encode in reverse, output 16-bit renorm words.
	var outBuf []byte

	for i := n - 1; i >= 0; i-- {
		si := i % nx
		sym := data[i]
		freq := normFreqs[sym]
		cf := cumFreq[sym]

		if freq == 0 {
			continue
		}

		// Renormalize: output 16-bit word when state is too large.
		// Output high byte first, then low byte, so after full buffer reversal
		// the pair ends up as low,high (little-endian) for the decoder.
		maxState := freq << (31 - ransNx16TFShift)
		for state[si] >= maxState {
			outBuf = append(outBuf, byte((state[si]>>8)&0xff), byte(state[si]&0xff))
			state[si] >>= 16
		}

		// Encode step.
		state[si] = ((state[si] / freq) << ransNx16TFShift) + (state[si] % freq) + cf
	}

	// Write final states.
	stateBytes := make([]byte, nx*4)
	for i := 0; i < nx; i++ {
		binary.LittleEndian.PutUint32(stateBytes[i*4:], state[i])
	}

	// Reverse output buffer.
	for i, j := 0, len(outBuf)-1; i < j; i, j = i+1, j-1 {
		outBuf[i], outBuf[j] = outBuf[j], outBuf[i]
	}

	var result []byte
	result = append(result, freqTable...)
	result = append(result, stateBytes...)
	result = append(result, outBuf...)
	return result
}

// nx16NormalizeFreqs normalizes frequencies to sum to ransNx16TotFreq (4096).
func nx16NormalizeFreqs(raw []uint32, total uint32) [256]uint32 {
	var norm [256]uint32
	if total == 0 {
		return norm
	}

	type sf struct {
		sym  int
		freq uint32
	}
	var active []sf
	for i, f := range raw {
		if f > 0 {
			active = append(active, sf{i, f})
		}
	}
	if len(active) == 0 {
		return norm
	}

	// Scale proportionally, minimum 1.
	for _, s := range active {
		f := uint64(s.freq) * uint64(ransNx16TotFreq) / uint64(total)
		if f == 0 {
			f = 1
		}
		norm[s.sym] = uint32(f)
	}

	// Fix sum.
	sort.Slice(active, func(i, j int) bool { return active[i].freq > active[j].freq })
	sum := uint32(0)
	for _, s := range active {
		sum += norm[s.sym]
	}
	if sum > ransNx16TotFreq {
		excess := sum - ransNx16TotFreq
		for i := 0; excess > 0 && i < len(active); i++ {
			s := active[i].sym
			if norm[s] > 1 {
				take := excess
				if take > norm[s]-1 {
					take = norm[s] - 1
				}
				norm[s] -= take
				excess -= take
			}
		}
	} else if sum < ransNx16TotFreq {
		norm[active[0].sym] += ransNx16TotFreq - sum
	}

	return norm
}

// nx16EncodeFreq writes the Nx16 frequency table.
// Format: alphabet (RLE) + frequencies (7-bit varint).
func nx16EncodeFreq(freqs []uint32) []byte {
	var out []byte

	// Write alphabet.
	out = nx16EncodeAlphabet(out, freqs)

	// Write frequencies.
	for j := 0; j < 256; j++ {
		if freqs[j] > 0 {
			out = varPutU32Slice(out, freqs[j])
		}
	}

	return out
}

// nx16EncodeAlphabet writes the RLE-compressed alphabet.
func nx16EncodeAlphabet(out []byte, freqs []uint32) []byte {
	first := -1
	for i := 0; i < 256; i++ {
		if freqs[i] > 0 {
			first = i
			break
		}
	}
	if first < 0 {
		out = append(out, 0)
		return out
	}

	out = append(out, byte(first))
	rle := 0

	for i := first; i < 256; {
		if rle > 0 {
			rle--
			i++
			continue
		}

		// Check for consecutive run.
		if i+1 < 256 && freqs[i+1] > 0 {
			out = append(out, byte(i+1))
			run := 0
			for j := i + 2; j < 256 && freqs[j] > 0; j++ {
				run++
			}
			out = append(out, byte(run))
			rle = run
			i++
			continue
		}

		// Find next non-zero.
		next := -1
		for j := i + 1; j < 256; j++ {
			if freqs[j] > 0 {
				next = j
				break
			}
		}
		if next < 0 {
			break
		}
		out = append(out, byte(next))
		i = next
	}

	out = append(out, 0) // terminator
	return out
}

// packNx16Encode applies the PACK transform: maps symbols to 0..nsyms-1
// and packs multiple symbols per byte when possible.
func packNx16Encode(data []byte, nsyms int) ([]byte, []byte) {
	// Build symbol map.
	var fwdMap [256]byte
	var revMap [256]byte
	idx := byte(0)
	for i := 0; i < 256; i++ {
		found := false
		for _, b := range data {
			if b == byte(i) {
				found = true
				break
			}
		}
		if found {
			fwdMap[i] = idx
			revMap[idx] = byte(i)
			idx++
		}
	}

	// Pack metadata: nsyms + symbol map.
	var meta []byte
	meta = append(meta, byte(nsyms))
	for i := 0; i < nsyms; i++ {
		meta = append(meta, revMap[i])
	}

	// Map symbols.
	mapped := make([]byte, len(data))
	for i, b := range data {
		mapped[i] = fwdMap[b]
	}

	// Pack bits based on nsyms. Uses LSB-first ordering to match the decoder.
	switch {
	case nsyms <= 1:
		// All same symbol, output empty.
		return []byte{}, meta
	case nsyms <= 2:
		// 1 bit per symbol, 8 per byte, LSB first.
		packed := make([]byte, (len(mapped)+7)/8)
		for i, b := range mapped {
			packed[i/8] |= (b & 1) << uint(i%8)
		}
		return packed, meta
	case nsyms <= 4:
		// 2 bits per symbol, 4 per byte, LSB first.
		packed := make([]byte, (len(mapped)+3)/4)
		for i, b := range mapped {
			packed[i/4] |= (b & 3) << (2 * uint(i%4))
		}
		return packed, meta
	case nsyms <= 16:
		// 4 bits per symbol, 2 per byte, LSB first.
		packed := make([]byte, (len(mapped)+1)/2)
		for i, b := range mapped {
			packed[i/2] |= (b & 0xf) << (4 * uint(i%2))
		}
		return packed, meta
	default:
		return mapped, meta
	}
}
