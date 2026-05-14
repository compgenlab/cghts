package codec

import (
	"encoding/binary"
	"sort"
)

// rANS Nx16 encoder for CRAM v3.1.
// Supports order-0/order-1 with optional PACK, RLE, STRIPE, and CAT transforms.

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

// EncodeRansNx16 encodes data using rANS Nx16 codec, trying all methods
// competitively (including STRIPE) and picking the smallest output.
func EncodeRansNx16(data []byte) []byte {
	if len(data) == 0 {
		return []byte{0, 0}
	}

	best := encodeRansNx16NoStripe(data)

	// Note: STRIPE is available via encodeRansNx16Stripe but not included in
	// competitive selection — it rarely wins and has htslib compatibility issues.

	return best
}

// encodeRansNx16NoStripe tries all methods except STRIPE to avoid recursion.
func encodeRansNx16NoStripe(data []byte) []byte {
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

	// Try multiple strategies and pick the smallest.
	best := encodeRansNx16Order0(data)

	// Try order-1.
	if len(data) > 16 {
		o1 := encodeRansNx16Order1(data)
		if len(o1) < len(best) {
			best = o1
		}
	}

	// Try PACK transform if few symbols.
	if nsyms <= 16 {
		packed, packMeta := packNx16Encode(data, nsyms)
		encoded := encodeRansNx16Order0Core(packed)
		flags := byte(ransOrderPack)

		var out []byte
		out = append(out, flags)
		out = varPutU32Slice(out, uint32(len(data)))
		out = append(out, packMeta...)
		out = varPutU32Slice(out, uint32(len(packed)))
		out = append(out, encoded...)

		if len(out) < len(best) {
			best = out
		}
	}

	// Try RLE + order-0.
	if rle := encodeRansNx16WithRLE(data, 0); rle != nil && len(rle) < len(best) {
		best = rle
	}

	// Try RLE + order-1 for larger inputs.
	if len(data) > 16 {
		if rle := encodeRansNx16WithRLE(data, 1); rle != nil && len(rle) < len(best) {
			best = rle
		}
	}

	// Try CAT (raw passthrough) as fallback for incompressible data.
	cat := encodeRansNx16Cat(data)
	if len(cat) < len(best) {
		best = cat
	}

	return best
}

// encodeRansNx16Order1 encodes data with order-1 rANS Nx16.
func encodeRansNx16Order1(data []byte) []byte {
	flags := byte(ransOrderMask) // order-1
	encoded := encodeRansNx16Order1Core(data)

	var out []byte
	out = append(out, flags)
	out = varPutU32Slice(out, uint32(len(data)))
	out = append(out, encoded...)
	return out
}

// encodeRansNx16Order1Core performs the core rANS Nx16 order-1 encoding.
// Returns freq tables + states + encoded data (no flags/size prefix).
func encodeRansNx16Order1Core(data []byte) []byte {
	if len(data) == 0 {
		return nil
	}

	nx := 4
	n := len(data)
	tfShift := uint(ransNx16TFShift)

	// Count per-context frequencies using NX-way split contexts (matching decoder).
	var freqs [256][256]uint32
	var totals [256]uint32

	// The decoder uses NX-way split: state z processes positions z*iszN .. (z+1)*iszN-1,
	// with the last state handling the remainder.
	iszN := n / nx
	for z := 0; z < nx; z++ {
		start := z * iszN
		end := start + iszN
		if z == nx-1 {
			end = n
		}
		prev := byte(0)
		for i := start; i < end; i++ {
			freqs[prev][data[i]]++
			totals[prev]++
			prev = data[i]
		}
	}

	// Build the alphabet (F0): includes all byte values that appear as EITHER
	// contexts or symbols, since the decoder uses F0 for both.
	var alpha [256]uint32
	for i := 0; i < 256; i++ {
		if totals[i] > 0 {
			alpha[i] = 1
		}
	}
	for ctx := 0; ctx < 256; ctx++ {
		if totals[ctx] == 0 {
			continue
		}
		for j := 0; j < 256; j++ {
			if freqs[ctx][j] > 0 {
				alpha[j] = 1
			}
		}
	}

	// Ensure every symbol in the alphabet also has a context entry (even if empty).
	// The decoder reads a freq table for every ctx where alpha[ctx] > 0.
	// For contexts that have no real data, we need a valid (all-zero) freq table
	// that sums to totFreq — give them a uniform distribution.
	for i := 0; i < 256; i++ {
		if alpha[i] > 0 && totals[i] == 0 {
			// This symbol appears in data but never as a context.
			// Give it a trivial frequency table (uniform over all alpha symbols).
			totals[i] = 0 // will be handled below
		}
	}

	// Normalize each context's frequencies.
	var normFreqs [256][256]uint32
	var cumFreqs [256][256]uint32
	for ctx := 0; ctx < 256; ctx++ {
		if alpha[ctx] == 0 {
			continue
		}
		if totals[ctx] > 0 {
			// Real context: normalize using only symbols in the alphabet.
			var filteredFreqs [256]uint32
			var filteredTotal uint32
			for j := 0; j < 256; j++ {
				if alpha[j] > 0 && freqs[ctx][j] > 0 {
					filteredFreqs[j] = freqs[ctx][j]
					filteredTotal += freqs[ctx][j]
				}
			}
			normFreqs[ctx] = nx16NormalizeFreqs(filteredFreqs[:], filteredTotal)
		} else {
			// Synthetic context: uniform distribution over all alpha symbols.
			var synthFreqs [256]uint32
			count := uint32(0)
			for j := 0; j < 256; j++ {
				if alpha[j] > 0 {
					synthFreqs[j] = 1
					count++
				}
			}
			normFreqs[ctx] = nx16NormalizeFreqs(synthFreqs[:], count)
		}
		cum := uint32(0)
		for j := 0; j < 256; j++ {
			cumFreqs[ctx][j] = cum
			cum += normFreqs[ctx][j]
		}
	}

	// Build the frequency table.
	freqTable := nx16EncodeFreqOrder1(normFreqs[:], alpha[:])

	// First byte: upper nibble = TF_SHIFT, bit 0 = compressed freq table flag.
	firstByte := byte(tfShift<<4) | 0 // uncompressed for now

	// Encode with NX-way split, in reverse.
	var state [4]uint32
	for i := 0; i < nx; i++ {
		state[i] = ransNx16L
	}

	var outBuf []byte

	// Build context arrays for reverse traversal per split.
	type posCtx struct {
		sym byte
		ctx byte
	}
	splits := make([][]posCtx, nx)
	for z := 0; z < nx; z++ {
		start := z * iszN
		end := start + iszN
		if z == nx-1 {
			end = n
		}
		s := make([]posCtx, end-start)
		prev := byte(0)
		for i := start; i < end; i++ {
			s[i-start] = posCtx{sym: data[i], ctx: prev}
			prev = data[i]
		}
		splits[z] = s
	}

	// Encode in reverse. The decoder processes: outer loop over positions,
	// inner loop over z from 0..nx-1. So the encoder must reverse that.
	// The last state handles the remainder first (in reverse).

	// First encode the remainder (last state, positions beyond iszN).
	lastSplit := splits[nx-1]
	for i := len(lastSplit) - 1; i >= iszN; i-- {
		sym := lastSplit[i].sym
		ctx := lastSplit[i].ctx
		freq := normFreqs[ctx][sym]
		cf := cumFreqs[ctx][sym]
		if freq == 0 {
			continue
		}
		maxState := freq << (31 - tfShift)
		for state[nx-1] >= maxState {
			outBuf = append(outBuf, byte((state[nx-1]>>8)&0xff), byte(state[nx-1]&0xff))
			state[nx-1] >>= 16
		}
		state[nx-1] = ((state[nx-1] / freq) << tfShift) + (state[nx-1] % freq) + cf
	}

	// Then encode the main body: positions iszN-1 down to 0, inner z from nx-1 to 0.
	for i := iszN - 1; i >= 0; i-- {
		for z := nx - 1; z >= 0; z-- {
			if i >= len(splits[z]) {
				continue
			}
			sym := splits[z][i].sym
			ctx := splits[z][i].ctx
			freq := normFreqs[ctx][sym]
			cf := cumFreqs[ctx][sym]
			if freq == 0 {
				continue
			}
			maxState := freq << (31 - tfShift)
			for state[z] >= maxState {
				outBuf = append(outBuf, byte((state[z]>>8)&0xff), byte(state[z]&0xff))
				state[z] >>= 16
			}
			state[z] = ((state[z] / freq) << tfShift) + (state[z] % freq) + cf
		}
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
	result = append(result, firstByte)
	result = append(result, freqTable...)
	result = append(result, stateBytes...)
	result = append(result, outBuf...)
	return result
}

// nx16EncodeFreqOrder1 writes order-1 frequency tables in the freqD format.
func nx16EncodeFreqOrder1(normFreqs [][256]uint32, totals []uint32) []byte {
	// Write the alphabet (which contexts are present).
	var alphaFreqs [256]uint32
	for ctx := 0; ctx < 256; ctx++ {
		if totals[ctx] > 0 {
			alphaFreqs[ctx] = 1
		}
	}
	out := nx16EncodeAlphabet(nil, alphaFreqs[:])

	// For each context, write freqD format.
	for ctx := 0; ctx < 256; ctx++ {
		if totals[ctx] == 0 {
			continue
		}
		out = nx16EncodeFreqD(out, normFreqs[ctx][:], alphaFreqs[:])
	}
	return out
}

// nx16EncodeFreqD writes a per-context frequency table in delta/zero-run format.
func nx16EncodeFreqD(out []byte, freqs []uint32, F0 []uint32) []byte {
	dz := 0
	for j := 0; j < 256; j++ {
		if F0[j] == 0 {
			continue
		}
		if dz > 0 {
			dz--
			continue
		}
		f := freqs[j]
		out = varPutU32Slice(out, f)
		if f == 0 {
			// Count consecutive zeros.
			run := 0
			for k := j + 1; k < 256; k++ {
				if F0[k] == 0 {
					continue
				}
				if freqs[k] != 0 {
					break
				}
				run++
			}
			out = append(out, byte(run))
			dz = run
		}
	}
	return out
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

// rleEncodeNx16 applies the RLE transform to data.
// Returns literals (for core rANS encoding), meta (nSyms + syms + run varints), and ok.
// If no runs are found, returns ok=false.
func rleEncodeNx16(data []byte) (literals []byte, meta []byte, ok bool) {
	// Find symbols with consecutive runs (2+ same byte in a row).
	var hasRun [256]bool
	for i := 1; i < len(data); i++ {
		if data[i] == data[i-1] {
			hasRun[data[i]] = true
		}
	}

	var rleSyms []byte
	for i := 0; i < 256; i++ {
		if hasRun[i] {
			rleSyms = append(rleSyms, byte(i))
		}
	}

	if len(rleSyms) == 0 {
		return nil, nil, false
	}

	var isRLE [256]bool
	for _, s := range rleSyms {
		isRLE[s] = true
	}

	// Build literals and run lengths.
	var runs []byte
	i := 0
	for i < len(data) {
		b := data[i]
		if isRLE[b] {
			runLen := 1
			for i+runLen < len(data) && data[i+runLen] == b {
				runLen++
			}
			literals = append(literals, b)
			runs = varPutU32Slice(runs, uint32(runLen-1))
			i += runLen
		} else {
			literals = append(literals, b)
			i++
		}
	}

	// Build meta: [nSyms] [syms...] [runs...]
	nSymByte := byte(len(rleSyms))
	if len(rleSyms) == 256 {
		nSymByte = 0
	}
	meta = append(meta, nSymByte)
	meta = append(meta, rleSyms...)
	meta = append(meta, runs...)

	return literals, meta, true
}

// encodeRansNx16WithRLE encodes data with RLE transform + core rANS.
// order is 0 or 1. Returns nil if RLE is not applicable.
func encodeRansNx16WithRLE(data []byte, order int) []byte {
	literals, rleMeta, ok := rleEncodeNx16(data)
	if !ok {
		return nil
	}

	// Compress literals with core rANS.
	var encoded []byte
	if order == 0 {
		encoded = encodeRansNx16Order0Core(literals)
	} else {
		encoded = encodeRansNx16Order1Core(literals)
	}

	flags := byte(ransOrderRLE)
	if order == 1 {
		flags |= ransOrderMask
	}

	// Use uncompressed meta (odd uMetaSize flag).
	uMetaSize := uint32(len(rleMeta))*2 + 1

	var out []byte
	out = append(out, flags)
	out = varPutU32Slice(out, uint32(len(data)))    // original size
	out = varPutU32Slice(out, uMetaSize)             // meta size (odd = uncompressed)
	out = varPutU32Slice(out, uint32(len(literals))) // rleLen (core decode output size)
	out = append(out, rleMeta...)
	out = append(out, encoded...)

	return out
}

// encodeRansNx16Cat encodes data as raw passthrough (CAT).
func encodeRansNx16Cat(data []byte) []byte {
	var out []byte
	out = append(out, byte(ransOrderCat))
	out = varPutU32Slice(out, uint32(len(data)))
	out = append(out, data...)
	return out
}

// encodeRansNx16Stripe encodes data with the STRIPE transform.
// Splits data into N byte-interleaved sub-streams, compresses each independently.
func encodeRansNx16Stripe(data []byte, N int) []byte {
	if N < 2 || len(data) < N {
		return nil
	}

	ulen := len(data)

	// Compute per-stripe sizes.
	ulenN := make([]int, N)
	for i := 0; i < N; i++ {
		ulenN[i] = ulen / N
		if ulen%N > i {
			ulenN[i]++
		}
	}

	// De-interleave into N sub-streams.
	substreams := make([][]byte, N)
	for i := 0; i < N; i++ {
		substreams[i] = make([]byte, ulenN[i])
		for j := 0; j < ulenN[i]; j++ {
			substreams[i][j] = data[j*N+i]
		}
	}

	// Compress each sub-stream (without STRIPE to avoid recursion).
	compressed := make([][]byte, N)
	for i := 0; i < N; i++ {
		compressed[i] = encodeRansNx16NoStripe(substreams[i])
	}

	// Build output.
	var out []byte
	out = append(out, byte(ransOrderStripe))
	out = varPutU32Slice(out, uint32(ulen))
	out = append(out, byte(N))
	for i := 0; i < N; i++ {
		out = varPutU32Slice(out, uint32(len(compressed[i])))
	}
	for i := 0; i < N; i++ {
		out = append(out, compressed[i]...)
	}

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
