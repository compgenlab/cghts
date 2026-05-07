package codec

import (
	"encoding/binary"
	"sort"
)

// rANS 4x8 encoder for CRAM v3.0+.
// Implements order-0 and order-1 encoding matching the htslib format.

// ransEncSymbol holds encoding info for one symbol.
type ransEncSymbol struct {
	freq    uint32
	cumFreq uint32
}

// encodeRans4x8 compresses data using rANS 4x8 codec.
// order is 0 or 1. Returns the compressed data including the order byte prefix.
func EncodeRans4x8(data []byte, order int) []byte {
	if order == 1 {
		return encodeRansOrder1(data)
	}
	return encodeRansOrder0(data)
}

func encodeRansOrder0(data []byte) []byte {
	if len(data) == 0 {
		return []byte{0, 0, 0, 0, 0, 0, 0, 0, 0}
	}

	// Count frequencies.
	var freqs [256]uint32
	for _, b := range data {
		freqs[b]++
	}

	// Normalize frequencies to sum to ransTotFreq.
	syms := normalizeFreqs(freqs[:], uint32(len(data)))

	// Encode frequency table.
	freqTable := writeFreqTableOrder0(syms[:])

	// Encode data using rANS with 4 interleaved states.
	encoded := ransEncode4x8(data, syms[:])

	// Build output: order(1) + compSize(4) + uncompSize(4) + freqTable + encoded
	compSize := uint32(len(freqTable) + len(encoded))
	uncompSize := uint32(len(data))

	var out []byte
	out = append(out, 0) // order = 0
	out = binary.LittleEndian.AppendUint32(out, compSize)
	out = binary.LittleEndian.AppendUint32(out, uncompSize)
	out = append(out, freqTable...)
	out = append(out, encoded...)
	return out
}

func encodeRansOrder1(data []byte) []byte {
	if len(data) == 0 {
		return []byte{1, 0, 0, 0, 0, 0, 0, 0, 0}
	}

	// Count context-dependent frequencies using per-stream contexts
	// (matching the decoder's 4-way interleaving).
	var freqs [256][256]uint32
	var totals [256]uint32
	var prevSym [4]byte
	for i, b := range data {
		si := i & 3
		ctx := prevSym[si]
		freqs[ctx][b]++
		totals[ctx]++
		prevSym[si] = b
	}

	// Normalize each context's frequencies.
	var syms [256][256]ransEncSymbol
	for ctx := 0; ctx < 256; ctx++ {
		if totals[ctx] == 0 {
			continue
		}
		syms[ctx] = normalizeFreqs(freqs[ctx][:], totals[ctx])
	}

	// Encode frequency tables.
	freqTable := writeFreqTableOrder1(syms[:], totals[:])

	// Encode data.
	encoded := ransEncode4x8Order1(data, syms[:])

	compSize := uint32(len(freqTable) + len(encoded))
	uncompSize := uint32(len(data))

	var out []byte
	out = append(out, 1) // order = 1
	out = binary.LittleEndian.AppendUint32(out, compSize)
	out = binary.LittleEndian.AppendUint32(out, uncompSize)
	out = append(out, freqTable...)
	out = append(out, encoded...)
	return out
}

// normalizeFreqs normalizes raw frequencies to sum to ransTotFreq (4096).
func normalizeFreqs(raw []uint32, total uint32) [256]ransEncSymbol {
	var syms [256]ransEncSymbol

	if total == 0 {
		return syms
	}

	// Find non-zero symbols.
	type symFreq struct {
		sym  int
		freq uint32
	}
	var active []symFreq
	for i, f := range raw {
		if f > 0 {
			active = append(active, symFreq{i, f})
		}
	}

	if len(active) == 0 {
		return syms
	}

	// Scale frequencies to sum to ransTotFreq.
	// First pass: proportional scaling, ensuring minimum of 1.
	var scaled [256]uint32
	remaining := uint32(ransTotFreq)
	for _, sf := range active {
		f := uint64(sf.freq) * uint64(ransTotFreq) / uint64(total)
		if f == 0 {
			f = 1
		}
		scaled[sf.sym] = uint32(f)
		remaining -= uint32(f)
	}

	// Adjust to hit exact total. Sort by descending frequency for stability.
	sort.Slice(active, func(i, j int) bool {
		return active[i].freq > active[j].freq
	})

	// Fix surplus/deficit by adjusting the most frequent symbols.
	sum := uint32(0)
	for _, sf := range active {
		sum += scaled[sf.sym]
	}

	if sum > ransTotFreq {
		// Remove excess from largest frequencies.
		excess := sum - ransTotFreq
		for i := 0; excess > 0 && i < len(active); i++ {
			s := active[i].sym
			if scaled[s] > 1 {
				take := excess
				if take > scaled[s]-1 {
					take = scaled[s] - 1
				}
				scaled[s] -= take
				excess -= take
			}
		}
	} else if sum < ransTotFreq {
		// Add deficit to largest frequencies.
		deficit := ransTotFreq - sum
		for i := 0; deficit > 0 && i < len(active); i++ {
			s := active[i].sym
			scaled[s] += deficit
			deficit = 0
		}
	}

	// Build cumulative frequencies.
	cumFreq := uint32(0)
	for i := 0; i < 256; i++ {
		syms[i] = ransEncSymbol{freq: scaled[i], cumFreq: cumFreq}
		cumFreq += scaled[i]
	}

	return syms
}

// writeFreqTableOrder0 encodes the frequency table in htslib's format.
func writeFreqTableOrder0(syms []ransEncSymbol) []byte {
	var out []byte

	// Find first non-zero symbol.
	first := -1
	for i := 0; i < 256; i++ {
		if syms[i].freq > 0 {
			first = i
			break
		}
	}
	if first < 0 {
		out = append(out, 0) // empty table terminator
		return out
	}

	out = append(out, byte(first))

	rle := 0
	for i := first; i < 256; {
		// Write frequency (1 or 2 bytes).
		f := syms[i].freq
		if f >= 128 {
			out = append(out, byte((f>>8)+128), byte(f&0xff))
		} else {
			out = append(out, byte(f))
		}

		// Check for consecutive run.
		if rle > 0 {
			rle--
			i++
			continue
		}

		// Look for next non-zero symbol.
		if i+1 < 256 && syms[i+1].freq > 0 {
			// Next symbol is consecutive — check for RLE.
			out = append(out, byte(i+1))
			// Count run length.
			run := 0
			for j := i + 2; j < 256 && syms[j].freq > 0; j++ {
				run++
			}
			out = append(out, byte(run))
			rle = run
			i++
			continue
		}

		// Find next non-zero symbol.
		next := -1
		for j := i + 1; j < 256; j++ {
			if syms[j].freq > 0 {
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

// writeFreqTableOrder1 encodes order-1 frequency tables.
func writeFreqTableOrder1(syms [][256]ransEncSymbol, totals []uint32) []byte {
	var out []byte

	// Find first context with non-zero total.
	first := -1
	for i := 0; i < 256; i++ {
		if totals[i] > 0 {
			first = i
			break
		}
	}
	if first < 0 {
		out = append(out, 0)
		return out
	}

	out = append(out, byte(first))

	rleCtx := 0
	for ctx := first; ctx < 256; {
		if totals[ctx] == 0 {
			ctx++
			continue
		}

		// Write inner frequency table for this context.
		inner := writeFreqTableOrder0(syms[ctx][:])
		out = append(out, inner...)

		if rleCtx > 0 {
			rleCtx--
			ctx++
			continue
		}

		// Look for next context with data.
		if ctx+1 < 256 && totals[ctx+1] > 0 {
			out = append(out, byte(ctx+1))
			run := 0
			for j := ctx + 2; j < 256 && totals[j] > 0; j++ {
				run++
			}
			out = append(out, byte(run))
			rleCtx = run
			ctx++
			continue
		}

		next := -1
		for j := ctx + 1; j < 256; j++ {
			if totals[j] > 0 {
				next = j
				break
			}
		}
		if next < 0 {
			break
		}
		out = append(out, byte(next))
		ctx = next
	}

	out = append(out, 0) // terminator
	return out
}

// ransEncode4x8 encodes data using 4 interleaved rANS states (order 0).
func ransEncode4x8(data []byte, syms []ransEncSymbol) []byte {
	n := len(data)
	if n == 0 {
		return nil
	}

	// Encode in reverse.
	var state [4]uint32
	for i := 0; i < 4; i++ {
		state[i] = ransL
	}

	// Output buffer (built in reverse).
	outBuf := make([]byte, 0, n)

	// Process from the end. The interleaving pattern is: byte i is processed by state (i & 3).
	// We process in reverse: last byte first, working backward.
	for i := n - 1; i >= 0; i-- {
		si := i & 3
		sym := data[i]
		freq := syms[sym].freq
		cumFreq := syms[sym].cumFreq

		if freq == 0 {
			continue
		}

		// Renormalize: while state >= freq * (ransL >> ransTFShift), output a byte.
		maxState := freq << (31 - ransTFShift)
		for state[si] >= maxState {
			outBuf = append(outBuf, byte(state[si]&0xff))
			state[si] >>= 8
		}

		// Encode step: state = (state / freq) * ransTotFreq + (state % freq) + cumFreq
		state[si] = ((state[si] / freq) << ransTFShift) + (state[si] % freq) + cumFreq
	}

	// Write final states.
	var stateBytes [16]byte
	for i := 0; i < 4; i++ {
		binary.LittleEndian.PutUint32(stateBytes[i*4:], state[i])
	}

	// Reverse the output buffer.
	for i, j := 0, len(outBuf)-1; i < j; i, j = i+1, j-1 {
		outBuf[i], outBuf[j] = outBuf[j], outBuf[i]
	}

	// Final: states + reversed data
	var result []byte
	result = append(result, stateBytes[:]...)
	result = append(result, outBuf...)
	return result
}

// ransEncode4x8Order1 encodes data using 4 interleaved rANS states (order 1).
func ransEncode4x8Order1(data []byte, syms [][256]ransEncSymbol) []byte {
	n := len(data)
	if n == 0 {
		return nil
	}

	var state [4]uint32
	for i := 0; i < 4; i++ {
		state[i] = ransL
	}

	outBuf := make([]byte, 0, n)

	// Track context per stream (previous symbol for each of the 4 interleaved streams).
	var lastSym [4]byte

	// Build context array for reverse traversal.
	// Forward pass to determine contexts.
	contexts := make([]byte, n)
	var fwdCtx [4]byte
	for i := 0; i < n; i++ {
		si := i & 3
		contexts[i] = fwdCtx[si]
		fwdCtx[si] = data[i]
	}
	_ = lastSym

	for i := n - 1; i >= 0; i-- {
		si := i & 3
		sym := data[i]
		ctx := contexts[i]
		freq := syms[ctx][sym].freq
		cumFreq := syms[ctx][sym].cumFreq

		if freq == 0 {
			continue
		}

		maxState := freq << (31 - ransTFShift)
		for state[si] >= maxState {
			outBuf = append(outBuf, byte(state[si]&0xff))
			state[si] >>= 8
		}

		state[si] = ((state[si] / freq) << ransTFShift) + (state[si] % freq) + cumFreq
	}

	var stateBytes [16]byte
	for i := 0; i < 4; i++ {
		binary.LittleEndian.PutUint32(stateBytes[i*4:], state[i])
	}

	for i, j := 0, len(outBuf)-1; i < j; i, j = i+1, j-1 {
		outBuf[i], outBuf[j] = outBuf[j], outBuf[i]
	}

	var result []byte
	result = append(result, stateBytes[:]...)
	result = append(result, outBuf...)
	return result
}
