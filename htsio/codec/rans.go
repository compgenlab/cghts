package codec

import (
	"encoding/binary"
	"fmt"
)

// rANS 4x8 decoder.
// This implements the order-0 and order-1 rANS codec used in CRAM v3.0+.
// Based on the htslib reference implementation.

const (
	ransTFShift = 12                  // frequency table precision (bottom bits of state)
	ransTotFreq = 1 << ransTFShift    // 4096: frequencies sum to this
	ransL       = 1 << 23            // 8388608: renormalization threshold (RANS_BYTE_L)
)

// ransDecSymbol holds precomputed decode info for one symbol.
type ransDecSymbol struct {
	cumFreq uint32 // cumulative frequency
	freq    uint32 // symbol frequency
}

func DecodeRans4x8(data []byte) ([]byte, error) {
	if len(data) < 1 {
		return nil, fmt.Errorf("rans: empty data")
	}

	order := data[0]
	data = data[1:]

	switch order {
	case 0:
		return decodeRansOrder0(data)
	case 1:
		return decodeRansOrder1(data)
	default:
		return nil, fmt.Errorf("rans: unsupported order: %d", order)
	}
}

func decodeRansOrder0(data []byte) ([]byte, error) {
	if len(data) < 8 {
		return nil, fmt.Errorf("rans0: data too short")
	}

	// Read compressed and uncompressed sizes.
	compSize := binary.LittleEndian.Uint32(data[0:4])
	uncompSize := binary.LittleEndian.Uint32(data[4:8])
	data = data[8:]

	if int(compSize) > len(data) {
		return nil, fmt.Errorf("rans0: compressed size %d exceeds data length %d", compSize, len(data))
	}
	data = data[:compSize]

	// Read frequency table.
	syms, cumTotal, pos, err := readFreqTableOrder0(data)
	if err != nil {
		return nil, err
	}
	data = data[pos:]

	// Build decode LUT (symbol, freq, base at each cumulative position).
	var lutSym [ransTotFreq]byte
	var lutFreq [ransTotFreq]uint32
	var lutBase [ransTotFreq]uint32
	for s := 0; s < 256; s++ {
		ds := syms[s]
		for j := uint32(0); j < ds.freq && ds.cumFreq+j < ransTotFreq; j++ {
			lutSym[ds.cumFreq+j] = byte(s)
			lutFreq[ds.cumFreq+j] = ds.freq
			lutBase[ds.cumFreq+j] = j
		}
	}
	// Handle the case where frequencies sum to ransTotFreq-1 (4095).
	if cumTotal == ransTotFreq-1 {
		lutSym[ransTotFreq-1] = lutSym[ransTotFreq-2]
		lutFreq[ransTotFreq-1] = lutFreq[ransTotFreq-2]
		lutBase[ransTotFreq-1] = lutBase[ransTotFreq-2] + 1
	}

	// Initialize 4 rANS states.
	if len(data) < 16 {
		return nil, fmt.Errorf("rans0: not enough data for initial states")
	}
	var state [4]uint32
	for i := 0; i < 4; i++ {
		state[i] = binary.LittleEndian.Uint32(data[i*4:])
	}
	data = data[16:]

	// Decode in batches of 4, matching htslib's interleaving.
	// All 4 states are advanced first, then renormalized in order,
	// because they share the same compressed byte stream.
	out := make([]byte, uncompSize)
	dataPos := 0
	outEnd := int(uncompSize) &^ 3 // round down to multiple of 4

	renorm := func(si int) {
		if state[si] >= ransL || dataPos >= len(data) {
			return
		}
		state[si] = (state[si] << 8) | uint32(data[dataPos])
		dataPos++
		if state[si] < ransL && dataPos < len(data) {
			state[si] = (state[si] << 8) | uint32(data[dataPos])
			dataPos++
		}
	}

	for i := 0; i < outEnd; i += 4 {
		// Decode and advance all 4 states.
		var m [4]uint32
		for j := 0; j < 4; j++ {
			m[j] = state[j] & (ransTotFreq - 1)
			state[j] = lutFreq[m[j]]*(state[j]>>ransTFShift) + lutBase[m[j]]
		}

		// Renormalize all 4 states (reads from shared byte stream).
		renorm(0)
		renorm(1)
		renorm(2)
		renorm(3)

		// Write output symbols.
		out[i+0] = lutSym[m[0]]
		out[i+1] = lutSym[m[1]]
		out[i+2] = lutSym[m[2]]
		out[i+3] = lutSym[m[3]]
	}

	// Handle remaining 0-3 bytes.
	for i := outEnd; i < int(uncompSize); i++ {
		si := i & 3
		m := state[si] & (ransTotFreq - 1)
		out[i] = lutSym[m]
		state[si] = lutFreq[m]*(state[si]>>ransTFShift) + lutBase[m]
		renorm(si)
	}

	return out, nil
}

// readFreqTableOrder0 reads the rANS order-0 frequency table.
// Format follows htslib's rANS_static.c encoding with run-length compression
// for consecutive symbols.
func readFreqTableOrder0(data []byte) ([256]ransDecSymbol, uint32, int, error) {
	var syms [256]ransDecSymbol
	pos := 0

	if pos >= len(data) {
		return syms, 0, pos, fmt.Errorf("rans0: freq table truncated")
	}

	j := int(data[pos])
	pos++
	cumTotal := uint32(0)
	rle := 0

	for {
		if pos >= len(data) {
			return syms, 0, pos, fmt.Errorf("rans0: freq table truncated at sym %d", j)
		}

		// Read frequency (1 or 2 bytes).
		freq := uint32(data[pos])
		pos++
		if freq >= 128 {
			if pos >= len(data) {
				return syms, 0, pos, fmt.Errorf("rans0: freq table truncated")
			}
			freq = ((freq - 128) << 8) | uint32(data[pos])
			pos++
		}

		syms[j] = ransDecSymbol{cumFreq: cumTotal, freq: freq}
		cumTotal += freq
		// Frequencies are read from untrusted input; their running total must
		// not exceed ransTotFreq or buildCumLUT would index past the LUT.
		if cumTotal > ransTotFreq {
			return syms, 0, pos, fmt.Errorf("rans0: cumulative frequency %d exceeds %d", cumTotal, ransTotFreq)
		}

		// Advance to next symbol.
		if rle == 0 && j < 255 && pos < len(data) && data[pos] == byte(j+1) {
			// Consecutive symbol — read run length.
			j = int(data[pos])
			pos++
			if pos >= len(data) {
				return syms, 0, pos, fmt.Errorf("rans0: freq table truncated at rle")
			}
			rle = int(data[pos])
			pos++
		} else if rle > 0 {
			rle--
			j++
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

	return syms, cumTotal, pos, nil
}

func buildCumLUT(syms [256]ransDecSymbol) [ransTotFreq]byte {
	var lut [ransTotFreq]byte
	for s := 0; s < 256; s++ {
		ds := syms[s]
		for j := uint32(0); j < ds.freq; j++ {
			// Defensive bound: the frequency-table readers reject totals beyond
			// ransTotFreq, but never index past the LUT even if a caller does
			// not validate first.
			if ds.cumFreq+j >= ransTotFreq {
				break
			}
			lut[ds.cumFreq+j] = byte(s)
		}
	}
	return lut
}

func decodeRansOrder1(data []byte) ([]byte, error) {
	if len(data) < 8 {
		return nil, fmt.Errorf("rans1: data too short")
	}

	compSize := binary.LittleEndian.Uint32(data[0:4])
	uncompSize := binary.LittleEndian.Uint32(data[4:8])
	data = data[8:]

	if int(compSize) > len(data) {
		return nil, fmt.Errorf("rans1: compressed size %d exceeds data length %d", compSize, len(data))
	}
	data = data[:compSize]

	// Read order-1 frequency tables: for each context symbol, a frequency table.
	var syms [256][256]ransDecSymbol
	var luts [256][ransTotFreq]byte
	pos := 0

	if pos >= len(data) {
		return nil, fmt.Errorf("rans1: frequency table truncated")
	}

	ctxI := int(data[pos])
	pos++
	rleCtx := 0

	for {
		// Read frequency table for this context using order-0 inner format.
		if pos >= len(data) {
			break
		}

		j := int(data[pos])
		pos++
		cumTotal := uint32(0)
		rle := 0

		for {
			if pos >= len(data) {
				return nil, fmt.Errorf("rans1: freq table truncated at ctx=%d sym=%d", ctxI, j)
			}

			freq := uint32(data[pos])
			pos++
			if freq >= 128 {
				if pos >= len(data) {
					return nil, fmt.Errorf("rans1: freq table truncated")
				}
				freq = ((freq - 128) << 8) | uint32(data[pos])
				pos++
			}

			syms[ctxI][j] = ransDecSymbol{cumFreq: cumTotal, freq: freq}
			cumTotal += freq
			// Bound the running total (see order-0 path above): an unchecked
			// total lets buildCumLUT index past the LUT on malformed input.
			if cumTotal > ransTotFreq {
				return nil, fmt.Errorf("rans1: cumulative frequency %d exceeds %d at ctx=%d", cumTotal, ransTotFreq, ctxI)
			}

			if rle == 0 && j < 255 && pos < len(data) && data[pos] == byte(j+1) {
				j = int(data[pos])
				pos++
				if pos >= len(data) {
					break
				}
				rle = int(data[pos])
				pos++
			} else if rle > 0 {
				rle--
				j++
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

		// Build LUT for this context.
		luts[ctxI] = buildCumLUT(syms[ctxI])

		// Advance to next context.
		if rleCtx == 0 && ctxI < 255 && pos < len(data) && data[pos] == byte(ctxI+1) {
			ctxI = int(data[pos])
			pos++
			if pos >= len(data) {
				break
			}
			rleCtx = int(data[pos])
			pos++
		} else if rleCtx > 0 {
			rleCtx--
			ctxI++
		} else {
			if pos >= len(data) {
				break
			}
			ctxI = int(data[pos])
			pos++
			if ctxI == 0 {
				break
			}
		}
	}

	data = data[pos:]

	// Initialize 4 rANS states.
	if len(data) < 16 {
		return nil, fmt.Errorf("rans1: not enough data for initial states")
	}
	var state [4]uint32
	for i := 0; i < 4; i++ {
		state[i] = binary.LittleEndian.Uint32(data[i*4:])
	}
	data = data[16:]

	// Decode interleaved.
	out := make([]byte, uncompSize)
	dataPos := 0
	var lastSym [4]byte

	for i := 0; i < int(uncompSize); i++ {
		si := i & 3
		context := lastSym[si]
		freq := state[si] & (ransTotFreq - 1)
		sym := luts[context][freq]
		out[i] = sym
		lastSym[si] = sym

		ds := syms[context][sym]
		state[si] = ds.freq*(state[si]>>ransTFShift) + freq - ds.cumFreq

		// Renormalize.
		for state[si] < ransL && dataPos < len(data) {
			state[si] = (state[si] << 8) | uint32(data[dataPos])
			dataPos++
		}
	}

	return out, nil
}
