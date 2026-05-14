package codec

// Range coder encoder (reverse of the decoder in range_coder.go).
// Based on Eugene Shelwien's public domain range coder, as used in htscodecs.
//
// Uses uint64 for low to detect carry. After shifting, low is truncated
// to uint32 (matching the htscodecs RC_ShiftLow implementation).

type rangeEncoder struct {
	low     uint64
	rng     uint32
	buf     []byte
	cache   byte // pending output byte (starts at 0)
	ffCount int  // number of pending 0xFF bytes
}

func newRangeEncoder() *rangeEncoder {
	return &rangeEncoder{
		rng:   0xFFFFFFFF,
		cache: 0,
	}
}

func (re *rangeEncoder) encode(cumFreq, freq, totFreq uint32) {
	re.rng /= totFreq
	re.low += uint64(cumFreq) * uint64(re.rng)
	re.rng *= freq
	for re.rng < rcTop {
		re.shiftLow()
		re.rng <<= 8
	}
}

func (re *rangeEncoder) shiftLow() {
	low32 := uint32(re.low)
	carry := byte(re.low >> 32)

	if low32 < 0xFF000000 || carry != 0 {
		re.buf = append(re.buf, re.cache+carry)
		for re.ffCount > 0 {
			if carry != 0 {
				re.buf = append(re.buf, 0x00)
			} else {
				re.buf = append(re.buf, 0xFF)
			}
			re.ffCount--
		}
		re.cache = byte(low32 >> 24)
	} else {
		re.ffCount++
	}
	// Truncate to uint32, then shift left 8.
	re.low = uint64(low32 << 8)
}

func (re *rangeEncoder) finish() []byte {
	for i := 0; i < 5; i++ {
		re.shiftLow()
	}
	// Flush the final cached byte and any pending FFs.
	re.buf = append(re.buf, re.cache)
	for re.ffCount > 0 {
		re.buf = append(re.buf, 0xFF)
		re.ffCount--
	}
	return re.buf
}
