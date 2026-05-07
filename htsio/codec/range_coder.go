package codec

// Range coder for adaptive arithmetic coding (used by arith_dynamic and fqzcomp).
// Based on Eugene Shelwien's public domain range coder, as used in htscodecs.

const rcTop = 1 << 24

type rangeCoder struct {
	low, code, rng uint32
	buf            []byte
	pos            int
	err            bool
}

func newRangeDecoder(data []byte) *rangeCoder {
	rc := &rangeCoder{
		rng: 0xFFFFFFFF,
		buf: data,
	}
	if len(data) < 5 {
		rc.err = true
		return rc
	}
	for i := 0; i < 5; i++ {
		rc.code = (rc.code << 8) | uint32(data[i])
	}
	rc.pos = 5
	return rc
}

func (rc *rangeCoder) getFreq(totFreq uint32) uint32 {
	if totFreq == 0 || rc.rng < totFreq {
		return 0
	}
	rc.rng /= totFreq
	return rc.code / rc.rng
}

func (rc *rangeCoder) decode(cumFreq, freq, totFreq uint32) {
	rc.code -= cumFreq * rc.rng
	rc.rng *= freq
	for rc.rng < rcTop {
		if rc.pos >= len(rc.buf) {
			rc.err = true
			return
		}
		rc.code = (rc.code << 8) | uint32(rc.buf[rc.pos])
		rc.pos++
		rc.rng <<= 8
	}
}

func (rc *rangeCoder) finish() bool {
	return !rc.err
}
