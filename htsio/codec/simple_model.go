package codec

// Adaptive simple frequency model for use with the range coder.
// Based on the htscodecs c_simple_model.h reference implementation.
//
// Maintains a list of symbols approximately sorted by frequency.
// When a symbol is decoded, its frequency is bumped and it may
// swap with its neighbor to maintain approximate sort order.

const (
	modelMaxFreq = (1 << 16) - 17
	modelStep    = 16
)

type symFreq struct {
	freq   uint16
	symbol uint16
}

type simpleModel struct {
	totFreq  uint32
	sentinel symFreq
	f        []symFreq // len = nsym + 1 (extra zero-freq terminator)
}

func newSimpleModel(nsym, maxSym int) *simpleModel {
	m := &simpleModel{
		f: make([]symFreq, nsym+1),
	}
	for i := 0; i < maxSym && i < nsym; i++ {
		m.f[i] = symFreq{freq: 1, symbol: uint16(i)}
	}
	for i := maxSym; i < nsym; i++ {
		m.f[i] = symFreq{freq: 0, symbol: uint16(i)}
	}
	m.f[nsym] = symFreq{freq: 0, symbol: 0} // terminator for normalize loop
	m.totFreq = uint32(maxSym)
	m.sentinel = symFreq{freq: modelMaxFreq, symbol: 0}
	return m
}

func (m *simpleModel) normalize() {
	m.totFreq = 0
	for i := range m.f {
		if m.f[i].freq == 0 {
			break
		}
		m.f[i].freq -= m.f[i].freq >> 1
		m.totFreq += uint32(m.f[i].freq)
	}
}

func (m *simpleModel) decodeSymbol(rc *rangeCoder) uint16 {
	freq := rc.getFreq(m.totFreq)
	if freq > modelMaxFreq {
		return 0
	}

	var accFreq uint32
	i := 0
	for accFreq+uint32(m.f[i].freq) <= freq {
		accFreq += uint32(m.f[i].freq)
		i++
		if i >= len(m.f) {
			return 0
		}
	}

	rc.decode(accFreq, uint32(m.f[i].freq), m.totFreq)
	m.f[i].freq += modelStep
	m.totFreq += modelStep

	if m.totFreq > modelMaxFreq {
		m.normalize()
	}

	// Keep approximately sorted: swap with predecessor if needed.
	if i > 0 && m.f[i].freq > m.f[i-1].freq {
		t := m.f[i]
		m.f[i] = m.f[i-1]
		m.f[i-1] = t
		return t.symbol
	}

	return m.f[i].symbol
}
