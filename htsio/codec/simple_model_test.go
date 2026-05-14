package codec

import "testing"

func TestSimpleModelRoundtrip(t *testing.T) {
	syms := []uint16{30, 31, 32, 33, 34, 30}
	
	// Encode.
	m1 := newSimpleModel(256, 35)
	re := newRangeEncoder()
	for _, s := range syms {
		m1.encodeSymbol(re, s)
	}
	encoded := re.finish()
	t.Logf("encoded %d symbols into %d bytes", len(syms), len(encoded))
	
	// Decode.
	m2 := newSimpleModel(256, 35)
	rc := newRangeDecoder(encoded)
	for i, want := range syms {
		got := m2.decodeSymbol(rc)
		if got != want {
			t.Errorf("symbol %d: got %d, want %d", i, got, want)
		}
	}
}
