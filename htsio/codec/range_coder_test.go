package codec

import (
	"testing"
)

func TestRangeCoderRoundtrip(t *testing.T) {
	// Encode a few symbols with a simple model.
	m1 := newSimpleModel(256, 10)
	re := newRangeEncoder()
	
	symbols := []uint16{5, 3, 7, 1, 5, 5, 3, 0, 9, 2}
	for _, s := range symbols {
		m1.encodeSymbol(re, s)
	}
	encoded := re.finish()
	t.Logf("encoded %d symbols into %d bytes: %v", len(symbols), len(encoded), encoded)
	
	// Decode.
	m2 := newSimpleModel(256, 10)
	rc := newRangeDecoder(encoded)
	
	for i, want := range symbols {
		got := m2.decodeSymbol(rc)
		if got != want {
			t.Errorf("symbol %d: got %d, want %d", i, got, want)
		}
	}
}

func TestFqzArrayRoundtrip(t *testing.T) {
	// Build a ptab-like array using positionBucket.
	var array [1024]uint
	for i := 0; i < 1024; i++ {
		array[i] = uint(positionBucket(i))
	}
	
	encoded := fqzWriteArray(array[:], 1024)
	t.Logf("ptab encoded: %d bytes", len(encoded))
	
	var decoded [1024]uint
	used := fqzReadArray(encoded, decoded[:], 1024)
	if used < 0 {
		t.Fatalf("fqzReadArray failed")
	}
	if used != len(encoded) {
		t.Errorf("consumed %d bytes, encoded %d", used, len(encoded))
	}
	
	for i := 0; i < 1024; i++ {
		if decoded[i] != array[i] {
			t.Errorf("pos %d: got %d, want %d", i, decoded[i], array[i])
			if i > 5 {
				break
			}
		}
	}
}
