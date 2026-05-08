package codec

import (
	"bytes"
	"math/rand"
	"testing"
)

func TestRansNx16Order0Roundtrip(t *testing.T) {
	tests := []struct {
		name string
		data []byte
	}{
		{"simple", []byte("hello world hello world hello")},
		{"single byte", []byte{42}},
		{"all same", bytes.Repeat([]byte{77}, 100)},
		{"random", func() []byte {
			r := rand.New(rand.NewSource(42))
			b := make([]byte, 1000)
			for i := range b {
				b[i] = byte(r.Intn(256))
			}
			return b
		}()},
		{"dna", bytes.Repeat([]byte("ACGTACGTNNACGT"), 100)},
		{"quality", func() []byte {
			r := rand.New(rand.NewSource(42))
			b := make([]byte, 500)
			for i := range b {
				b[i] = byte(33 + r.Intn(40))
			}
			return b
		}()},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			encoded := EncodeRansNx16(tt.data)
			decoded, err := DecodeRansNx16(encoded)
			if err != nil {
				t.Fatalf("decode error: %v\n  encoded (%d bytes): %v", err, len(encoded), encoded[:min16(50, len(encoded))])
			}
			if !bytes.Equal(decoded, tt.data) {
				t.Errorf("roundtrip mismatch: got %d bytes, want %d", len(decoded), len(tt.data))
				if len(decoded) < 50 && len(tt.data) < 50 {
					t.Errorf("  got:  %v", decoded)
					t.Errorf("  want: %v", tt.data)
				}
			}
		})
	}
}

func TestRansNx16PackRoundtrip(t *testing.T) {
	// Test data with few unique symbols (triggers PACK).
	tests := []struct {
		name string
		data []byte
	}{
		{"binary", bytes.Repeat([]byte{0, 1, 0, 1, 0, 0, 1}, 100)},
		{"dna4", bytes.Repeat([]byte("ACGT"), 250)},
		{"dna5", bytes.Repeat([]byte("ACGTN"), 200)},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			encoded := EncodeRansNx16(tt.data)
			decoded, err := DecodeRansNx16(encoded)
			if err != nil {
				t.Fatalf("decode error: %v", err)
			}
			if !bytes.Equal(decoded, tt.data) {
				t.Errorf("roundtrip mismatch")
			}

			t.Logf("%s: %d -> %d (%.1f%%)", tt.name, len(tt.data), len(encoded),
				100*float64(len(encoded))/float64(len(tt.data)))
		})
	}
}

func TestRansNx16Order1Roundtrip(t *testing.T) {
	tests := []struct {
		name string
		data []byte
	}{
		{"simple", []byte("hello world hello world hello")},
		{"all same", bytes.Repeat([]byte{77}, 100)},
		{"random", func() []byte {
			r := rand.New(rand.NewSource(42))
			b := make([]byte, 1000)
			for i := range b {
				b[i] = byte(r.Intn(256))
			}
			return b
		}()},
		{"repetitive", bytes.Repeat([]byte("ABCABCABC"), 200)},
		{"quality", func() []byte {
			r := rand.New(rand.NewSource(99))
			b := make([]byte, 2000)
			for i := range b {
				b[i] = byte(33 + r.Intn(40))
			}
			return b
		}()},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Encode with order-1 specifically.
			encoded := encodeRansNx16Order1(tt.data)
			decoded, err := DecodeRansNx16(encoded)
			if err != nil {
				t.Fatalf("decode error: %v", err)
			}
			if !bytes.Equal(decoded, tt.data) {
				t.Errorf("roundtrip mismatch: got %d bytes, want %d", len(decoded), len(tt.data))
				if len(decoded) < 100 && len(tt.data) < 100 {
					t.Errorf("  got:  %v", decoded)
					t.Errorf("  want: %v", tt.data)
				}
			}
			t.Logf("%s: %d -> %d (%.1f%%)", tt.name, len(tt.data), len(encoded),
				100*float64(len(encoded))/float64(len(tt.data)))
		})
	}
}

func TestRansNx16CompetitiveRoundtrip(t *testing.T) {
	// Test that EncodeRansNx16 (competitive) always roundtrips correctly,
	// regardless of which method it picks.
	tests := []struct {
		name string
		data []byte
	}{
		{"short", []byte("test")},
		{"dna", bytes.Repeat([]byte("ACGTACGT"), 500)},
		{"random256", func() []byte {
			r := rand.New(rand.NewSource(123))
			b := make([]byte, 5000)
			for i := range b {
				b[i] = byte(r.Intn(256))
			}
			return b
		}()},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			encoded := EncodeRansNx16(tt.data)
			decoded, err := DecodeRansNx16(encoded)
			if err != nil {
				t.Fatalf("decode error: %v", err)
			}
			if !bytes.Equal(decoded, tt.data) {
				t.Fatalf("roundtrip mismatch")
			}
			t.Logf("%s: %d -> %d (%.1f%%) flags=0x%02x", tt.name, len(tt.data), len(encoded),
				100*float64(len(encoded))/float64(len(tt.data)), encoded[0])
		})
	}
}

func min16(a, b int) int {
	if a < b {
		return a
	}
	return b
}
