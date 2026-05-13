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

func TestRansNx16RLERoundtrip(t *testing.T) {
	tests := []struct {
		name string
		data []byte
	}{
		{"runs", func() []byte {
			// Data with lots of consecutive runs.
			var b []byte
			for i := 0; i < 50; i++ {
				for j := 0; j < 10+i%5; j++ {
					b = append(b, byte(i%4))
				}
			}
			return b
		}()},
		{"quality_runs", func() []byte {
			// Quality scores with repeated values (common in real data).
			var b []byte
			for i := 0; i < 200; i++ {
				q := byte(33 + (i/10)%40)
				for j := 0; j < 3+i%8; j++ {
					b = append(b, q)
				}
			}
			return b
		}()},
		{"mixed", func() []byte {
			// Mix of runs and non-runs.
			var b []byte
			for i := 0; i < 100; i++ {
				b = append(b, byte(i%10))
				for j := 0; j < i%5; j++ {
					b = append(b, byte(i%10))
				}
			}
			return b
		}()},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Test RLE + order-0 directly.
			rle0 := encodeRansNx16WithRLE(tt.data, 0)
			if rle0 == nil {
				t.Fatal("RLE encoding returned nil")
			}
			if rle0[0]&ransOrderRLE == 0 {
				t.Fatal("RLE flag not set")
			}
			decoded, err := DecodeRansNx16(rle0)
			if err != nil {
				t.Fatalf("RLE order-0 decode error: %v", err)
			}
			if !bytes.Equal(decoded, tt.data) {
				t.Errorf("RLE order-0 roundtrip mismatch: got %d bytes, want %d", len(decoded), len(tt.data))
			}
			t.Logf("%s RLE+O0: %d -> %d (%.1f%%)", tt.name, len(tt.data), len(rle0),
				100*float64(len(rle0))/float64(len(tt.data)))

			// Test RLE + order-1 directly.
			rle1 := encodeRansNx16WithRLE(tt.data, 1)
			if rle1 == nil {
				t.Fatal("RLE order-1 encoding returned nil")
			}
			decoded, err = DecodeRansNx16(rle1)
			if err != nil {
				t.Fatalf("RLE order-1 decode error: %v", err)
			}
			if !bytes.Equal(decoded, tt.data) {
				t.Errorf("RLE order-1 roundtrip mismatch: got %d bytes, want %d", len(decoded), len(tt.data))
			}
			t.Logf("%s RLE+O1: %d -> %d (%.1f%%)", tt.name, len(tt.data), len(rle1),
				100*float64(len(rle1))/float64(len(tt.data)))
		})
	}
}

func TestRansNx16CatRoundtrip(t *testing.T) {
	data := []byte("incompressible random data xyz")
	encoded := encodeRansNx16Cat(data)
	if encoded[0]&ransOrderCat == 0 {
		t.Fatal("CAT flag not set")
	}
	decoded, err := DecodeRansNx16(encoded)
	if err != nil {
		t.Fatalf("CAT decode error: %v", err)
	}
	if !bytes.Equal(decoded, data) {
		t.Errorf("CAT roundtrip mismatch")
	}
}

func TestRansNx16StripeRoundtrip(t *testing.T) {
	tests := []struct {
		name string
		data []byte
	}{
		{"interleaved_quality", func() []byte {
			// Simulates byte-interleaved structure.
			b := make([]byte, 1000)
			for i := range b {
				b[i] = byte(33 + (i%4)*10 + i/100)
			}
			return b
		}()},
		{"random", func() []byte {
			r := rand.New(rand.NewSource(77))
			b := make([]byte, 2000)
			for i := range b {
				b[i] = byte(r.Intn(256))
			}
			return b
		}()},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			encoded := encodeRansNx16Stripe(tt.data, 4)
			if encoded == nil {
				t.Fatal("STRIPE encoding returned nil")
			}
			if encoded[0]&ransOrderStripe == 0 {
				t.Fatal("STRIPE flag not set")
			}
			decoded, err := DecodeRansNx16(encoded)
			if err != nil {
				t.Fatalf("STRIPE decode error: %v", err)
			}
			if !bytes.Equal(decoded, tt.data) {
				t.Errorf("STRIPE roundtrip mismatch: got %d bytes, want %d", len(decoded), len(tt.data))
			}
			t.Logf("%s STRIPE: %d -> %d (%.1f%%)", tt.name, len(tt.data), len(encoded),
				100*float64(len(encoded))/float64(len(tt.data)))
		})
	}
}

func min16(a, b int) int {
	if a < b {
		return a
	}
	return b
}
