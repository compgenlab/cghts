package cram

import (
	"bytes"
	"math/rand"
	"testing"
)

func TestRansOrder0Roundtrip(t *testing.T) {
	tests := []struct {
		name string
		data []byte
	}{
		{"simple", []byte("hello world hello world hello")},
		{"single byte", []byte{42}},
		{"two bytes", []byte{0, 1}},
		{"all same", bytes.Repeat([]byte{77}, 100)},
		{"sequential", func() []byte {
			b := make([]byte, 256)
			for i := range b {
				b[i] = byte(i)
			}
			return b
		}()},
		{"random", func() []byte {
			r := rand.New(rand.NewSource(42))
			b := make([]byte, 1000)
			for i := range b {
				b[i] = byte(r.Intn(256))
			}
			return b
		}()},
		{"skewed", func() []byte {
			// Heavy bias toward a few symbols.
			r := rand.New(rand.NewSource(42))
			b := make([]byte, 1000)
			for i := range b {
				v := r.Intn(100)
				if v < 50 {
					b[i] = 'A'
				} else if v < 80 {
					b[i] = 'C'
				} else if v < 95 {
					b[i] = 'G'
				} else {
					b[i] = 'T'
				}
			}
			return b
		}()},
		{"quality scores", func() []byte {
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
			encoded := encodeRans4x8(tt.data, 0)

			// Decode using existing decoder (skip the order byte prefix).
			decoded, err := decodeRans4x8(encoded)
			if err != nil {
				t.Fatalf("decode error: %v", err)
			}

			if !bytes.Equal(decoded, tt.data) {
				t.Errorf("roundtrip mismatch: got %d bytes, want %d bytes", len(decoded), len(tt.data))
				if len(decoded) < 50 && len(tt.data) < 50 {
					t.Errorf("  got:  %v", decoded)
					t.Errorf("  want: %v", tt.data)
				}
			}
		})
	}
}

func TestRansOrder1Roundtrip(t *testing.T) {
	tests := []struct {
		name string
		data []byte
	}{
		{"simple", []byte("hello world hello world hello")},
		{"all same", bytes.Repeat([]byte{77}, 100)},
		{"sequential", func() []byte {
			b := make([]byte, 256)
			for i := range b {
				b[i] = byte(i)
			}
			return b
		}()},
		{"random", func() []byte {
			r := rand.New(rand.NewSource(42))
			b := make([]byte, 1000)
			for i := range b {
				b[i] = byte(r.Intn(256))
			}
			return b
		}()},
		{"quality scores", func() []byte {
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
			encoded := encodeRans4x8(tt.data, 1)

			decoded, err := decodeRans4x8(encoded)
			if err != nil {
				t.Fatalf("decode error: %v", err)
			}

			if !bytes.Equal(decoded, tt.data) {
				t.Errorf("roundtrip mismatch: got %d bytes, want %d bytes", len(decoded), len(tt.data))
				if len(decoded) < 50 && len(tt.data) < 50 {
					t.Errorf("  got:  %v", decoded)
					t.Errorf("  want: %v", tt.data)
				}
			}
		})
	}
}

func TestRansCompression(t *testing.T) {
	// Test that rANS actually compresses well-suited data.
	data := bytes.Repeat([]byte("ACGT"), 1000)
	encoded := encodeRans4x8(data, 0)
	ratio := float64(len(encoded)) / float64(len(data))
	t.Logf("order-0: %d -> %d (%.1f%%)", len(data), len(encoded), ratio*100)

	if ratio > 0.5 {
		t.Errorf("order-0 compression ratio too high: %.1f%%", ratio*100)
	}

	encoded1 := encodeRans4x8(data, 1)
	ratio1 := float64(len(encoded1)) / float64(len(data))
	t.Logf("order-1: %d -> %d (%.1f%%)", len(data), len(encoded1), ratio1*100)
}
