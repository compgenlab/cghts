package cram

import (
	"bytes"
	"os"
	"testing"

	"github.com/compgen-io/cgltk/htsio/codec"
)

// stripNewlines removes all 0x0a bytes from data.
func stripNewlines(data []byte) []byte {
	out := make([]byte, 0, len(data))
	for _, b := range data {
		if b != '\n' {
			out = append(out, b)
		}
	}
	return out
}

func TestArithDecompress(t *testing.T) {
	tests := []struct {
		name     string
		compFile string
		rawFile  string
	}{
		{"q4.0", "/tmp/htscodecs_test/arith_q4.0", "/tmp/htscodecs_test/q4"},
		{"q4.1", "/tmp/htscodecs_test/arith_q4.1", "/tmp/htscodecs_test/q4"},
		{"q8.0", "/tmp/htscodecs_test/arith_q8.0", "/tmp/htscodecs_test/q8"},
		{"q8.1", "/tmp/htscodecs_test/arith_q8.1", "/tmp/htscodecs_test/q8"},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			comp, err := os.ReadFile(tt.compFile)
			if err != nil {
				t.Skipf("test data not found: %v", err)
			}
			raw, err := os.ReadFile(tt.rawFile)
			if err != nil {
				t.Skipf("raw data not found: %v", err)
			}
			raw = stripNewlines(raw)

			decoded, err := codec.DecodeArithDynamic(comp)
			if err != nil {
				t.Fatalf("codec.DecodeArithDynamic failed: %v", err)
			}

			if !bytes.Equal(decoded, raw) {
				t.Errorf("decoded data mismatch: got %d bytes, want %d bytes", len(decoded), len(raw))
				for i := 0; i < len(decoded) && i < len(raw); i++ {
					if decoded[i] != raw[i] {
						t.Errorf("first diff at byte %d: got 0x%02x, want 0x%02x", i, decoded[i], raw[i])
						break
					}
				}
			}
		})
	}
}

func TestFqzcompDecompress(t *testing.T) {
	tests := []struct {
		name     string
		compFile string
		rawFile  string
	}{
		{"q4.0", "/tmp/htscodecs_test/fqz_q4.0", "/tmp/htscodecs_test/q4"},
		{"q8.0", "/tmp/htscodecs_test/fqz_q8.0", "/tmp/htscodecs_test/q8"},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			comp, err := os.ReadFile(tt.compFile)
			if err != nil {
				t.Skipf("test data not found: %v", err)
			}
			raw, err := os.ReadFile(tt.rawFile)
			if err != nil {
				t.Skipf("raw data not found: %v", err)
			}
			// fqzcomp stores raw Phred values (0-based).
			// The test vectors were compressed from ASCII quality (Phred+33).
			// Strip newlines and subtract 33 to get raw values for comparison.
			raw = stripNewlines(raw)
			for i := range raw {
				raw[i] -= 33
			}

			decoded, err := codec.DecodeFqzcomp(comp)
			if err != nil {
				t.Fatalf("codec.DecodeFqzcomp failed: %v", err)
			}

			if !bytes.Equal(decoded, raw) {
				t.Errorf("decoded data mismatch: got %d bytes, want %d bytes", len(decoded), len(raw))
				for i := 0; i < len(decoded) && i < len(raw); i++ {
					if decoded[i] != raw[i] {
						t.Errorf("first diff at byte %d: got 0x%02x, want 0x%02x", i, decoded[i], raw[i])
						break
					}
				}
			}
		})
	}
}
