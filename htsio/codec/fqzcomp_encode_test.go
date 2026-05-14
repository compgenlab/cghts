package codec

import (
	"bytes"
	"math/rand"
	"testing"
)

func TestFqzcompRoundtrip(t *testing.T) {
	tests := []struct {
		name        string
		quals       []byte
		readLengths []int
	}{
		{
			"single read",
			[]byte{30, 30, 30, 25, 25, 20, 20, 15, 10, 5},
			[]int{10},
		},
		{
			"two reads same length",
			[]byte{30, 30, 28, 25, 22, 20, 18, 15, 30, 30, 28, 25, 22, 20, 18, 15},
			[]int{8, 8},
		},
		{
			"variable lengths",
			[]byte{30, 25, 20, 15, 10, 30, 28, 26, 24, 22, 20, 18},
			[]int{5, 7},
		},
		{
			"constant quality",
			bytes.Repeat([]byte{30}, 100),
			[]int{50, 50},
		},
		{
			"typical Illumina",
			func() []byte {
				r := rand.New(rand.NewSource(42))
				quals := make([]byte, 300)
				for i := range quals {
					// Typical quality: mostly 30-40, drops at ends.
					base := 35
					if i%100 < 5 || i%100 > 90 {
						base = 20
					}
					quals[i] = byte(base + r.Intn(5) - 2)
				}
				return quals
			}(),
			[]int{100, 100, 100},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			encoded := EncodeFqzcomp(tt.quals, tt.readLengths)
			if encoded == nil {
				t.Fatal("EncodeFqzcomp returned nil")
			}

			decoded, err := DecodeFqzcomp(encoded)
			if err != nil {
				t.Fatalf("DecodeFqzcomp error: %v\n  encoded len=%d", err, len(encoded))
			}

			if !bytes.Equal(decoded, tt.quals) {
				t.Errorf("roundtrip mismatch: got %d bytes, want %d", len(decoded), len(tt.quals))
				if len(decoded) < 30 && len(tt.quals) < 30 {
					t.Errorf("  got:  %v", decoded)
					t.Errorf("  want: %v", tt.quals)
				} else {
					// Find first difference.
					for i := 0; i < len(decoded) && i < len(tt.quals); i++ {
						if decoded[i] != tt.quals[i] {
							t.Errorf("  first diff at %d: got %d, want %d", i, decoded[i], tt.quals[i])
							break
						}
					}
				}
			}

			t.Logf("%s: %d -> %d (%.1f%%)", tt.name, len(tt.quals), len(encoded),
				100*float64(len(encoded))/float64(len(tt.quals)))
		})
	}
}

func TestFqzcompRoundtripLargerFixed(t *testing.T) {
	// Test with deterministic quality values to isolate context issues.
	// Pattern: repeating 30,31,32,33,34 for 200 bases across 2 reads.
	quals := make([]byte, 200)
	for i := range quals {
		quals[i] = byte(30 + i%5)
	}
	readLengths := []int{100, 100}
	
	encoded := EncodeFqzcomp(quals, readLengths)
	decoded, err := DecodeFqzcomp(encoded)
	if err != nil {
		t.Fatalf("decode error: %v", err)
	}
	
	for i := 0; i < len(quals); i++ {
		if decoded[i] != quals[i] {
			t.Errorf("mismatch at %d: got %d, want %d (read %d, pos %d)", 
				i, decoded[i], quals[i], i/100, i%100)
			if i > 20 {
				t.FailNow()
			}
		}
	}
	t.Logf("200 -> %d (%.1f%%)", len(encoded), 100*float64(len(encoded))/200.0)
}
