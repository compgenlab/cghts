package codec

import (
	"bytes"
	"fmt"
	"strings"
	"testing"
)

func TestNameTokenizerRoundtrip(t *testing.T) {
	tests := []struct {
		name  string
		names []string
	}{
		{
			"single name",
			[]string{"read1"},
		},
		{
			"two identical",
			[]string{"read1", "read1"},
		},
		{
			"simple incrementing",
			[]string{"read1", "read2", "read3"},
		},
		{
			"illumina-style",
			[]string{
				"INSTRUMENT:1:FLOWCELL:1:1101:1000:2000",
				"INSTRUMENT:1:FLOWCELL:1:1101:1001:2000",
				"INSTRUMENT:1:FLOWCELL:1:1101:1002:2001",
				"INSTRUMENT:1:FLOWCELL:1:1101:1003:2002",
			},
		},
		{
			"mixed separators",
			[]string{
				"SRR123456.1",
				"SRR123456.2",
				"SRR123456.3",
			},
		},
		{
			"zero-padded",
			[]string{
				"READ_0001_A",
				"READ_0002_A",
				"READ_0003_B",
			},
		},
		{
			"all different structure",
			[]string{
				"alpha123",
				"beta456",
				"gamma789",
			},
		},
		{
			"many reads",
			func() []string {
				names := make([]string, 100)
				for i := range names {
					names[i] = fmt.Sprintf("INSTRUMENT:1:FLOWCELL:1:1101:%d:%d", 1000+i, 2000+i%10)
				}
				return names
			}(),
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			encoded := EncodeNameTokenizer(tt.names)
			if encoded == nil {
				t.Fatal("EncodeNameTokenizer returned nil")
			}

			decoded, err := DecodeNameTokenizer(encoded)
			if err != nil {
				t.Fatalf("DecodeNameTokenizer error: %v\n  encoded len=%d", err, len(encoded))
			}

			// Build expected output: names joined by null terminators.
			expected := []byte(strings.Join(tt.names, "\x00") + "\x00")

			if !bytes.Equal(decoded, expected) {
				t.Errorf("roundtrip mismatch")
				// Parse decoded names for comparison.
				decNames := strings.Split(string(decoded), "\x00")
				if len(decNames) > 0 && decNames[len(decNames)-1] == "" {
					decNames = decNames[:len(decNames)-1]
				}
				for i := 0; i < len(tt.names) && i < len(decNames); i++ {
					if decNames[i] != tt.names[i] {
						t.Errorf("  name %d: got %q, want %q", i, decNames[i], tt.names[i])
						if i > 3 {
							break
						}
					}
				}
				if len(decNames) != len(tt.names) {
					t.Errorf("  got %d names, want %d", len(decNames), len(tt.names))
				}
			}

			t.Logf("%s: %d names, %d -> %d bytes (%.1f%%)", tt.name, len(tt.names),
				len(expected), len(encoded), 100*float64(len(encoded))/float64(len(expected)))
		})
	}
}
