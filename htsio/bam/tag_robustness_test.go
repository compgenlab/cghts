package bam

import (
	"strings"
	"testing"

	"github.com/compgen-io/cgkit/htsio"
)

// TestDecodeTagsTruncatedNoPanic feeds deliberately truncated and malformed aux
// data to decodeTags/decodeArrayTag and asserts they never panic or hang. These
// inputs would otherwise be reachable from a corrupt BAM record.
func TestDecodeTagsTruncatedNoPanic(t *testing.T) {
	cases := map[string][]byte{
		"empty":                 {},
		"tag only":              {'X', 'Y'},
		"type but no value c":   {'X', 'Y', 'c'},
		"type but no value i":   {'X', 'Y', 'i'},
		"unterminated Z":        {'X', 'Y', 'Z', 'a', 'b', 'c'},
		"unknown type":          {'X', 'Y', 'Q', 1, 2, 3},
		"array no count":        {'X', 'Y', 'B', 'i'},
		"array truncated elems": {'X', 'Y', 'B', 'i', 0x03, 0x00, 0x00, 0x00, 0x01},
		// Array claiming a huge count with an UNKNOWN element type: must bail,
		// not spin building 2^31 commas.
		"array huge count unknown elem": {'X', 'Y', 'B', 'Q', 0xFF, 0xFF, 0xFF, 0x7F},
		// Array claiming a huge count with a known type but no data: bounded by
		// the (absent) input bytes.
		"array huge count no data": {'X', 'Y', 'B', 'C', 0xFF, 0xFF, 0xFF, 0x7F},
	}
	for name, data := range cases {
		t.Run(name, func(t *testing.T) {
			// Must return without panicking.
			_, _ = decodeTags(data)
		})
	}
}

// TestEncodeOneTagErrors asserts the writer rejects unsupported tag types and
// unparseable values instead of silently emitting a corrupt or zeroed tag.
func TestEncodeOneTagErrors(t *testing.T) {
	cases := []struct {
		name    string
		tag     htsio.SamTag
		wantErr bool
	}{
		{"valid int", htsio.SamTag{Type: 'i', Value: "42"}, false},
		{"valid float", htsio.SamTag{Type: 'f', Value: "1.5"}, false},
		{"valid string", htsio.SamTag{Type: 'Z', Value: "hello"}, false},
		{"valid array", htsio.SamTag{Type: 'B', Value: "i,1,2,3"}, false},
		{"empty array", htsio.SamTag{Type: 'B', Value: "C"}, false},
		{"bad int", htsio.SamTag{Type: 'i', Value: "notanumber"}, true},
		{"bad float", htsio.SamTag{Type: 'f', Value: "xyz"}, true},
		{"unknown type", htsio.SamTag{Type: 'Q', Value: "1"}, true},
		{"array unknown elem", htsio.SamTag{Type: 'B', Value: "Q,1,2"}, true},
		{"array bad elem", htsio.SamTag{Type: 'B', Value: "i,1,bad"}, true},
		{"array empty value", htsio.SamTag{Type: 'B', Value: ""}, true},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			_, err := encodeOneTag(nil, "XY", tc.tag)
			if (err != nil) != tc.wantErr {
				t.Errorf("encodeOneTag(%+v) err = %v, wantErr = %v", tc.tag, err, tc.wantErr)
			}
		})
	}
}

// TestEncodeDecodeTagRoundTrip checks a few representative tags survive an
// encode/decode round trip through the BAM aux format.
func TestEncodeDecodeTagRoundTrip(t *testing.T) {
	in := map[string]htsio.SamTag{
		"NM": {Type: 'i', Value: "5"},
		"AS": {Type: 'i', Value: "-3"},
		"RG": {Type: 'Z', Value: "group1"},
		"XF": {Type: 'f', Value: "2.5"},
		"BC": {Type: 'B', Value: "C,1,2,3"},
	}
	order := []string{"NM", "AS", "RG", "XF", "BC"}
	encoded, err := encodeAuxTags(in, order)
	if err != nil {
		t.Fatalf("encodeAuxTags: %v", err)
	}
	out, gotOrder := decodeTags(encoded)
	if strings.Join(gotOrder, ",") != strings.Join(order, ",") {
		t.Errorf("tag order = %v, want %v", gotOrder, order)
	}
	for k, want := range in {
		got, ok := out[k]
		if !ok {
			t.Errorf("tag %s missing after round trip", k)
			continue
		}
		if got.Value != want.Value {
			t.Errorf("tag %s value = %q, want %q", k, got.Value, want.Value)
		}
	}
}
