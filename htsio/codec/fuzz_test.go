package codec

import (
	"bytes"
	"strings"
	"testing"
)

// The decode-only fuzz targets assert that the decoders never panic, hang, or
// allocate without bound on arbitrary (including hostile) input. The Go fuzzing
// engine reports a failure automatically if the target panics or times out, so
// the body only needs to call the decoder and discard the result.

func FuzzDecodeRans4x8(f *testing.F) {
	f.Add([]byte{})
	f.Add([]byte{0, 1, 2, 3, 4, 5, 6, 7})
	f.Fuzz(func(t *testing.T, data []byte) {
		_, _ = DecodeRans4x8(data)
	})
}

func FuzzDecodeRansNx16(f *testing.F) {
	f.Add([]byte{})
	f.Add([]byte{0, 1, 2, 3, 4, 5, 6, 7})
	f.Fuzz(func(t *testing.T, data []byte) {
		_, _ = DecodeRansNx16(data)
	})
}

func FuzzDecodeArithDynamic(f *testing.F) {
	f.Add([]byte{})
	f.Add([]byte{0, 1, 2, 3, 4, 5, 6, 7})
	f.Fuzz(func(t *testing.T, data []byte) {
		_, _ = DecodeArithDynamic(data)
	})
}

func FuzzDecodeFqzcomp(f *testing.F) {
	f.Add([]byte{})
	f.Add(EncodeFqzcomp([]byte{30, 28, 25, 20, 10}, []int{5}))
	f.Fuzz(func(t *testing.T, data []byte) {
		_, _ = DecodeFqzcomp(data)
	})
}

func FuzzDecodeNameTokenizer(f *testing.F) {
	f.Add([]byte{})
	f.Add(EncodeNameTokenizer([]string{"read1", "read2"}))
	f.Fuzz(func(t *testing.T, data []byte) {
		_, _ = DecodeNameTokenizer(data)
	})
}

// Round-trip fuzz targets assert encode-then-decode is the identity for the
// encoders that have a matching decoder.

func FuzzRoundTripRans4x8(f *testing.F) {
	f.Add([]byte("hello world"), 0)
	f.Add([]byte("AAAAAAAAAABBBBBBBBBB"), 1)
	f.Fuzz(func(t *testing.T, data []byte, order int) {
		if len(data) == 0 {
			return
		}
		order &= 1 // only orders 0 and 1 are supported
		encoded := EncodeRans4x8(data, order)
		if encoded == nil {
			return
		}
		decoded, err := DecodeRans4x8(encoded)
		if err != nil {
			t.Fatalf("DecodeRans4x8 (order %d) failed on round trip: %v", order, err)
		}
		if !bytes.Equal(decoded, data) {
			t.Fatalf("rans4x8 order %d round-trip mismatch: got %d bytes, want %d", order, len(decoded), len(data))
		}
	})
}

func FuzzRoundTripRansNx16(f *testing.F) {
	f.Add([]byte("hello world"))
	f.Add(bytes.Repeat([]byte{7}, 64))
	f.Fuzz(func(t *testing.T, data []byte) {
		if len(data) == 0 {
			return
		}
		encoded := EncodeRansNx16(data)
		if encoded == nil {
			return
		}
		decoded, err := DecodeRansNx16(encoded)
		if err != nil {
			t.Fatalf("DecodeRansNx16 failed on round trip: %v", err)
		}
		if !bytes.Equal(decoded, data) {
			t.Fatalf("ransNx16 round-trip mismatch: got %d bytes, want %d", len(decoded), len(data))
		}
	})
}

func FuzzRoundTripFqzcomp(f *testing.F) {
	f.Add([]byte{30, 28, 25, 20, 10})
	f.Fuzz(func(t *testing.T, quals []byte) {
		if len(quals) == 0 {
			return
		}
		// Constrain to the valid Phred quality range (0..93). The byte 0xFF is
		// the BAM "missing quality" sentinel rather than a literal per-base
		// score, so it is outside fqzcomp's input domain and intentionally not
		// round-tripped here.
		quals = append([]byte(nil), quals...)
		for i := range quals {
			quals[i] %= 94
		}
		// Encode as a single read spanning all the quality bytes.
		encoded := EncodeFqzcomp(quals, []int{len(quals)})
		if encoded == nil {
			return
		}
		decoded, err := DecodeFqzcomp(encoded)
		if err != nil {
			t.Fatalf("DecodeFqzcomp failed on round trip: %v", err)
		}
		if !bytes.Equal(decoded, quals) {
			t.Fatalf("fqzcomp round-trip mismatch: got %d bytes, want %d", len(decoded), len(quals))
		}
	})
}

func FuzzRoundTripNameTokenizer(f *testing.F) {
	f.Add("read1\nread2\nread3")
	f.Fuzz(func(t *testing.T, joined string) {
		// Names are NUL-delimited in the decoded output, so a name may not
		// contain a NUL byte. Derive names from newline-separated input.
		var names []string
		for _, n := range strings.Split(joined, "\n") {
			if n != "" && !strings.ContainsRune(n, 0) {
				names = append(names, n)
			}
		}
		if len(names) == 0 {
			return
		}
		encoded := EncodeNameTokenizer(names)
		if encoded == nil {
			return
		}
		decoded, err := DecodeNameTokenizer(encoded)
		if err != nil {
			t.Fatalf("DecodeNameTokenizer failed on round trip: %v", err)
		}
		want := []byte(strings.Join(names, "\x00") + "\x00")
		if !bytes.Equal(decoded, want) {
			t.Fatalf("name tokenizer round-trip mismatch for %d names", len(names))
		}
	})
}
