package codec_test

import (
	"fmt"

	"github.com/compgenlab/hts/htsio/codec"
)

// ExampleEncodeRans4x8 shows a self-contained rANS 4x8 order-0 encode/decode
// round trip on a small byte slice.
func ExampleEncodeRans4x8() {
	original := []byte("the quick brown fox jumps over the lazy dog")

	compressed := codec.EncodeRans4x8(original, codec.Order0)

	decoded, err := codec.DecodeRans4x8(compressed)
	if err != nil {
		fmt.Println("decode error:", err)
		return
	}

	fmt.Println(string(decoded))
	fmt.Println(string(decoded) == string(original))
	// Output:
	// the quick brown fox jumps over the lazy dog
	// true
}
