package bgzf_test

import (
	"bytes"
	"fmt"
	"io"

	"github.com/compgenlab/cghts/htsio/bgzf"
)

// Example demonstrates a BGZF write/read round trip over an in-memory buffer.
func Example() {
	var buf bytes.Buffer

	// Compress some data into BGZF blocks.
	w := bgzf.NewWriter(&buf)
	if _, err := io.WriteString(w, "hello, bgzf"); err != nil {
		panic(err)
	}
	if err := w.Close(); err != nil {
		panic(err)
	}

	// Decompress it back.
	r := bgzf.NewReader(&buf)
	got, err := io.ReadAll(r)
	if err != nil {
		panic(err)
	}

	fmt.Printf("%s\n", got)
	// Output: hello, bgzf
}
