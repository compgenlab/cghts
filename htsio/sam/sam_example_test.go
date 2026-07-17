package sam_test

import (
	"fmt"
	"io"
	"strings"

	"github.com/compgenlab/cghts/htsio/sam"
)

// ExampleTextReader demonstrates reading a small SAM document from an
// in-memory string: parsing the header and iterating over the records.
func ExampleTextReader() {
	const doc = "@HD\tVN:1.6\tSO:coordinate\n" +
		"@SQ\tSN:chr1\tLN:1000\n" +
		"read1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tNM:i:0\n" +
		"read2\t16\tchr1\t200\t60\t10M\t*\t0\t0\tTTTTGGGGCC\tIIIIIIIIII\n"

	r, err := sam.NewTextReaderFromReader(io.NopCloser(strings.NewReader(doc)), nil)
	if err != nil {
		panic(err)
	}
	defer r.Close()

	hdr, _ := r.Header()
	fmt.Println("references:", len(hdr.References()))

	for rec, err := range r.Records() {
		if err != nil {
			panic(err)
		}
		fmt.Printf("%s %s:%d flag=%d\n", rec.ReadName, rec.RefName, rec.Pos, rec.Flag)
	}

	// Output:
	// references: 1
	// read1 chr1:100 flag=0
	// read2 chr1:200 flag=16
}
