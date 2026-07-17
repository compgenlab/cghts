package bam_test

import (
	"bytes"
	"fmt"
	"io"

	"github.com/compgenlab/cghts/htsio"
	"github.com/compgenlab/cghts/htsio/bam"
)

// Example_roundTrip writes a single alignment to an in-memory BAM stream and
// then reads it back, demonstrating that the native writer and reader agree on
// the encoding without any external files.
func Example_roundTrip() {
	// Build a minimal header with one reference sequence.
	header := htsio.NewSamHeader()
	header.AddLine("@HD\tVN:1.6\tSO:unsorted")
	header.AddLine("@SQ\tSN:chr1\tLN:1000")

	// Write one record to an in-memory buffer (BGZF-compressed BAM).
	var buf bytes.Buffer
	w := bam.NewWriterFromWriter(&buf, header)
	rec := &htsio.SamRecord{
		ReadName: "read1",
		Flag:     0,
		RefName:  "chr1",
		Pos:      100, // 1-based, as in SAM
		MapQ:     60,
		Cigar:    "4M",
		RefNext:  "*",
		Seq:      "ACGT",
		Qual:     "IIII",
	}
	if err := w.Write(rec); err != nil {
		panic(err)
	}
	if err := w.Close(); err != nil {
		panic(err)
	}

	// Read the records back.
	r, err := bam.NewReader(io.NopCloser(&buf), "", nil)
	if err != nil {
		panic(err)
	}
	defer r.Close()

	for got, err := range r.Records() {
		if err != nil {
			panic(err)
		}
		fmt.Printf("%s %s:%d %s %s\n", got.ReadName, got.RefName, got.Pos, got.Cigar, got.Seq)
	}

	// Output:
	// read1 chr1:100 4M ACGT
}
