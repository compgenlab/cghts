package htsio_test

import (
	"fmt"

	"github.com/compgenlab/hts/htsio"

	// Blank imports register the BAM/SAM/CRAM readers with the auto-detection
	// registry used by NewSamReader.
	_ "github.com/compgenlab/hts/htsio/bam"
	_ "github.com/compgenlab/hts/htsio/sam"
)

// ExampleNewSamReader opens a BAM file (format auto-detected from its magic
// bytes) and iterates over every record using the range-over-func API.
func ExampleNewSamReader() {
	r, err := htsio.NewSamReader("testdata/test.bam")
	if err != nil {
		fmt.Println("open:", err)
		return
	}
	defer r.Close()

	var n int
	for rec, err := range r.Records() {
		if err != nil {
			fmt.Println("read:", err)
			return
		}
		_ = rec
		n++
	}
	fmt.Printf("read %d records\n", n)
}

// ExampleParseRegion shows how a samtools-style region string is converted into
// 0-based half-open coordinates suitable for SamReader.Query.
func ExampleParseRegion() {
	ref, start, end, err := htsio.ParseRegion("chr1:1,000-2,000")
	if err != nil {
		fmt.Println(err)
		return
	}
	fmt.Printf("ref=%s start=%d end=%d\n", ref, start, end)
	// Output: ref=chr1 start=999 end=2000
}
