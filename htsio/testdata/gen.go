//go:build ignore

package main

import (
	"fmt"
	"os"
	"os/exec"

	"github.com/compgen-io/cgkit/htsio"
	"github.com/compgen-io/cgkit/htsio/bam"
)

func main() {
	dir := "htsio/testdata"
	unsorted := dir + "/unsorted.bam"
	sorted := dir + "/test.bam"

	header := htsio.NewSamHeader()
	header.AddLine("@HD\tVN:1.6\tSO:coordinate")
	header.AddLine("@SQ\tSN:chr1\tLN:100000")
	header.AddLine("@SQ\tSN:chr2\tLN:50000")

	seq := "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC"

	records := []*htsio.SamRecord{
		{ReadName: "read1", Flag: 0, RefName: "chr1", Pos: 100, MapQ: 60, Cigar: "50M", RefNext: "*", PosNext: 0, InsertLen: 0, Seq: seq, Qual: "*", Tags: map[string]htsio.SamTag{"RG": {Type: 'Z', Value: "sample1"}}},
		{ReadName: "read2", Flag: 0, RefName: "chr1", Pos: 500, MapQ: 60, Cigar: "50M", RefNext: "*", PosNext: 0, InsertLen: 0, Seq: seq, Qual: "*", Tags: map[string]htsio.SamTag{}},
		{ReadName: "read3", Flag: 0, RefName: "chr1", Pos: 1000, MapQ: 60, Cigar: "50M", RefNext: "*", PosNext: 0, InsertLen: 0, Seq: seq, Qual: "*", Tags: map[string]htsio.SamTag{}},
		{ReadName: "read4", Flag: 0, RefName: "chr1", Pos: 5000, MapQ: 60, Cigar: "50M", RefNext: "*", PosNext: 0, InsertLen: 0, Seq: seq, Qual: "*", Tags: map[string]htsio.SamTag{}},
		{ReadName: "read5", Flag: 0, RefName: "chr2", Pos: 200, MapQ: 60, Cigar: "50M", RefNext: "*", PosNext: 0, InsertLen: 0, Seq: seq, Qual: "*", Tags: map[string]htsio.SamTag{}},
	}

	w, err := bam.NewWriter(unsorted, header)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bam.NewWriter: %v\n", err)
		os.Exit(1)
	}
	for _, rec := range records {
		if err := w.Write(rec); err != nil {
			fmt.Fprintf(os.Stderr, "Write: %v\n", err)
			os.Exit(1)
		}
	}
	if err := w.Close(); err != nil {
		fmt.Fprintf(os.Stderr, "Close: %v\n", err)
		os.Exit(1)
	}

	cmd := exec.Command("samtools", "sort", "-o", sorted, unsorted)
	if out, err := cmd.CombinedOutput(); err != nil {
		fmt.Fprintf(os.Stderr, "samtools sort: %v\n%s\n", err, out)
		os.Exit(1)
	}

	cmd = exec.Command("samtools", "index", sorted)
	if out, err := cmd.CombinedOutput(); err != nil {
		fmt.Fprintf(os.Stderr, "samtools index: %v\n%s\n", err, out)
		os.Exit(1)
	}

	os.Remove(unsorted)
	fmt.Println("Generated test.bam and test.bam.bai")
}
