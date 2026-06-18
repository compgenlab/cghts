package seqio_test

import (
	"fmt"
	"strings"

	"github.com/compgenlab/hts/seqio"
)

// ExampleNewFastaReader streams the records of a FASTA file, printing each
// record's name and full sequence.
func ExampleNewFastaReader() {
	const fasta = ">seq1 first record\n" +
		"ACGTACGT\n" +
		"ACGT\n" +
		">seq2\n" +
		"TTTTGGGG\n"

	r, err := seqio.NewFastaReader(strings.NewReader(fasta))
	if err != nil {
		panic(err)
	}

	for {
		rec, err := r.NextSeq()
		if err != nil {
			break // io.EOF marks the end of the stream
		}
		fmt.Printf("%s: %s\n", rec.Name(), rec.FullSeq().Seq())
	}

	// Output:
	// seq1: ACGTACGTACGT
	// seq2: TTTTGGGG
}

// ExampleNewFastqReader streams a FASTQ file, printing the name, sequence, and
// quality of each record.
func ExampleNewFastqReader() {
	const fastq = "@read1 sample\n" +
		"ACGTACGT\n" +
		"+\n" +
		"IIIIIIII\n" +
		"@read2\n" +
		"TTGG\n" +
		"+\n" +
		"!!!!\n"

	r, err := seqio.NewFastqReader(strings.NewReader(fastq))
	if err != nil {
		panic(err)
	}

	for {
		rec, err := r.NextSeq()
		if err != nil {
			break
		}
		sq := rec.FullSeq()
		fmt.Printf("%s\t%s\t%s\n", rec.Name(), sq.Seq(), sq.Qual())
	}

	// Output:
	// read1	ACGTACGT	IIIIIIII
	// read2	TTGG	!!!!
}
