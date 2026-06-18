package sequtils_test

import (
	"fmt"

	"github.com/compgenlab/hts/support/sequtils"
)

func ExampleReverseComplement() {
	fmt.Println(sequtils.ReverseComplement("ACGTN"))
	// Output: NACGT
}

func ExampleDNAMatches() {
	// 'R' is the IUPAC code for "A or G", so it matches 'A'.
	fmt.Println(sequtils.DNAMatches('R', 'A'))
	// 'A' and 'C' have no nucleotide in common.
	fmt.Println(sequtils.DNAMatches('A', 'C'))
	// Output:
	// true
	// false
}
