package seqanalysis_test

import (
	"fmt"
	"strings"

	seqanalysis "github.com/compgenlab/hts/analysis/seq"
	"github.com/compgenlab/hts/seqio"
)

func ExampleCalcGC() {
	r, err := seqio.NewFastaReader(strings.NewReader(">seq1\nGGCCAATT\n"))
	if err != nil {
		panic(err)
	}
	defer r.Close()

	rec, err := r.NextSeq()
	if err != nil {
		panic(err)
	}

	// 4 of the 8 bases are G or C.
	fmt.Println(seqanalysis.CalcGC(rec))
	// Output: 0.5
}
