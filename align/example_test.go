package align_test

import (
	"fmt"

	"github.com/compgenlab/hts/align"
	"github.com/compgenlab/hts/seqio"
)

// ExampleNewLocalAligner demonstrates a local alignment of two short DNA
// sequences using the Illumina short-read scoring defaults. The query is fully
// contained in the target, so it aligns as a single run of matches and the
// reported coordinates locate it within the target.
func ExampleNewLocalAligner() {
	query := seqio.NewStringSeq("ACGTACGT", "query").FullSeq()
	target := seqio.NewStringSeq("TTACGTACGTTT", "target").FullSeq()

	aligner := align.NewLocalAligner(align.DnaAlignmentDefaults())
	aln := aligner.Align(query, target)

	fmt.Println("CIGAR:", aln.CIGAR)
	fmt.Printf("target region: %d-%d\n", aln.TargetStart, aln.TargetEnd)
	fmt.Println("matches:", aln.Matches())
	// Output:
	// CIGAR: 8M
	// target region: 2-10
	// matches: 8
}

// ExampleCigarCondense shows converting an expanded, per-base CIGAR string into
// the run-length-encoded form.
func ExampleCigarCondense() {
	fmt.Println(align.CigarCondense("IIMMMMMDMM"))
	// Output: 2I5M1D2M
}

// ExampleCigarExpand shows converting a run-length-encoded CIGAR string back
// into its expanded, per-base form.
func ExampleCigarExpand() {
	expanded, err := align.CigarExpand("2I5M1D2M")
	if err != nil {
		panic(err)
	}
	fmt.Println(expanded)
	// Output: IIMMMMMDMM
}
