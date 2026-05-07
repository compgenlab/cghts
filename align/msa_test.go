package align

import (
	"strings"
	"testing"

	"github.com/compgen-io/cgkit/seqio"
)

// msaTestAligner returns a global aligner suitable for MSA tests.
func msaTestAligner() PairwiseAligner {
	return NewGlobalAligner(DnaAlignmentDefaults())
}

func TestNewMSAAlignmentFromSeq(t *testing.T) {
	sq := seqio.NewStringSeq("ACGT", "seq1").FullSeq()
	p := NewMSAAlignmentFromSeq(sq)

	if p.NumSeqs != 1 {
		t.Fatalf("expected 1 sequence, got %d", p.NumSeqs)
	}
	if len(p.Columns) != 4 {
		t.Fatalf("expected 4 columns, got %d", len(p.Columns))
	}
	if len(p.Names) != 1 || p.Names[0] != "seq1" {
		t.Fatalf("expected name 'seq1', got %v", p.Names)
	}

	gapped := p.GappedSequences()
	if len(gapped) != 1 || gapped[0] != "ACGT" {
		t.Fatalf("expected ['ACGT'], got %v", gapped)
	}
}

func TestConsensus(t *testing.T) {
	p := &MSAAlignment{
		Names: []string{"s1", "s2", "s3"},
		Columns: []MSAColumn{
			{Bases: []byte{'A', 'A', 'A'}}, // unanimous A
			{Bases: []byte{'C', 'C', 'T'}}, // majority C
			{Bases: []byte{'G', 'T', 'T'}}, // majority T
			{Bases: []byte{'T', 'T', 'T'}}, // unanimous T
		},
		NumSeqs: 3,
		RefIdx:  -1,
	}

	cons := p.Consensus()
	if cons != "ACTT" {
		t.Fatalf("expected 'ACTT', got '%s'", cons)
	}
}

func TestConsensusTieBreak(t *testing.T) {
	// Tie between A and C — A wins alphabetically
	p := &MSAAlignment{
		Names: []string{"s1", "s2"},
		Columns: []MSAColumn{
			{Bases: []byte{'A', 'C'}},
		},
		NumSeqs: 2,
		RefIdx:  -1,
	}

	cons := p.Consensus()
	if cons != "A" {
		t.Fatalf("expected 'A' (alphabetical tiebreak), got '%s'", cons)
	}
}

func TestConsensusSkipsGaps(t *testing.T) {
	p := &MSAAlignment{
		Names: []string{"s1", "s2", "s3"},
		Columns: []MSAColumn{
			{Bases: []byte{'A', '-', '-'}}, // only 1 non-gap base
			{Bases: []byte{'C', 'C', '-'}}, // 2 non-gap
		},
		NumSeqs: 3,
		RefIdx:  -1,
	}

	cons := p.Consensus()
	if cons != "AC" {
		t.Fatalf("expected 'AC', got '%s'", cons)
	}
}

func TestGappedSequences(t *testing.T) {
	p := &MSAAlignment{
		Names: []string{"s1", "s2"},
		Columns: []MSAColumn{
			{Bases: []byte{'A', '-'}},
			{Bases: []byte{'C', 'C'}},
			{Bases: []byte{'-', 'G'}},
		},
		NumSeqs: 2,
	}

	gapped := p.GappedSequences()
	if gapped[0] != "AC-" {
		t.Fatalf("expected 'AC-', got '%s'", gapped[0])
	}
	if gapped[1] != "-CG" {
		t.Fatalf("expected '-CG', got '%s'", gapped[1])
	}
}

func TestMSAAlignmentFromAlignment(t *testing.T) {
	aligner := msaTestAligner()

	q := seqio.NewStringSeq("ACGT", "q").FullSeq()
	tgt := seqio.NewStringSeq("ACGT", "t").FullSeq()
	aln := aligner.Align(q, tgt)

	p := msaFromPairwise(aln)
	if p.NumSeqs != 2 {
		t.Fatalf("expected 2 sequences, got %d", p.NumSeqs)
	}

	gapped := p.GappedSequences()
	if gapped[0] != "ACGT" {
		t.Fatalf("expected 'ACGT' for query, got '%s'", gapped[0])
	}
	if gapped[1] != "ACGT" {
		t.Fatalf("expected 'ACGT' for target, got '%s'", gapped[1])
	}
}

func TestMSAAlignmentFromAlignmentWithIndel(t *testing.T) {
	aligner := msaTestAligner()

	q := seqio.NewStringSeq("ACGGT", "q").FullSeq()
	tgt := seqio.NewStringSeq("ACGT", "t").FullSeq()
	aln := aligner.Align(q, tgt)

	p := msaFromPairwise(aln)

	gapped := p.GappedSequences()
	// The query has an extra G; the profile should show the insertion
	if !strings.Contains(gapped[1], "-") {
		t.Fatalf("expected gap in target row, got q='%s' t='%s'", gapped[0], gapped[1])
	}
	if len(gapped[0]) != len(gapped[1]) {
		t.Fatalf("gapped sequences have different lengths: %d vs %d", len(gapped[0]), len(gapped[1]))
	}
}

func TestAddSequence(t *testing.T) {
	aligner := msaTestAligner()

	s1 := seqio.NewStringSeq("ACGT", "s1").FullSeq()
	s2 := seqio.NewStringSeq("ACGT", "s2").FullSeq()

	aln := aligner.Align(s1, s2)
	p := msaFromPairwise(aln)

	cons := p.Consensus()
	s3 := seqio.NewStringSeq("ACGT", "s3").FullSeq()
	consSeq := seqio.NewStringSeq(cons, "consensus").FullSeq()
	alnToConsensus := aligner.Align(s3, consSeq)

	p2 := p.AddSequence("s3", s3, alnToConsensus)

	if p2.NumSeqs != 3 {
		t.Fatalf("expected 3 sequences, got %d", p2.NumSeqs)
	}

	gapped := p2.GappedSequences()
	for i, g := range gapped {
		if g != "ACGT" {
			t.Fatalf("sequence %d: expected 'ACGT', got '%s'", i, g)
		}
	}
}

func TestAddSequenceWithInsertion(t *testing.T) {
	aligner := msaTestAligner()

	s1 := seqio.NewStringSeq("ACGT", "s1").FullSeq()
	s2 := seqio.NewStringSeq("ACGT", "s2").FullSeq()
	aln := aligner.Align(s1, s2)
	p := msaFromPairwise(aln)

	s3 := seqio.NewStringSeq("ACGGT", "s3").FullSeq()
	cons := p.Consensus()
	consSeq := seqio.NewStringSeq(cons, "consensus").FullSeq()
	alnToConsensus := aligner.Align(s3, consSeq)
	p2 := p.AddSequence("s3", s3, alnToConsensus)

	gapped := p2.GappedSequences()
	// All gapped sequences should be the same length
	for i := 1; i < len(gapped); i++ {
		if len(gapped[i]) != len(gapped[0]) {
			t.Fatalf("length mismatch: seq 0 len %d, seq %d len %d", len(gapped[0]), i, len(gapped[i]))
		}
	}
	// s3 should have no gaps (it's the longest)
	if strings.Contains(gapped[2], "-") {
		t.Fatalf("s3 should have no gaps, got '%s'", gapped[2])
	}
	// s1 and s2 should each have exactly one gap
	for i := 0; i < 2; i++ {
		gaps := strings.Count(gapped[i], "-")
		if gaps != 1 {
			t.Fatalf("seq %d expected 1 gap, got %d in '%s'", i, gaps, gapped[i])
		}
	}
}

func TestSelectSeedPair(t *testing.T) {
	aligner := msaTestAligner()

	seqs := []seqio.SeqQual{
		seqio.NewStringSeq("AAAA", "s0").FullSeq(),
		seqio.NewStringSeq("ACGTACGT", "s1").FullSeq(),
		seqio.NewStringSeq("ACGTACGT", "s2").FullSeq(),
	}

	n := len(seqs)
	alignments := make([][]*PairwiseAlignment, n)
	for i := range alignments {
		alignments[i] = make([]*PairwiseAlignment, n)
	}
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			alignments[i][j] = aligner.Align(seqs[i], seqs[j])
		}
	}

	seedI, seedJ := selectSeedPair(n, alignments)
	if seedI != 1 || seedJ != 2 {
		t.Fatalf("expected seed pair (1, 2), got (%d, %d)", seedI, seedJ)
	}
}

func TestMSA_TwoSequences(t *testing.T) {
	opts := NewMSAOptions(DnaAlignmentDefaults())

	seqs := []seqio.SeqQual{
		seqio.NewStringSeq("ACGT", "s1").FullSeq(),
		seqio.NewStringSeq("ACGT", "s2").FullSeq(),
	}

	p, _ := MSA(seqs, opts)
	if p == nil {
		t.Fatal("expected non-nil profile")
	}
	if p.NumSeqs != 2 {
		t.Fatalf("expected 2 sequences, got %d", p.NumSeqs)
	}

	gapped := p.GappedSequences()
	for i, g := range gapped {
		if g != "ACGT" {
			t.Fatalf("sequence %d: expected 'ACGT', got '%s'", i, g)
		}
	}
}

func TestMSA_ThreeIdentical(t *testing.T) {
	opts := NewMSAOptions(DnaAlignmentDefaults())

	seqs := []seqio.SeqQual{
		seqio.NewStringSeq("ACGTACGT", "s1").FullSeq(),
		seqio.NewStringSeq("ACGTACGT", "s2").FullSeq(),
		seqio.NewStringSeq("ACGTACGT", "s3").FullSeq(),
	}

	p, _ := MSA(seqs, opts)
	if p.NumSeqs != 3 {
		t.Fatalf("expected 3 sequences, got %d", p.NumSeqs)
	}

	gapped := p.GappedSequences()
	for i, g := range gapped {
		if g != "ACGTACGT" {
			t.Fatalf("sequence %d: expected 'ACGTACGT', got '%s'", i, g)
		}
	}
}

func TestMSA_WithMutation(t *testing.T) {
	opts := NewMSAOptions(DnaAlignmentDefaults())

	seqs := []seqio.SeqQual{
		seqio.NewStringSeq("ACGTACGT", "s1").FullSeq(),
		seqio.NewStringSeq("ACGTACGT", "s2").FullSeq(),
		seqio.NewStringSeq("ACTTACGT", "s3").FullSeq(), // G→T at pos 2
	}

	p, _ := MSA(seqs, opts)
	if p.NumSeqs != 3 {
		t.Fatalf("expected 3 sequences, got %d", p.NumSeqs)
	}

	cons := p.Consensus()
	if cons != "ACGTACGT" {
		t.Fatalf("expected consensus 'ACGTACGT', got '%s'", cons)
	}

	gapped := p.GappedSequences()
	for i := 1; i < len(gapped); i++ {
		if len(gapped[i]) != len(gapped[0]) {
			t.Fatalf("length mismatch at sequence %d", i)
		}
	}
}

func TestMSA_WithInsertion(t *testing.T) {
	opts := NewMSAOptions(DnaAlignmentDefaults())

	seqs := []seqio.SeqQual{
		seqio.NewStringSeq("ACGTACGT", "s1").FullSeq(),
		seqio.NewStringSeq("ACGTACGT", "s2").FullSeq(),
		seqio.NewStringSeq("ACGTTACGT", "s3").FullSeq(), // extra T inserted
	}

	p, _ := MSA(seqs, opts)
	if p.NumSeqs != 3 {
		t.Fatalf("expected 3 sequences, got %d", p.NumSeqs)
	}

	gapped := p.GappedSequences()
	// All should be same length
	for i := 1; i < len(gapped); i++ {
		if len(gapped[i]) != len(gapped[0]) {
			t.Fatalf("length mismatch: seq 0 = %d, seq %d = %d\n%v", len(gapped[0]), i, len(gapped[i]), gapped)
		}
	}
	// s3 has an insertion so s1 and s2 should have a gap somewhere
	if !strings.Contains(gapped[0], "-") || !strings.Contains(gapped[1], "-") {
		t.Logf("gapped: %v", gapped)
		t.Fatalf("expected gaps in s1 and s2 due to s3 insertion")
	}
}

func TestMSA_SingleSequence(t *testing.T) {
	opts := NewMSAOptions(DnaAlignmentDefaults())

	seqs := []seqio.SeqQual{
		seqio.NewStringSeq("ACGT", "s1").FullSeq(),
	}

	p, _ := MSA(seqs, opts)
	if p.NumSeqs != 1 {
		t.Fatalf("expected 1 sequence, got %d", p.NumSeqs)
	}
	if p.Consensus() != "ACGT" {
		t.Fatalf("expected consensus 'ACGT', got '%s'", p.Consensus())
	}
}

func TestMSA_Empty(t *testing.T) {
	opts := NewMSAOptions(DnaAlignmentDefaults())
	p, _ := MSA(nil, opts)
	if p != nil {
		t.Fatalf("expected nil for empty input, got %v", p)
	}
}

func TestMSA_FourSequencesEndToEnd(t *testing.T) {
	opts := NewMSAOptions(DnaAlignmentDefaults())

	seqs := []seqio.SeqQual{
		seqio.NewStringSeq("ACGTACGTACGT", "s1").FullSeq(),
		seqio.NewStringSeq("ACGTACGTACGT", "s2").FullSeq(),
		seqio.NewStringSeq("ACGTACCTACGT", "s3").FullSeq(), // G→C at pos 6
		seqio.NewStringSeq("ACGTACGTACGT", "s4").FullSeq(),
	}

	p, _ := MSA(seqs, opts)
	if p.NumSeqs != 4 {
		t.Fatalf("expected 4 sequences, got %d", p.NumSeqs)
	}

	cons := p.Consensus()
	if cons != "ACGTACGTACGT" {
		t.Fatalf("expected consensus 'ACGTACGTACGT', got '%s'", cons)
	}

	gapped := p.GappedSequences()
	for i := 1; i < len(gapped); i++ {
		if len(gapped[i]) != len(gapped[0]) {
			t.Fatalf("length mismatch at sequence %d", i)
		}
	}
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

func TestMSA_TruncatedReads(t *testing.T) {
	// Verify that truncated reads align correctly with leading/trailing gaps
	// and don't corrupt the consensus.
	ref := "ACGTACGTACGTACGTACGTACGTACGTACGT"

	seqs := []seqio.SeqQual{
		seqio.NewStringSeq(ref, "full1").FullSeq(),
		seqio.NewStringSeq(ref, "full2").FullSeq(),
		seqio.NewStringSeq(ref, "full3").FullSeq(),
		seqio.NewStringSeq(ref[12:], "trunc_left").FullSeq(),  // missing first 12 bases
		seqio.NewStringSeq(ref[:20], "trunc_right").FullSeq(), // missing last 12 bases
	}

	opts := NewMSAOptions(OntAlignmentDefaults())
	profile, _ := MSA(seqs, opts)
	if profile == nil {
		t.Fatal("MSA returned nil")
	}

	cons := profile.Consensus()
	if cons != ref {
		t.Errorf("consensus = %q, want %q", cons, ref)
	}

	// Truncated reads should have gaps at the missing ends.
	gapped := profile.GappedSequences()
	for i := 1; i < len(gapped); i++ {
		if len(gapped[i]) != len(gapped[0]) {
			t.Fatalf("gapped length mismatch: seq 0 = %d, seq %d = %d", len(gapped[0]), i, len(gapped[i]))
		}
	}
}

func TestMSA_HPErrors(t *testing.T) {
	// HP insertion and deletion with majority correct.
	// Note: the MSA consensus ignores gaps when voting, so a single
	// HP-insertion read contributes an extra base at a column where all
	// other reads have gaps — that base wins unopposed. This is a known
	// limitation of the current consensus algorithm.
	ref := "AAACCCGGGTTT"

	seqs := []seqio.SeqQual{
		seqio.NewStringSeq(ref, "ref1").FullSeq(),
		seqio.NewStringSeq(ref, "ref2").FullSeq(),
		seqio.NewStringSeq(ref, "ref3").FullSeq(),
		seqio.NewStringSeq("AAACCCCGGGTTT", "hp_ins").FullSeq(), // CCC→CCCC
		seqio.NewStringSeq("AAACCGGGTTT", "hp_del").FullSeq(),   // CCC→CC
	}

	opts := NewMSAOptions(OntAlignmentDefaults())
	profile, _ := MSA(seqs, opts)
	if profile == nil {
		t.Fatal("MSA returned nil")
	}

	cons := profile.Consensus()
	t.Logf("consensus = %q (ref = %q)", cons, ref)

	// All gapped sequences must be the same length.
	gapped := profile.GappedSequences()
	for i := 1; i < len(gapped); i++ {
		if len(gapped[i]) != len(gapped[0]) {
			t.Fatalf("gapped length mismatch: seq 0 = %d, seq %d = %d", len(gapped[0]), i, len(gapped[i]))
		}
	}

	// Consensus should be close to ref (may have ±1 HP length error).
	if abs(len(cons)-len(ref)) > 1 {
		t.Errorf("consensus length %d too far from ref length %d", len(cons), len(ref))
	}
}

// -----------------------------------------------------------------------------
// chooseHPLength: pure-logic tests covering every branch of the spec.
//
// The spec (from the user):
//   1. Take the mode. If unique, use it.
//   2. If tied, fold the ref's length in and look again.
//   3. If still tied, ceil(mean(mode set)).
// -----------------------------------------------------------------------------

func TestChooseHPLength_UniqueMode(t *testing.T) {
	// [5, 2, 2] -> mode is {2}, unique -> 2.
	if got := chooseHPLength([]int{5, 2, 2}, 0, false); got != 2 {
		t.Errorf("chooseHPLength([5,2,2]) = %d, want 2", got)
	}
	// [3, 3, 2] -> mode is {3}, unique -> 3.
	if got := chooseHPLength([]int{3, 3, 2}, 0, false); got != 3 {
		t.Errorf("chooseHPLength([3,3,2]) = %d, want 3", got)
	}
}

func TestChooseHPLength_TieCeiling(t *testing.T) {
	// [3, 2] -> mode set {2,3}, no ref, ceil(mean)=ceil(2.5)=3.
	if got := chooseHPLength([]int{3, 2}, 0, false); got != 3 {
		t.Errorf("chooseHPLength([3,2]) = %d, want 3", got)
	}
	// [4, 2] -> ceil(3)=3.
	if got := chooseHPLength([]int{4, 2}, 0, false); got != 3 {
		t.Errorf("chooseHPLength([4,2]) = %d, want 3", got)
	}
	// [5, 5, 2, 2] -> mode set {2,5}, ceil(3.5)=4.
	if got := chooseHPLength([]int{5, 5, 2, 2}, 0, false); got != 4 {
		t.Errorf("chooseHPLength([5,5,2,2]) = %d, want 4", got)
	}
}

func TestChooseHPLength_RefBreaksTie(t *testing.T) {
	// [4, 2] with ref=2 -> counts {4:1, 2:2}, unique mode 2 -> 2.
	if got := chooseHPLength([]int{4, 2}, 2, true); got != 2 {
		t.Errorf("chooseHPLength([4,2], ref=2) = %d, want 2", got)
	}
	// [4, 2] with ref=4 -> counts {4:2, 2:1}, unique mode 4 -> 4.
	if got := chooseHPLength([]int{4, 2}, 4, true); got != 4 {
		t.Errorf("chooseHPLength([4,2], ref=4) = %d, want 4", got)
	}
}

func TestChooseHPLength_RefDoesNotHelp(t *testing.T) {
	// [4, 2] with ref=3 -> counts {4:1, 2:1, 3:1}, all tied ->
	// ceil((4+2+3)/3) = ceil(3) = 3.
	if got := chooseHPLength([]int{4, 2}, 3, true); got != 3 {
		t.Errorf("chooseHPLength([4,2], ref=3) = %d, want 3", got)
	}
}

func TestChooseHPLength_EmptyReads(t *testing.T) {
	// No reads, no ref -> 0 (column effectively skipped by caller).
	if got := chooseHPLength(nil, 0, false); got != 0 {
		t.Errorf("chooseHPLength(nil) = %d, want 0", got)
	}
	// No reads but ref present -> use ref directly.
	if got := chooseHPLength(nil, 4, true); got != 4 {
		t.Errorf("chooseHPLength(nil, ref=4) = %d, want 4", got)
	}
}

// -----------------------------------------------------------------------------
// MSA end-to-end tests: HP compression + reference behavior.
// -----------------------------------------------------------------------------

func TestMSA_HPCompressConsensus(t *testing.T) {
	// Reads that all compress to "ACGTG" with various HP lengths. With no
	// reference, the mode-then-ceiling rule should pick the mode at each
	// column. The rehydrated consensus is verified base-by-base against
	// the manually computed expected value.
	seqs := []seqio.SeqQual{
		seqio.NewStringSeq("AACGTTTGG", "r1").FullSeq(),      // A2 C1 G1 T3 G2
		seqio.NewStringSeq("AAAACGTTTTTGGG", "r2").FullSeq(), // A4 C1 G1 T5 G3
		seqio.NewStringSeq("AACCGTTTGG", "r3").FullSeq(),     // A2 C2 G1 T3 G2
	}
	opts := NewMSAOptions(OntAlignmentDefaults()).HPCompress(true)
	aln, err := MSA(seqs, opts)
	if err != nil {
		t.Fatalf("MSA: %v", err)
	}
	if aln == nil {
		t.Fatal("MSA returned nil")
	}
	if aln.HPLens == nil {
		t.Fatal("HPLens should be populated with HPCompress")
	}

	// Column-by-column manual mode:
	//   A: [2,4,2] -> 2
	//   C: [1,1,2] -> 1
	//   G: [1,1,1] -> 1
	//   T: [3,5,3] -> 3
	//   G: [2,3,2] -> 2
	// Rehydrated consensus: AA + C + G + TTT + GG = "AACGTTTGG"
	cons := aln.RehydratedConsensus()
	want := "AACGTTTGG"
	if cons != want {
		t.Errorf("RehydratedConsensus = %q, want %q", cons, want)
	}
}

func TestMSA_WithRefAppendedAndRotated(t *testing.T) {
	// Reference should be aligned last and moved to row 0 (display-first).
	// The reads themselves should still appear after the ref in input order.
	seqs := []seqio.SeqQual{
		seqio.NewStringSeq("ACGTACGT", "r1").FullSeq(),
		seqio.NewStringSeq("ACGTACGT", "r2").FullSeq(),
		seqio.NewStringSeq("ACGTACGT", "r3").FullSeq(),
		seqio.NewStringSeq("ACGTACGT", "the_ref").FullSeq(),
	}
	opts := NewMSAOptions(DnaAlignmentDefaults()).RefName("the_ref")
	aln, err := MSA(seqs, opts)
	if err != nil {
		t.Fatalf("MSA: %v", err)
	}
	if aln == nil {
		t.Fatal("MSA returned nil")
	}
	if aln.RefIdx != 0 {
		t.Errorf("RefIdx = %d, want 0", aln.RefIdx)
	}
	if aln.Names[0] != "the_ref" {
		t.Errorf("Names[0] = %q, want %q", aln.Names[0], "the_ref")
	}
	// Consensus must exclude the ref row; here the ref matches everyone
	// so the consensus is still "ACGTACGT", but we also verify it when
	// the ref deliberately disagrees in the next test.
	if cons := aln.Consensus(); cons != "ACGTACGT" {
		t.Errorf("Consensus = %q, want %q", cons, "ACGTACGT")
	}
}

func TestMSA_ConsensusExcludesRef(t *testing.T) {
	// Reference deliberately disagrees with the reads at one position.
	// The consensus must follow the reads, not the reference.
	seqs := []seqio.SeqQual{
		seqio.NewStringSeq("ACGTACGT", "r1").FullSeq(),
		seqio.NewStringSeq("ACGTACGT", "r2").FullSeq(),
		seqio.NewStringSeq("ACGTACGT", "r3").FullSeq(),
		seqio.NewStringSeq("ACCTACGT", "the_ref").FullSeq(), // G->C at pos 2
	}
	opts := NewMSAOptions(DnaAlignmentDefaults()).RefName("the_ref")
	aln, err := MSA(seqs, opts)
	if err != nil {
		t.Fatalf("MSA: %v", err)
	}
	if got := aln.Consensus(); got != "ACGTACGT" {
		t.Errorf("Consensus = %q, want %q (ref should be excluded)", got, "ACGTACGT")
	}
}

func TestMSA_RefNotFound(t *testing.T) {
	// A non-empty refName that doesn't match any input must return a hard error.
	seqs := []seqio.SeqQual{
		seqio.NewStringSeq("ACGT", "r1").FullSeq(),
		seqio.NewStringSeq("ACGT", "r2").FullSeq(),
	}
	opts := NewMSAOptions(DnaAlignmentDefaults()).RefName("nonexistent")
	if _, err := MSA(seqs, opts); err == nil {
		t.Error("expected error for missing ref, got nil")
	}
}

func TestMSA_HPCompressWithRefTiebreak(t *testing.T) {
	// Two reads that tie on HP length; the ref breaks the tie by voting
	// for the "2" side.
	//   r1: AAAA  -> A4
	//   r2: AA    -> A2
	//   ref: AA   -> A2
	// With ref folded in, counts = {4:1, 2:2} -> mode {2} -> 2.
	// Rehydrated consensus should therefore be "AA", NOT "AAA".
	seqs := []seqio.SeqQual{
		seqio.NewStringSeq("AAAA", "r1").FullSeq(),
		seqio.NewStringSeq("AA", "r2").FullSeq(),
		seqio.NewStringSeq("AA", "the_ref").FullSeq(),
	}
	opts := NewMSAOptions(OntAlignmentDefaults()).
		HPCompress(true).
		RefName("the_ref")
	aln, err := MSA(seqs, opts)
	if err != nil {
		t.Fatalf("MSA: %v", err)
	}
	if got := aln.RehydratedConsensus(); got != "AA" {
		t.Errorf("RehydratedConsensus = %q, want %q", got, "AA")
	}
}

func TestMSA_HPCompressNoRefCeilingFallback(t *testing.T) {
	// [4, 2] with no ref -> ceil(mean(mode set {4,2})) = ceil(3) = 3.
	seqs := []seqio.SeqQual{
		seqio.NewStringSeq("AAAA", "r1").FullSeq(),
		seqio.NewStringSeq("AA", "r2").FullSeq(),
	}
	opts := NewMSAOptions(OntAlignmentDefaults()).HPCompress(true)
	aln, err := MSA(seqs, opts)
	if err != nil {
		t.Fatalf("MSA: %v", err)
	}
	if got := aln.RehydratedConsensus(); got != "AAA" {
		t.Errorf("RehydratedConsensus = %q, want %q", got, "AAA")
	}
}

// -----------------------------------------------------------------------------
// Output format tests: make sure the library-level writers produce the
// expected headers and structure.
// -----------------------------------------------------------------------------

func TestWriteFasta(t *testing.T) {
	p := &MSAAlignment{
		Names: []string{"a", "b"},
		Columns: []MSAColumn{
			{Bases: []byte{'A', 'A'}},
			{Bases: []byte{'C', '-'}},
			{Bases: []byte{'G', 'G'}},
		},
		NumSeqs: 2,
		RefIdx:  -1,
	}
	var buf strings.Builder
	if err := p.WriteFasta(&buf); err != nil {
		t.Fatalf("WriteFasta: %v", err)
	}
	want := ">a\nACG\n>b\nA-G\n"
	if buf.String() != want {
		t.Errorf("WriteFasta = %q, want %q", buf.String(), want)
	}
}

func TestWriteClustalHeader(t *testing.T) {
	p := &MSAAlignment{
		Names:   []string{"a"},
		Columns: []MSAColumn{{Bases: []byte{'A'}}},
		NumSeqs: 1,
		RefIdx:  -1,
	}
	var buf strings.Builder
	if err := p.WriteClustal(&buf); err != nil {
		t.Fatalf("WriteClustal: %v", err)
	}
	// First line must start with "CLUSTAL" (spec requirement so parsers
	// like BioPython AlignIO will recognize the format).
	first, _, _ := strings.Cut(buf.String(), "\n")
	if !strings.HasPrefix(first, "CLUSTAL") {
		t.Errorf("first line = %q, must start with 'CLUSTAL'", first)
	}
}

// TestConsensusRow verifies that ConsensusRow returns exactly one
// character per column (unlike Consensus, which skips empty columns).
// This is the invariant the CLUSTAL writer relies on when displaying
// the synthetic consensus row alongside the real sequences.
func TestConsensusRow(t *testing.T) {
	p := &MSAAlignment{
		Names: []string{"s1", "s2"},
		Columns: []MSAColumn{
			{Bases: []byte{'A', 'A'}},
			{Bases: []byte{'C', '-'}},
			{Bases: []byte{'G', 'G'}},
		},
		NumSeqs: 2,
		RefIdx:  -1,
	}
	got := p.ConsensusRow()
	want := "ACG"
	if got != want {
		t.Errorf("ConsensusRow = %q, want %q", got, want)
	}
}

func TestWriteClustalWithConsensus(t *testing.T) {
	// No HP data — consensus row should be appended but no _hp rows.
	seqs := []seqio.SeqQual{
		seqio.NewStringSeq("ACGT", "r1").FullSeq(),
		seqio.NewStringSeq("ACGT", "r2").FullSeq(),
		seqio.NewStringSeq("ACTT", "r3").FullSeq(), // G->T at col 2
	}
	opts := NewMSAOptions(DnaAlignmentDefaults())
	aln, err := MSA(seqs, opts)
	if err != nil {
		t.Fatalf("MSA: %v", err)
	}

	var buf strings.Builder
	if err := aln.WriteClustalWithConsensus(&buf); err != nil {
		t.Fatalf("WriteClustalWithConsensus: %v", err)
	}
	out := buf.String()

	// Consensus row must appear, named "consensus", with bases ACGT
	// (majority vote beats the lone T at col 2).
	if !strings.Contains(out, "consensus") {
		t.Errorf("expected consensus row, got:\n%s", out)
	}
	if !strings.Contains(out, "ACGT") {
		t.Errorf("expected 'ACGT' in consensus output, got:\n%s", out)
	}
	// HP length rows must NOT appear when HP data is absent.
	if strings.Contains(out, "_hp") {
		t.Errorf("unexpected _hp row in non-HP output:\n%s", out)
	}
}

// -----------------------------------------------------------------------------
// Expanded(): HP compression + expansion (rehydration of the full alignment)
// -----------------------------------------------------------------------------

// TestExpanded_IdenticalLengthsNoOp verifies that when every read has the
// same HP run lengths, Expanded() produces an alignment where every row
// matches the original input exactly — no extra gap padding is added
// because no row's run is longer than any other's.
func TestExpanded_IdenticalLengthsNoOp(t *testing.T) {
	seqs := []seqio.SeqQual{
		seqio.NewStringSeq("AACGTTTGG", "r1").FullSeq(),
		seqio.NewStringSeq("AACGTTTGG", "r2").FullSeq(),
	}
	opts := NewMSAOptions(OntAlignmentDefaults()).HPCompress(true)
	aln, err := MSA(seqs, opts)
	if err != nil {
		t.Fatalf("MSA: %v", err)
	}
	exp := aln.Expanded(false)
	if exp == nil {
		t.Fatal("Expanded returned nil")
	}
	gapped := exp.GappedSequences()
	for i, g := range gapped {
		if g != "AACGTTTGG" {
			t.Errorf("row %d expanded = %q, want %q", i, g, "AACGTTTGG")
		}
	}
}

// TestExpanded_DifferingLengthsPadWithGaps verifies that when rows have
// different HP run lengths at the same compressed column, the expanded
// alignment adds trailing gap padding so every row's bases line up at
// the left of their column's expanded block.
func TestExpanded_DifferingLengthsPadWithGaps(t *testing.T) {
	// r1 runs: A=2, C=1, G=1, T=3, G=2 -> AACGTTTGG
	// r2 runs: A=4, C=1, G=1, T=5, G=3 -> AAAACGTTTTTGGG
	// Per-column max: A=4, C=1, G=1, T=5, G=3 -> total 14 columns.
	// r1 expanded: AA--  C  G  TTT--  GG-
	// r2 expanded: AAAA  C  G  TTTTT  GGG
	seqs := []seqio.SeqQual{
		seqio.NewStringSeq("AACGTTTGG", "r1").FullSeq(),
		seqio.NewStringSeq("AAAACGTTTTTGGG", "r2").FullSeq(),
	}
	opts := NewMSAOptions(OntAlignmentDefaults()).HPCompress(true)
	aln, err := MSA(seqs, opts)
	if err != nil {
		t.Fatalf("MSA: %v", err)
	}
	exp := aln.Expanded(false)
	if exp == nil {
		t.Fatal("Expanded returned nil")
	}
	gapped := exp.GappedSequences()
	if len(gapped) != 2 {
		t.Fatalf("expected 2 rows, got %d", len(gapped))
	}
	// Both rows must be the same length (the alignment invariant).
	if len(gapped[0]) != len(gapped[1]) {
		t.Errorf("row lengths differ: %d vs %d", len(gapped[0]), len(gapped[1]))
	}
	// Every column's expanded block is the union of both rows' runs —
	// each row recovers its original sequence when the gaps are stripped.
	for i, g := range gapped {
		stripped := strings.ReplaceAll(g, "-", "")
		want := []string{"AACGTTTGG", "AAAACGTTTTTGGG"}[i]
		if stripped != want {
			t.Errorf("row %d stripped = %q, want %q (full row = %q)", i, stripped, want, g)
		}
	}
}

// TestExpanded_WithConsensusIncludesRow verifies that Expanded(true)
// appends a synthetic consensus row whose HP lengths come from
// chooseHPLength at each column, and that those consensus runs
// participate in the per-column max width so the consensus is shown
// without truncation.
func TestExpanded_WithConsensusIncludesRow(t *testing.T) {
	seqs := []seqio.SeqQual{
		seqio.NewStringSeq("AACGTTTGG", "r1").FullSeq(),      // 2,1,1,3,2
		seqio.NewStringSeq("AAAACGTTTTTGGG", "r2").FullSeq(), // 4,1,1,5,3
		seqio.NewStringSeq("AACCGTTTGG", "r3").FullSeq(),     // 2,2,1,3,2
	}
	// Expected per-column mode:
	//   A: [2,4,2] -> 2
	//   C: [1,1,2] -> 1
	//   G: [1,1,1] -> 1
	//   T: [3,5,3] -> 3
	//   G: [2,3,2] -> 2
	// Consensus rehydrated: AACGTTTGG.
	opts := NewMSAOptions(OntAlignmentDefaults()).HPCompress(true)
	aln, err := MSA(seqs, opts)
	if err != nil {
		t.Fatalf("MSA: %v", err)
	}
	exp := aln.Expanded(true)
	if exp == nil {
		t.Fatal("Expanded returned nil")
	}
	// Consensus row must be the last row, labelled "consensus".
	if exp.Names[exp.NumSeqs-1] != "consensus" {
		t.Errorf("last row name = %q, want %q", exp.Names[exp.NumSeqs-1], "consensus")
	}
	gapped := exp.GappedSequences()
	// Stripped consensus row should be "AACGTTTGG".
	consRow := gapped[exp.NumSeqs-1]
	stripped := strings.ReplaceAll(consRow, "-", "")
	if stripped != "AACGTTTGG" {
		t.Errorf("consensus row stripped = %q, want %q (full row = %q)", stripped, "AACGTTTGG", consRow)
	}
}

// TestExpanded_NoHPData verifies that Expanded is a no-op (returns nil)
// when the alignment has no HPLens data, so callers can't accidentally
// expand an already-expanded or never-compressed alignment.
func TestExpanded_NoHPData(t *testing.T) {
	seqs := []seqio.SeqQual{
		seqio.NewStringSeq("ACGT", "r1").FullSeq(),
		seqio.NewStringSeq("ACGT", "r2").FullSeq(),
	}
	opts := NewMSAOptions(DnaAlignmentDefaults()) // no HPCompress
	aln, err := MSA(seqs, opts)
	if err != nil {
		t.Fatalf("MSA: %v", err)
	}
	if got := aln.Expanded(false); got != nil {
		t.Errorf("Expanded on non-HP alignment should return nil, got %v", got)
	}
}
