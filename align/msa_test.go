package align

import (
	"strings"
	"testing"

	"github.com/compgen-io/cgltk/seqio"
)

// msaTestAligner returns a global aligner suitable for MSA tests.
func msaTestAligner() PairwiseAligner {
	return NewGlobalAligner(DnaAlignmentDefaults())
}

func TestNewProfileFromSeq(t *testing.T) {
	sq := seqio.NewStringSeq("ACGT", "seq1").FullSeq()
	p := NewProfileFromSeq(sq)

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
	p := &Profile{
		Names: []string{"s1", "s2", "s3"},
		Columns: []ProfileColumn{
			{Bases: []byte{'A', 'A', 'A'}}, // unanimous A
			{Bases: []byte{'C', 'C', 'T'}}, // majority C
			{Bases: []byte{'G', 'T', 'T'}}, // majority T
			{Bases: []byte{'T', 'T', 'T'}}, // unanimous T
		},
		NumSeqs: 3,
	}

	cons := p.Consensus()
	if cons != "ACTT" {
		t.Fatalf("expected 'ACTT', got '%s'", cons)
	}
}

func TestConsensusTieBreak(t *testing.T) {
	// Tie between A and C — A wins alphabetically
	p := &Profile{
		Names: []string{"s1", "s2"},
		Columns: []ProfileColumn{
			{Bases: []byte{'A', 'C'}},
		},
		NumSeqs: 2,
	}

	cons := p.Consensus()
	if cons != "A" {
		t.Fatalf("expected 'A' (alphabetical tiebreak), got '%s'", cons)
	}
}

func TestConsensusSkipsGaps(t *testing.T) {
	p := &Profile{
		Names: []string{"s1", "s2", "s3"},
		Columns: []ProfileColumn{
			{Bases: []byte{'A', '-', '-'}}, // only 1 non-gap base
			{Bases: []byte{'C', 'C', '-'}}, // 2 non-gap
		},
		NumSeqs: 3,
	}

	cons := p.Consensus()
	if cons != "AC" {
		t.Fatalf("expected 'AC', got '%s'", cons)
	}
}

func TestGappedSequences(t *testing.T) {
	p := &Profile{
		Names: []string{"s1", "s2"},
		Columns: []ProfileColumn{
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

func TestProfileFromAlignment(t *testing.T) {
	aligner := msaTestAligner()

	q := seqio.NewStringSeq("ACGT", "q").FullSeq()
	tgt := seqio.NewStringSeq("ACGT", "t").FullSeq()
	aln := aligner.Align(q, tgt)

	p := profileFromAlignment(aln)
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

func TestProfileFromAlignmentWithIndel(t *testing.T) {
	aligner := msaTestAligner()

	q := seqio.NewStringSeq("ACGGT", "q").FullSeq()
	tgt := seqio.NewStringSeq("ACGT", "t").FullSeq()
	aln := aligner.Align(q, tgt)

	p := profileFromAlignment(aln)

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
	p := profileFromAlignment(aln)

	cons := p.Consensus()
	s3 := seqio.NewStringSeq("ACGT", "s3").FullSeq()
	consSeq := seqio.NewStringSeq(cons, "consensus").FullSeq()
	alnToConsensus := aligner.Align(s3, consSeq)

	p2 := p.addSequence("s3", s3, alnToConsensus)

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
	p := profileFromAlignment(aln)

	s3 := seqio.NewStringSeq("ACGGT", "s3").FullSeq()
	cons := p.Consensus()
	consSeq := seqio.NewStringSeq(cons, "consensus").FullSeq()
	alnToConsensus := aligner.Align(s3, consSeq)
	p2 := p.addSequence("s3", s3, alnToConsensus)

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

	p := MSA(seqs, opts)
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

	p := MSA(seqs, opts)
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

	p := MSA(seqs, opts)
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

	p := MSA(seqs, opts)
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

	p := MSA(seqs, opts)
	if p.NumSeqs != 1 {
		t.Fatalf("expected 1 sequence, got %d", p.NumSeqs)
	}
	if p.Consensus() != "ACGT" {
		t.Fatalf("expected consensus 'ACGT', got '%s'", p.Consensus())
	}
}

func TestMSA_Empty(t *testing.T) {
	opts := NewMSAOptions(DnaAlignmentDefaults())
	p := MSA(nil, opts)
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

	p := MSA(seqs, opts)
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
