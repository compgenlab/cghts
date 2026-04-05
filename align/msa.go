package align

import (
	"math"
	"strings"
	"sync"

	"github.com/compgen-io/cgltk/seqio"
	"github.com/compgen-io/cgltk/support/utils"
)

// ProfileColumn holds the bases (or gap '-') for each sequence at one
// alignment column.
type ProfileColumn struct {
	Bases    []byte
	HPRunLen []int // nil for now; future HP-compressed mode
}

// Profile represents a multiple sequence alignment.
type Profile struct {
	Names   []string
	Columns []ProfileColumn
	NumSeqs int
}

// NewProfileFromSeq creates a single-sequence profile from a SeqQual.
func NewProfileFromSeq(sq seqio.SeqQual) *Profile {
	seq := sq.Seq()
	cols := make([]ProfileColumn, len(seq))
	for i := range seq {
		cols[i] = ProfileColumn{Bases: []byte{seq[i]}}
	}
	return &Profile{
		Names:   []string{sq.Name()},
		Columns: cols,
		NumSeqs: 1,
	}
}

// GappedSequences returns the aligned sequences with gap characters.
func (p *Profile) GappedSequences() []string {
	builders := make([]strings.Builder, p.NumSeqs)
	for _, col := range p.Columns {
		for i, b := range col.Bases {
			builders[i].WriteByte(b)
		}
	}
	result := make([]string, p.NumSeqs)
	for i := range builders {
		result[i] = builders[i].String()
	}
	return result
}

// Consensus returns the majority-vote consensus sequence.
// At each column, the most common non-gap base is chosen.
// Ties are broken alphabetically. Columns with no non-gap bases are skipped.
func (p *Profile) Consensus() string {
	var buf strings.Builder
	for _, col := range p.Columns {
		b := consensusBase(col)
		if b != 0 {
			buf.WriteByte(b)
		}
	}
	return buf.String()
}

// consensusBase returns the majority non-gap base for a column, or 0 if all gaps.
func consensusBase(col ProfileColumn) byte {
	var counts [256]int
	total := 0
	for _, b := range col.Bases {
		if b != '-' {
			counts[b]++
			total++
		}
	}
	if total == 0 {
		return 0
	}
	bestBase := byte(0)
	bestCount := 0
	// Check ACGT first for deterministic tie-breaking
	for _, b := range []byte("ACGT") {
		if counts[b] > bestCount {
			bestBase = b
			bestCount = counts[b]
		}
	}
	// Fall back to any other base
	if bestBase == 0 {
		for b := byte('A'); b <= byte('Z'); b++ {
			if counts[b] > bestCount {
				bestBase = b
				bestCount = counts[b]
			}
		}
	}
	return bestBase
}

// extendColumn appends a base to a copy of the column.
func extendColumn(col ProfileColumn, base byte) ProfileColumn {
	newBases := make([]byte, len(col.Bases)+1)
	copy(newBases, col.Bases)
	newBases[len(col.Bases)] = base
	return ProfileColumn{Bases: newBases}
}

// newInsertColumn creates a column with gaps for all existing sequences and a base for the new one.
func newInsertColumn(existingSeqs int, base byte) ProfileColumn {
	bases := make([]byte, existingSeqs+1)
	for i := 0; i < existingSeqs; i++ {
		bases[i] = '-'
	}
	bases[existingSeqs] = base
	return ProfileColumn{Bases: bases}
}

// profileFromAlignment converts a pairwise alignment into a 2-sequence profile.
// The query is row 0 and the target is row 1.
func profileFromAlignment(aln *PairwiseAlignment) *Profile {
	cigar := aln.cigarExpanded
	qSeq := aln.Query.Seq()
	tSeq := aln.Target.Seq()

	cols := make([]ProfileColumn, 0, len(cigar)+aln.TargetStart+(len(tSeq)-aln.TargetEnd))
	qPos := 0
	tPos := 0
	cigarIdx := 0

	// Leading soft clips (query bases before aligned region)
	for cigarIdx < len(cigar) && cigar[cigarIdx] == 'S' {
		cols = append(cols, ProfileColumn{Bases: []byte{qSeq[qPos], '-'}})
		qPos++
		cigarIdx++
	}

	// Target bases before aligned region
	for tPos < aln.TargetStart {
		cols = append(cols, ProfileColumn{Bases: []byte{'-', tSeq[tPos]}})
		tPos++
	}

	// Aligned region (M, I, D)
	for cigarIdx < len(cigar) && cigar[cigarIdx] != 'S' {
		switch cigar[cigarIdx] {
		case 'M':
			cols = append(cols, ProfileColumn{Bases: []byte{qSeq[qPos], tSeq[tPos]}})
			qPos++
			tPos++
		case 'I':
			cols = append(cols, ProfileColumn{Bases: []byte{qSeq[qPos], '-'}})
			qPos++
		case 'D':
			cols = append(cols, ProfileColumn{Bases: []byte{'-', tSeq[tPos]}})
			tPos++
		}
		cigarIdx++
	}

	// Target bases after aligned region
	for tPos < len(tSeq) {
		cols = append(cols, ProfileColumn{Bases: []byte{'-', tSeq[tPos]}})
		tPos++
	}

	// Trailing soft clips
	for cigarIdx < len(cigar) && cigar[cigarIdx] == 'S' {
		cols = append(cols, ProfileColumn{Bases: []byte{qSeq[qPos], '-'}})
		qPos++
		cigarIdx++
	}

	return &Profile{
		Names:   []string{aln.Query.Name(), aln.Target.Name()},
		Columns: cols,
		NumSeqs: 2,
	}
}

// addSequence adds a new sequence to the profile given its alignment to the consensus.
// The alignment must have the new sequence as query and the profile's Consensus() as target.
// This relies on the invariant that every column has at least one non-gap base,
// so len(Consensus()) == len(Columns).
func (p *Profile) addSequence(name string, seq seqio.SeqQual, aln *PairwiseAlignment) *Profile {
	cigar := aln.cigarExpanded
	qSeq := seq.Seq()

	newCols := make([]ProfileColumn, 0, len(p.Columns)+seq.Len())
	qPos := 0
	tPos := 0
	cigarIdx := 0

	// Leading soft clips — new columns at the start
	for cigarIdx < len(cigar) && cigar[cigarIdx] == 'S' {
		newCols = append(newCols, newInsertColumn(p.NumSeqs, qSeq[qPos]))
		qPos++
		cigarIdx++
	}

	// Existing columns before the aligned region
	for tPos < aln.TargetStart {
		newCols = append(newCols, extendColumn(p.Columns[tPos], '-'))
		tPos++
	}

	// Aligned region
	for cigarIdx < len(cigar) && cigar[cigarIdx] != 'S' {
		switch cigar[cigarIdx] {
		case 'M':
			newCols = append(newCols, extendColumn(p.Columns[tPos], qSeq[qPos]))
			qPos++
			tPos++
		case 'D':
			newCols = append(newCols, extendColumn(p.Columns[tPos], '-'))
			tPos++
		case 'I':
			newCols = append(newCols, newInsertColumn(p.NumSeqs, qSeq[qPos]))
			qPos++
		}
		cigarIdx++
	}

	// Existing columns after the aligned region
	for tPos < len(p.Columns) {
		newCols = append(newCols, extendColumn(p.Columns[tPos], '-'))
		tPos++
	}

	// Trailing soft clips — new columns at the end
	for cigarIdx < len(cigar) && cigar[cigarIdx] == 'S' {
		newCols = append(newCols, newInsertColumn(p.NumSeqs, qSeq[qPos]))
		qPos++
		cigarIdx++
	}

	newNames := make([]string, len(p.Names)+1)
	copy(newNames, p.Names)
	newNames[len(p.Names)] = name

	return &Profile{
		Names:   newNames,
		Columns: newCols,
		NumSeqs: p.NumSeqs + 1,
	}
}

// alignedLength returns the number of aligned positions (M, I, D — excluding soft clips).
func alignedLength(aln *PairwiseAlignment) int {
	count := 0
	for i := range aln.cigarExpanded {
		if aln.cigarExpanded[i] != 'S' {
			count++
		}
	}
	return count
}

// selectSeedPair finds the pair with highest alignment score.
// Tiebreak: longest aligned length. Still tied: first pair found.
func selectSeedPair(n int, alignments [][]*PairwiseAlignment) (int, int) {
	bestI, bestJ := 0, 1
	bestScore := alignments[0][1].Score
	bestLen := alignedLength(alignments[0][1])

	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			aln := alignments[i][j]
			alnLen := alignedLength(aln)
			if aln.Score > bestScore || (aln.Score == bestScore && alnLen > bestLen) {
				bestI, bestJ = i, j
				bestScore = aln.Score
				bestLen = alnLen
			}
		}
	}
	return bestI, bestJ
}

// MSAOptions controls the MSA algorithm.
type MSAOptions struct {
	alignmentOpts *alignmentOptions
	maxWorkers    int
	verbose       bool
}

// NewMSAOptions creates MSA options with the given alignment settings.
func NewMSAOptions(alignOpts *alignmentOptions) *MSAOptions {
	return &MSAOptions{
		alignmentOpts: alignOpts,
		maxWorkers:    1,
	}
}

// MaxWorkers sets the maximum number of parallel workers for the initial
// all-pairs alignment phase.
func (o *MSAOptions) MaxWorkers(n int) *MSAOptions {
	o.maxWorkers = n
	return o
}

// Verbose enables debug output during alignment.
func (o *MSAOptions) Verbose(v bool) *MSAOptions {
	o.verbose = v
	return o
}

// MSA performs multiple sequence alignment using incremental consensus.
//
// Algorithm:
//  1. Compute all-pairs pairwise alignments (parallelizable).
//  2. Select the seed pair (highest score, tiebreak longest aligned length).
//  3. Build an initial 2-sequence profile from the seed pair.
//  4. Repeatedly: compute consensus, align all remaining sequences to it,
//     pick the best-scoring one, and add it to the profile.
func MSA(seqs []seqio.SeqQual, opts *MSAOptions) *Profile {
	n := len(seqs)
	if n == 0 {
		return nil
	}
	if n == 1 {
		return NewProfileFromSeq(seqs[0])
	}

	globalAligner := NewGlobalAligner(opts.alignmentOpts)

	if n == 2 {
		aln := globalAligner.Align(seqs[0], seqs[1])
		return profileFromAlignment(aln)
	}

	// All-pairs pairwise alignment (upper triangle only) — global alignment
	alignments := make([][]*PairwiseAlignment, n)
	for i := range alignments {
		alignments[i] = make([]*PairwiseAlignment, n)
	}

	sem := utils.NewSemaphore(max(opts.maxWorkers, 1))
	var wg sync.WaitGroup

	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			wg.Add(1)
			sem.Acquire()
			go func() {
				defer wg.Done()
				defer sem.Release()
				alignments[i][j] = globalAligner.Align(seqs[i], seqs[j])
			}()
		}
	}
	wg.Wait()

	// Select seed pair
	seedI, seedJ := selectSeedPair(n, alignments)

	// Build initial profile from seed pair
	profile := profileFromAlignment(alignments[seedI][seedJ])

	// Track which sequences have been incorporated
	incorporated := make([]bool, n)
	incorporated[seedI] = true
	incorporated[seedJ] = true

	// Semi-global aligner for incorporation: query (read) fully aligned,
	// target (consensus) free end gaps — handles truncated reads.
	semiGlobalAligner := NewSemiGlobalAligner(opts.alignmentOpts)

	// Incrementally add remaining sequences
	for added := 2; added < n; added++ {
		cons := profile.Consensus()
		consSeq := seqio.NewStringSeq(cons, "consensus").FullSeq()

		// Align all remaining sequences to the current consensus, pick best
		bestScore := float32(math.Inf(-1))
		bestIdx := -1
		var bestAln *PairwiseAlignment

		for k := 0; k < n; k++ {
			if incorporated[k] {
				continue
			}
			aln := semiGlobalAligner.Align(seqs[k], consSeq)
			if aln.Score > bestScore {
				bestScore = aln.Score
				bestIdx = k
				bestAln = aln
			}
		}

		profile = profile.addSequence(seqs[bestIdx].Name(), seqs[bestIdx], bestAln)
		incorporated[bestIdx] = true
	}

	return profile
}
