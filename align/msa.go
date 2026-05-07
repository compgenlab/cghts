package align

import (
	"fmt"
	"math"
	"strings"
	"sync"

	"github.com/compgen-io/cgkit/seqio"
	"github.com/compgen-io/cgkit/support/utils"
)

// MSAColumn holds the bases (or gap '-') for each sequence at one
// alignment column. Parallel to MSAAlignment.Names: Bases[i] is the base
// from row i of the alignment at this column.
type MSAColumn struct {
	Bases []byte
}

// MSAAlignment is the full multiple sequence alignment returned by MSA.
//
// The alignment is stored column-major: each MSAColumn lists one base per
// sequence at that column. Names is parallel to the Bases slices so that
// Names[i] is the row i identifier across every column.
//
// When a reference sequence was supplied via MSAOptions.Ref, RefIdx is the
// row index of that sequence within the alignment; otherwise RefIdx is -1.
// The reference is excluded from consensus voting and is only used for
// display ordering and HP-length tiebreaking in RehydratedConsensus.
//
// When homopolymer compression was enabled via MSAOptions.HPCompress, HPLens
// is populated: HPLens[i] holds the original run length for each compressed
// position of row i (before compression). HPLens is parallel to Names. Rows
// with no HP data have a nil HPLens entry.
type MSAAlignment struct {
	Names   []string
	Columns []MSAColumn
	NumSeqs int

	// RefIdx is the row index of the reference sequence, or -1 if no
	// reference was supplied.
	RefIdx int

	// HPLens holds the per-row homopolymer run lengths at each compressed
	// position. Populated only when MSAOptions.HPCompress was set.
	HPLens [][]int
}

// NewMSAAlignmentFromSeq creates a single-row alignment from a SeqQual.
// Useful as a starting point when you only have one sequence (trivial MSA)
// or when you want to construct an alignment manually for testing.
func NewMSAAlignmentFromSeq(sq seqio.SeqQual) *MSAAlignment {
	seq := sq.Seq()
	cols := make([]MSAColumn, len(seq))
	for i := range seq {
		cols[i] = MSAColumn{Bases: []byte{seq[i]}}
	}
	return &MSAAlignment{
		Names:   []string{sq.Name()},
		Columns: cols,
		NumSeqs: 1,
		RefIdx:  -1,
	}
}

// GappedSequences returns the aligned sequences with gap characters.
func (p *MSAAlignment) GappedSequences() []string {
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

// Consensus returns the majority-vote consensus sequence across the read
// rows of the alignment. When a reference sequence is present (RefIdx >= 0)
// it is excluded from the vote so the consensus represents the reads, not
// the reference.
//
// At each column, the most common non-gap base from the non-ref rows is
// chosen. Ties are broken by preferring ACGT in canonical order. Columns
// with no non-gap non-ref bases are skipped entirely (they do not appear in
// the output string), matching the invariant that len(Consensus()) equals
// the number of non-empty columns across the read rows.
func (p *MSAAlignment) Consensus() string {
	var buf strings.Builder
	for _, col := range p.Columns {
		b := consensusBase(col, p.RefIdx)
		if b != 0 {
			buf.WriteByte(b)
		}
	}
	return buf.String()
}

// consensusBase returns the majority non-gap base for a column, ignoring the
// row at excludeIdx (typically a reference row). Returns 0 if no non-gap
// bases exist outside the excluded row.
//
// A negative excludeIdx disables exclusion (all rows vote).
func consensusBase(col MSAColumn, excludeIdx int) byte {
	var counts [256]int
	total := 0
	for i, b := range col.Bases {
		if i == excludeIdx {
			continue
		}
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
	// Fall back to any other non-gap character
	if bestBase == 0 {
		for b := 1; b < 256; b++ {
			if byte(b) != '-' && counts[b] > bestCount {
				bestBase = byte(b)
				bestCount = counts[b]
			}
		}
	}
	return bestBase
}

// extendColumn appends a base to a copy of the column.
func extendColumn(col MSAColumn, base byte) MSAColumn {
	newBases := make([]byte, len(col.Bases)+1)
	copy(newBases, col.Bases)
	newBases[len(col.Bases)] = base
	return MSAColumn{Bases: newBases}
}

// newInsertColumn creates a column with gaps for all existing sequences and a base for the new one.
func newInsertColumn(existingSeqs int, base byte) MSAColumn {
	bases := make([]byte, existingSeqs+1)
	for i := 0; i < existingSeqs; i++ {
		bases[i] = '-'
	}
	bases[existingSeqs] = base
	return MSAColumn{Bases: bases}
}

// msaFromPairwise converts a pairwise alignment into a 2-sequence profile.
// The query is row 0 and the target is row 1.
func msaFromPairwise(aln *PairwiseAlignment) *MSAAlignment {
	cigar := aln.cigarExpanded
	qSeq := aln.Query.Seq()
	tSeq := aln.Target.Seq()

	cols := make([]MSAColumn, 0, len(cigar)+aln.TargetStart+(len(tSeq)-aln.TargetEnd))
	qPos := 0
	tPos := 0
	cigarIdx := 0

	// Leading soft clips (query bases before aligned region)
	for cigarIdx < len(cigar) && cigar[cigarIdx] == 'S' {
		cols = append(cols, MSAColumn{Bases: []byte{qSeq[qPos], '-'}})
		qPos++
		cigarIdx++
	}

	// Target bases before aligned region
	for tPos < aln.TargetStart {
		cols = append(cols, MSAColumn{Bases: []byte{'-', tSeq[tPos]}})
		tPos++
	}

	// Aligned region (M, I, D)
	for cigarIdx < len(cigar) && cigar[cigarIdx] != 'S' {
		switch cigar[cigarIdx] {
		case 'M':
			cols = append(cols, MSAColumn{Bases: []byte{qSeq[qPos], tSeq[tPos]}})
			qPos++
			tPos++
		case 'I':
			cols = append(cols, MSAColumn{Bases: []byte{qSeq[qPos], '-'}})
			qPos++
		case 'D':
			cols = append(cols, MSAColumn{Bases: []byte{'-', tSeq[tPos]}})
			tPos++
		}
		cigarIdx++
	}

	// Target bases after aligned region
	for tPos < len(tSeq) {
		cols = append(cols, MSAColumn{Bases: []byte{'-', tSeq[tPos]}})
		tPos++
	}

	// Trailing soft clips
	for cigarIdx < len(cigar) && cigar[cigarIdx] == 'S' {
		cols = append(cols, MSAColumn{Bases: []byte{qSeq[qPos], '-'}})
		qPos++
		cigarIdx++
	}

	return &MSAAlignment{
		Names:   []string{aln.Query.Name(), aln.Target.Name()},
		Columns: cols,
		NumSeqs: 2,
		RefIdx:  -1,
	}
}

// AddSequence adds a new sequence to the profile given its alignment to the consensus.
// The alignment must have the new sequence as query and the profile's Consensus() as target.
// This relies on the invariant that every column has at least one non-gap base,
// so len(Consensus()) == len(Columns).
//
// Insertions in the new sequence (I operations in the CIGAR) expand the MSA
// with new columns in which all pre-existing rows receive gap bases. This
// makes AddSequence safe for appending a reference sequence that may carry
// bases the existing rows lacked.
func (p *MSAAlignment) AddSequence(name string, seq seqio.SeqQual, aln *PairwiseAlignment) *MSAAlignment {
	cigar := aln.cigarExpanded
	qSeq := seq.Seq()

	newCols := make([]MSAColumn, 0, len(p.Columns)+seq.Len())
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

	return &MSAAlignment{
		Names:   newNames,
		Columns: newCols,
		NumSeqs: p.NumSeqs + 1,
		RefIdx:  p.RefIdx, // propagate ref index unchanged (new row is appended at end)
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
//
// The options are built with a fluent API: start from NewMSAOptions and
// chain the setters. Unset options fall back to sensible defaults
// (no HP compression, no reference, single-threaded).
type MSAOptions struct {
	alignmentOpts *alignmentOptions
	maxWorkers    int
	verbose       bool

	// hpCompress enables homopolymer compression of every input sequence
	// before alignment. The resulting MSAAlignment retains the original
	// run lengths in HPLens, which RehydratedConsensus uses to reconstruct
	// a full-length consensus.
	hpCompress bool

	// refName is the name of the reference sequence within the input set.
	// The reference is removed from the read pool, the MSA is built over
	// the reads only, and the reference is globally aligned to the
	// consensus and appended afterwards. RefIdx on the resulting
	// MSAAlignment is set to 0 (the ref is rotated to the front).
	refName string
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

// HPCompress enables homopolymer compression. When set, each input sequence
// is collapsed to single bases before alignment and the original run lengths
// are retained on the returned MSAAlignment.HPLens so the consensus can be
// rehydrated with the per-column mode length.
func (o *MSAOptions) HPCompress(v bool) *MSAOptions {
	o.hpCompress = v
	return o
}

// RefName marks one of the input sequences (by name) as the reference
// sequence. The reference is removed from the read pool, the MSA is built
// without it, and the reference is globally aligned to the consensus and
// appended at row 0 (display-first position). Its bases are ignored by
// Consensus and are only used as a last-resort tiebreak in
// RehydratedConsensus.
//
// If refName is empty, no reference is used. If refName is non-empty but no
// input sequence matches, MSA returns an error.
func (o *MSAOptions) RefName(refName string) *MSAOptions {
	o.refName = refName
	return o
}

// MSA performs multiple sequence alignment on the input sequences and
// returns the resulting MSAAlignment. It handles all of the feature flags
// on MSAOptions — HP compression and reference handling — inside the
// library so callers get a single consistent result object regardless of
// which options are enabled.
//
// High-level algorithm:
//
//  1. If MSAOptions.Ref is set, split the reference out of the input so the
//     reads can be aligned on their own. Error out if the ref name is not
//     found.
//  2. If MSAOptions.HPCompress is set, collapse homopolymer runs on every
//     sequence (reads and ref). Keep the per-position run lengths on the
//     side so they can be restored later.
//  3. Run the incremental-consensus MSA over the read sequences:
//     a. Compute all-pairs pairwise alignments (parallelizable).
//     b. Pick the seed pair (highest score, tiebreak longest aligned length).
//     c. Build the 2-row starting alignment from the seed.
//     d. Repeatedly compute the consensus, align each remaining read to it
//        semi-globally, and append the best-scoring read.
//  4. If a reference was supplied, globally align it to the finished
//     consensus and append it with AddSequence. Global alignment is used
//     because the ref is expected to span the full locus; we want every
//     ref base mapped to a column. Rotate the ref row to index 0 so
//     downstream display code can emit it first.
//  5. Populate HPLens and RefIdx on the result so the caller can produce
//     rehydrated / reference-aware output without having to re-derive any
//     of this information.
//
// MSA returns (nil, nil) on empty input — callers should treat that as
// "no sequences to align" rather than an error.
func MSA(seqs []seqio.SeqQual, opts *MSAOptions) (*MSAAlignment, error) {
	if len(seqs) == 0 {
		return nil, nil
	}

	// Phase 1: split out the reference sequence if one was named.
	readSeqs, refSeq, hasRef, err := splitRefSeq(seqs, opts.refName)
	if err != nil {
		return nil, err
	}
	if len(readSeqs) == 0 {
		return nil, fmt.Errorf("no read sequences to align (reference cannot be the only input)")
	}

	// Phase 2: homopolymer-compress every sequence if requested. We keep
	// the HP length arrays in parallel slices keyed off the read's name
	// so we can look them up after align.MSA reorders the rows
	// internally (the seed-pair heuristic does not preserve input order).
	var hpLensByName map[string][]int
	var refHPLens []int
	if opts.hpCompress {
		readSeqs, hpLensByName = compressAll(readSeqs)
		if hasRef {
			refSeq, refHPLens = compressOne(refSeq)
		}
	}

	// Phase 3: run the reads-only incremental MSA.
	profile := runIncrementalMSA(readSeqs, opts)
	if profile == nil {
		return nil, nil
	}

	// Phase 4: optionally attach the reference and rotate it to row 0.
	refIdx := -1
	if hasRef {
		profile = attachReference(profile, refSeq, opts.alignmentOpts)
		profile = rotateRowToFront(profile, profile.NumSeqs-1)
		refIdx = 0
	}
	profile.RefIdx = refIdx

	// Phase 5: populate HPLens in row order if HP compression was on.
	if opts.hpCompress {
		profile.HPLens = make([][]int, profile.NumSeqs)
		for i, name := range profile.Names {
			if refIdx == i {
				profile.HPLens[i] = refHPLens
				continue
			}
			profile.HPLens[i] = hpLensByName[name]
		}
	}

	return profile, nil
}

// runIncrementalMSA is the core reads-only MSA loop. It runs over an already
// filtered sequence list (no reference, optionally HP-compressed) and
// returns the resulting alignment without populating HPLens or RefIdx.
func runIncrementalMSA(seqs []seqio.SeqQual, opts *MSAOptions) *MSAAlignment {
	n := len(seqs)
	if n == 0 {
		return nil
	}
	if n == 1 {
		return NewMSAAlignmentFromSeq(seqs[0])
	}

	globalAligner := NewGlobalAligner(opts.alignmentOpts)

	if n == 2 {
		aln := globalAligner.Align(seqs[0], seqs[1])
		return msaFromPairwise(aln)
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

	// Build initial alignment from seed pair
	profile := msaFromPairwise(alignments[seedI][seedJ])

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

		profile = profile.AddSequence(seqs[bestIdx].Name(), seqs[bestIdx], bestAln)
		incorporated[bestIdx] = true
	}

	return profile
}

// attachReference runs a global alignment of the reference sequence against
// the current alignment's consensus and appends it as a new row. Global
// alignment is used (rather than semi-global) because the reference is
// expected to cover the entire locus — every ref base should map to a
// column. Ref-side insertions (bases the reads lacked) create brand-new
// MSA columns via the standard AddSequence path.
func attachReference(profile *MSAAlignment, refSeq seqio.SeqQual, alignOpts *alignmentOptions) *MSAAlignment {
	cons := profile.Consensus()
	consSeq := seqio.NewStringSeq(cons, "consensus").FullSeq()
	globalAligner := NewGlobalAligner(alignOpts)
	aln := globalAligner.Align(refSeq, consSeq)
	return profile.AddSequence(refSeq.Name(), refSeq, aln)
}
