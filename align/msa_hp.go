package align

// This file holds the helper logic for MSA's homopolymer-compression and
// reference-handling features. The public entry point is align.MSA in msa.go;
// the helpers here are internal implementation details called during the
// MSA pipeline and kept in their own file to keep msa.go focused on the
// alignment algorithm itself.

import (
	"fmt"
	"io"
	"sort"
	"strings"

	"github.com/compgen-io/cgkit/seqio"
	"github.com/compgen-io/cgkit/support/sequtils"
)

// splitRefSeq pulls the reference sequence (by name) out of a sequence list.
// It returns the reads with the reference removed, the reference itself, and
// a flag indicating whether a reference was present.
//
//   - Empty refName: no-op. Returns the original slice and hasRef=false.
//   - Non-empty refName: the matching sequence is removed from the slice
//     (order of the remaining sequences is preserved). Returns hasRef=true.
//   - Non-empty refName with no match: returns a hard error.
func splitRefSeq(seqs []seqio.SeqQual, refName string) (reads []seqio.SeqQual, ref seqio.SeqQual, hasRef bool, err error) {
	if refName == "" {
		return seqs, seqio.SeqQual{}, false, nil
	}
	for i, sq := range seqs {
		if sq.Name() == refName {
			reads = make([]seqio.SeqQual, 0, len(seqs)-1)
			reads = append(reads, seqs[:i]...)
			reads = append(reads, seqs[i+1:]...)
			return reads, sq, true, nil
		}
	}
	return nil, seqio.SeqQual{}, false, fmt.Errorf("reference sequence %q not found in input", refName)
}

// compressOne homopolymer-compresses a single sequence. The returned SeqQual
// holds the collapsed bases (quality is dropped; per-base qualities would
// not line up with collapsed runs and seq-msa does not consume them). The
// returned []int holds the original run length at each compressed position,
// i.e. lens[k] == length of the homopolymer run that the k-th compressed
// base represents.
func compressOne(sq seqio.SeqQual) (seqio.SeqQual, []int) {
	compressed, lens := sequtils.HomopolymerCompress(sq.Seq())
	return seqio.NewStringSeq(compressed, sq.Name()).FullSeq(), lens
}

// compressAll runs compressOne on every sequence in the input slice and
// returns the compressed sequences plus a name-indexed map of run lengths.
// A map (rather than a parallel slice) is used because the incremental MSA
// heuristic reorders rows internally based on pairwise scores — looking up
// the run lengths by sequence name is the easiest way to stitch them back
// together after MSA is done.
//
// Assumes sequence names are unique; duplicates silently overwrite.
func compressAll(seqs []seqio.SeqQual) ([]seqio.SeqQual, map[string][]int) {
	outSeqs := make([]seqio.SeqQual, len(seqs))
	lensByName := make(map[string][]int, len(seqs))
	for i, sq := range seqs {
		outSeqs[i], lensByName[sq.Name()] = compressOne(sq)
	}
	return outSeqs, lensByName
}

// rotateRowToFront returns a new MSAAlignment with row `idx` moved to
// position 0. All other rows retain their relative order. Used after a
// reference sequence is appended to place it at the display-first position.
// HPLens is not reshuffled here because MSA populates HPLens after rotation
// (so this helper can stay HP-agnostic).
//
// A no-op for idx <= 0 or idx >= NumSeqs.
func rotateRowToFront(p *MSAAlignment, idx int) *MSAAlignment {
	if idx <= 0 || idx >= p.NumSeqs {
		return p
	}

	// newOrder[k] = the old row index that should live at new row k.
	// We put idx first, then walk 0..NumSeqs skipping idx.
	newOrder := make([]int, 0, p.NumSeqs)
	newOrder = append(newOrder, idx)
	for i := 0; i < p.NumSeqs; i++ {
		if i != idx {
			newOrder = append(newOrder, i)
		}
	}

	newNames := make([]string, p.NumSeqs)
	for newIdx, oldIdx := range newOrder {
		newNames[newIdx] = p.Names[oldIdx]
	}

	newCols := make([]MSAColumn, len(p.Columns))
	for c, col := range p.Columns {
		newBases := make([]byte, p.NumSeqs)
		for newIdx, oldIdx := range newOrder {
			newBases[newIdx] = col.Bases[oldIdx]
		}
		newCols[c] = MSAColumn{Bases: newBases}
	}

	return &MSAAlignment{
		Names:   newNames,
		Columns: newCols,
		NumSeqs: p.NumSeqs,
		RefIdx:  p.RefIdx, // caller updates RefIdx after rotation
	}
}

// forEachConsensusColumn walks the alignment columns once and calls fn for
// every column with the consensus base (0 if the column has no non-gap
// bases among the non-ref rows) and the chosen HP length for that column.
//
// When HPLens is populated, the chosen length is computed by chooseHPLength
// over the non-ref reads at that column, folding in the reference as a
// last-resort tiebreaker. When HPLens is nil, chosenLen is always 0 and fn
// should ignore it.
//
// This is the single-pass walker used by RehydratedConsensus, ConsensusRow,
// and the HP-length row rendering in WriteClustal — they all need the same
// column-by-column (base, length) derivation, so keeping it in one place
// avoids drift between them.
func (p *MSAAlignment) forEachConsensusColumn(fn func(colIdx int, base byte, chosenLen int)) {
	hasHP := p.HPLens != nil

	// rowCursors[i] tracks the compressed-sequence index into HPLens[i].
	// Every time row i has a non-gap base at a column, we consume
	// HPLens[i][rowCursors[i]] and advance the cursor. This maps
	// "alignment column k" back to "original HP run k' in row i" without
	// needing to pre-expand anything.
	rowCursors := make([]int, p.NumSeqs)

	for c, col := range p.Columns {
		base := consensusBase(col, p.RefIdx)

		chosen := 0
		if hasHP {
			var readLens []int
			var refLen int
			var hasRef bool
			for i, b := range col.Bases {
				if b == '-' {
					continue
				}
				cur := rowCursors[i]
				rowCursors[i]++
				if p.HPLens[i] == nil || cur >= len(p.HPLens[i]) {
					// Defensive guard — mismatched HPLens should not
					// walk off the end.
					continue
				}
				runLen := p.HPLens[i][cur]
				if i == p.RefIdx {
					refLen = runLen
					hasRef = true
					continue
				}
				readLens = append(readLens, runLen)
			}
			if base != 0 {
				chosen = chooseHPLength(readLens, refLen, hasRef)
			}
		}

		fn(c, base, chosen)
	}
}

// RehydratedConsensus reconstructs the full-length consensus by expanding
// each compressed-position consensus base by its chosen homopolymer run
// length. It requires HPLens to be populated (i.e., MSA must have been
// called with MSAOptions.HPCompress(true)).
//
// Length-selection rule at each column (see chooseHPLength):
//
//  1. Count HP length frequencies across the non-ref rows with a non-gap
//     base at this column.
//  2. If the mode is unique, use it.
//  3. If multiple lengths are tied, fold the reference's HP length in and
//     look again.
//  4. If still tied, return the ceiling of the mean of the tied mode set.
func (p *MSAAlignment) RehydratedConsensus() string {
	if p.HPLens == nil {
		// No HP data — fall back to the plain consensus. Callers who
		// really want rehydration should have enabled HPCompress.
		return p.Consensus()
	}
	var out strings.Builder
	p.forEachConsensusColumn(func(_ int, base byte, chosen int) {
		if base == 0 {
			return
		}
		for k := 0; k < chosen; k++ {
			out.WriteByte(base)
		}
	})
	return out.String()
}

// ConsensusRow returns the majority-vote consensus as a string aligned
// one-to-one with Columns: exactly len(Columns) characters, with '-' for
// any column that has no non-gap non-ref base. Unlike Consensus(), this
// preserves positional correspondence with the rest of the alignment so
// it can be shown as an extra row in CLUSTAL-style output.
//
// The reference row (if any) is excluded from the vote — the consensus
// reflects the reads only.
func (p *MSAAlignment) ConsensusRow() string {
	out := make([]byte, len(p.Columns))
	for c, col := range p.Columns {
		b := consensusBase(col, p.RefIdx)
		if b == 0 {
			out[c] = '-'
		} else {
			out[c] = b
		}
	}
	return string(out)
}

// chooseHPLength implements the per-column length-selection rule documented
// on RehydratedConsensus. Extracted as a pure function so the unit tests
// can exercise it independently of a real alignment.
func chooseHPLength(lens []int, refLen int, hasRef bool) int {
	if len(lens) == 0 {
		if hasRef {
			return refLen
		}
		return 0
	}

	counts := make(map[int]int, len(lens)+1)
	for _, L := range lens {
		counts[L]++
	}

	mset := modeSet(counts)
	if len(mset) == 1 {
		return mset[0]
	}

	// Tie: fold the reference in if present and look again.
	if hasRef {
		counts[refLen]++
		mset = modeSet(counts)
		if len(mset) == 1 {
			return mset[0]
		}
	}

	// Still tied: ceiling of the mean of the (possibly updated) mode set.
	// Integer ceiling via (sum + n - 1) / n. Valid because HP lengths are
	// always >= 1 in practice.
	sum := 0
	for _, v := range mset {
		sum += v
	}
	return (sum + len(mset) - 1) / len(mset)
}

// modeSet returns the sorted set of values that share the highest frequency
// in counts. A unique mode produces a single-element slice; a tie produces
// multiple values in ascending order (sorted for deterministic output).
func modeSet(counts map[int]int) []int {
	maxCount := 0
	for _, c := range counts {
		if c > maxCount {
			maxCount = c
		}
	}
	var out []int
	for L, c := range counts {
		if c == maxCount {
			out = append(out, L)
		}
	}
	sort.Ints(out)
	return out
}

// WriteFasta writes the alignment as a gapped multi-sequence FASTA. Each
// row appears as a single FASTA record with its gaps (`-`) intact. Row
// order is preserved from the MSAAlignment; if a reference was set it will
// already be at row 0.
func (p *MSAAlignment) WriteFasta(w io.Writer) error {
	gapped := p.GappedSequences()
	for i, seq := range gapped {
		if _, err := fmt.Fprintf(w, ">%s\n%s\n", p.Names[i], seq); err != nil {
			return err
		}
	}
	return nil
}

// consensusRowName is the label used for the synthetic consensus row
// appended by WriteClustalWithConsensus and Expanded(withConsensus=true).
const consensusRowName = "consensus"

// WriteClustal writes the alignment in CLUSTAL interleaved format.
//
// Blocks are 60 columns wide. A conservation line under each block marks
// columns where every row has the same non-gap base with '*' (amino-acid
// similarity groups are not emitted — this is a DNA format). Row names
// are padded to a consistent width so the sequence columns stay aligned
// across blocks.
func (p *MSAAlignment) WriteClustal(w io.Writer) error {
	return p.writeClustalImpl(w, false)
}

// WriteClustalWithConsensus is like WriteClustal but appends a synthetic
// "consensus" row to the bottom of every block, showing the majority-vote
// base at each column (with '-' for columns that have no non-gap non-ref
// bases). This is not strictly valid CLUSTAL — the consensus is a derived
// sequence, not an input — but it is useful when eyeballing PCR-duplicate
// read clusters to see how the consensus differs from the individual
// reads. The required "CLUSTAL" header line is preserved so standard
// parsers will still accept the file.
//
// The conservation line is computed from the input sequences only; the
// consensus row is not included in the vote.
func (p *MSAAlignment) WriteClustalWithConsensus(w io.Writer) error {
	return p.writeClustalImpl(w, true)
}

// writeClustalImpl is the shared CLUSTAL writer used by both WriteClustal
// and WriteClustalWithConsensus. showConsensus toggles whether a synthetic
// consensus row is appended to each block.
func (p *MSAAlignment) writeClustalImpl(w io.Writer, showConsensus bool) error {
	const blockWidth = 60

	gapped := p.GappedSequences()
	if len(gapped) == 0 {
		return nil
	}

	// When showConsensus is set, precompute the ConsensusRow once so the
	// block loop can slice it. Exactly len(Columns) characters.
	var consensusRow string
	if showConsensus {
		consensusRow = p.ConsensusRow()
	}

	// Name-column width: max of (10, longest real name, consensus label)
	// plus a 6-space gutter to match the classical CLUSTAL layout.
	nameWidth := 10
	for _, n := range p.Names {
		if len(n) > nameWidth {
			nameWidth = len(n)
		}
	}
	if showConsensus && len(consensusRowName) > nameWidth {
		nameWidth = len(consensusRowName)
	}
	nameWidth += 6

	alnLen := len(gapped[0])
	if _, err := fmt.Fprintln(w, "CLUSTAL multiple sequence alignment format"); err != nil {
		return err
	}
	if _, err := fmt.Fprintln(w); err != nil {
		return err
	}
	if _, err := fmt.Fprintln(w); err != nil {
		return err
	}

	for start := 0; start < alnLen; start += blockWidth {
		end := start + blockWidth
		if end > alnLen {
			end = alnLen
		}

		// Input sequence rows.
		for i, seq := range gapped {
			if _, err := fmt.Fprintf(w, "%-*s%s\n", nameWidth, p.Names[i], seq[start:end]); err != nil {
				return err
			}
		}

		// Synthetic consensus row (optional). Emitted before the
		// conservation line so the block visually ends with the stars.
		if showConsensus {
			if _, err := fmt.Fprintf(w, "%-*s%s\n", nameWidth, consensusRowName, consensusRow[start:end]); err != nil {
				return err
			}
		}

		// Conservation line: '*' where every input row has the same
		// non-gap base, space otherwise. A single gap anywhere in the
		// column suppresses the mark. Intentionally computed from the
		// sequences only, not from the consensus row.
		var cons strings.Builder
		for j := start; j < end; j++ {
			first := gapped[0][j]
			match := first != '-'
			for _, seq := range gapped[1:] {
				if seq[j] != first {
					match = false
					break
				}
			}
			if match {
				cons.WriteByte('*')
			} else {
				cons.WriteByte(' ')
			}
		}
		if _, err := fmt.Fprintf(w, "%-*s%s\n", nameWidth, "", cons.String()); err != nil {
			return err
		}
		if _, err := fmt.Fprintln(w); err != nil {
			return err
		}
	}
	return nil
}

// Expanded returns a new MSAAlignment where each compressed column has
// been expanded by the per-row homopolymer run lengths recorded in
// HPLens. The expansion turns the collapsed-bases alignment into a
// full-length alignment that can be compared against the original input
// sequences.
//
// How a single compressed column expands:
//
//   - Determine the maximum run length at this column across the rows
//     that contribute (all non-gap rows, plus the consensus if
//     withConsensus is true and the column has a majority base).
//   - Emit that many sub-columns. In each sub-column, a row's base is
//     repeated while the sub-column index is still within that row's run
//     length; once exhausted, the row gets a gap for the remaining
//     sub-columns.
//
// This produces the "expand" half of HP compress/expand: every input
// sequence ends up back at its original length, with gap padding where
// other rows had longer homopolymer runs. The result's HPLens is nil
// (the bases are already fully expanded), and the reference row index
// is preserved.
//
// If withConsensus is true, a synthetic "consensus" row is added at the
// end of the alignment. Its bases are the ConsensusRow and its HP
// lengths are the per-column chosen lengths from chooseHPLength. These
// also contribute to the per-column max, so every run in the result is
// visible — the consensus row is not silently truncated.
//
// Returns nil if HPLens is not populated (no HP data to expand).
func (p *MSAAlignment) Expanded(withConsensus bool) *MSAAlignment {
	if p.HPLens == nil {
		return nil
	}

	numCols := len(p.Columns)

	// Walk the compressed columns once, recording per-row HP lengths,
	// the consensus base and chosen length, and the per-column max
	// width of the expanded block.
	rowLens := make([][]int, numCols) // rowLens[c][i] = original HP run length for row i at column c; 0 for gaps
	maxLens := make([]int, numCols)   // total sub-columns at expanded column c
	consBases := make([]byte, numCols) // consBases[c] = majority base at column c, 0 if none
	consLens := make([]int, numCols)   // consLens[c] = chosen HP length for the consensus at column c

	// Per-row cursor into HPLens[i]. Advances on every non-gap base
	// encountered, so we can map (row, column) -> original run length.
	rowCursors := make([]int, p.NumSeqs)

	for c, col := range p.Columns {
		rowLens[c] = make([]int, p.NumSeqs)
		maxLen := 0

		// Non-ref read lengths for the consensus-length calculation and
		// the ref length (if any) for the tiebreak.
		var readLens []int
		var refLen int
		var hasRef bool

		for i, b := range col.Bases {
			if b == '-' {
				// Row has a gap at this column — contributes nothing.
				continue
			}
			cur := rowCursors[i]
			rowCursors[i]++
			if p.HPLens[i] == nil || cur >= len(p.HPLens[i]) {
				// Defensive: mismatched HPLens, leave length at 0.
				continue
			}
			runLen := p.HPLens[i][cur]
			rowLens[c][i] = runLen
			if runLen > maxLen {
				maxLen = runLen
			}
			if i == p.RefIdx {
				refLen = runLen
				hasRef = true
				continue
			}
			readLens = append(readLens, runLen)
		}

		// Consensus base + chosen length for this column.
		base := consensusBase(col, p.RefIdx)
		if base != 0 {
			consBases[c] = base
			if withConsensus {
				chosen := chooseHPLength(readLens, refLen, hasRef)
				consLens[c] = chosen
				if chosen > maxLen {
					maxLen = chosen
				}
			}
		}

		maxLens[c] = maxLen
	}

	// Build the expanded columns.
	//
	// For each compressed column c, we emit maxLens[c] sub-columns. At
	// sub-column s of column c, a row's base is:
	//   - '-' if the row had a gap at c (rowLens[c][i] == 0)
	//   - col.Bases[i] if s < rowLens[c][i]
	//   - '-' otherwise (the row's run is exhausted)
	//
	// The consensus row (when requested) follows the same pattern but
	// uses consBases[c] and consLens[c].
	total := 0
	for _, m := range maxLens {
		total += m
	}

	numRows := p.NumSeqs
	if withConsensus {
		numRows++
	}

	newCols := make([]MSAColumn, 0, total)
	for c, col := range p.Columns {
		for s := 0; s < maxLens[c]; s++ {
			bases := make([]byte, numRows)
			for i := 0; i < p.NumSeqs; i++ {
				if s < rowLens[c][i] {
					bases[i] = col.Bases[i]
				} else {
					bases[i] = '-'
				}
			}
			if withConsensus {
				if consBases[c] != 0 && s < consLens[c] {
					bases[p.NumSeqs] = consBases[c]
				} else {
					bases[p.NumSeqs] = '-'
				}
			}
			newCols = append(newCols, MSAColumn{Bases: bases})
		}
	}

	newNames := make([]string, numRows)
	copy(newNames, p.Names)
	if withConsensus {
		newNames[p.NumSeqs] = consensusRowName
	}

	return &MSAAlignment{
		Names:   newNames,
		Columns: newCols,
		NumSeqs: numRows,
		RefIdx:  p.RefIdx, // ref row keeps its position; consensus is appended after
		// HPLens left nil — the bases are already expanded.
	}
}
