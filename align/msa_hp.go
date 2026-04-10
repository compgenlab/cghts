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

	"github.com/compgen-io/cgltk/seqio"
	"github.com/compgen-io/cgltk/support/sequtils"
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

// RehydratedConsensus reconstructs the full-length consensus by expanding
// each compressed-position consensus base by its chosen homopolymer run
// length. It requires HPLens to be populated (i.e., MSA must have been
// called with MSAOptions.HPCompress(true)).
//
// Length-selection rule at each column:
//
//  1. Count HP length frequencies across the non-ref rows with a non-gap
//     base at this column (ignore rows that have a gap — they contribute
//     nothing, matching the spec).
//  2. If the mode is unique, use it.
//  3. If multiple lengths are tied for the mode, fold the reference's HP
//     length into the tally (when a reference is present and non-gap at
//     this column) and look for a new unique mode.
//  4. If still tied, return the ceiling of the mean of the currently-tied
//     mode set. This guarantees a single integer output.
//
// Edge case: if no non-ref row has a non-gap base at a column but the
// reference does, the ref's length is used directly. If no row contributes,
// the column is skipped.
func (p *MSAAlignment) RehydratedConsensus() string {
	if p.HPLens == nil {
		// No HP data — fall back to the plain consensus. Callers who
		// really want rehydration should have enabled HPCompress.
		return p.Consensus()
	}

	// rowCursors[i] tracks the compressed-sequence index into HPLens[i].
	// Every time row i has a non-gap base at a column we consume
	// HPLens[i][rowCursors[i]] and advance the cursor. This is how we map
	// "alignment column k" back to "original HP run k' in row i" without
	// needing to pre-expand anything.
	rowCursors := make([]int, p.NumSeqs)

	var out strings.Builder
	for _, col := range p.Columns {
		// Majority-vote base, ignoring the ref row (see Consensus).
		base := consensusBase(col, p.RefIdx)

		// Gather per-row HP lengths. We still have to advance the cursors
		// for all non-gap rows even when base == 0 (shouldn't normally
		// happen because Columns only contain columns with >=1 non-gap),
		// so we run this loop unconditionally.
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
				// Defensive guard — if the caller hand-built an alignment
				// with mismatched HPLens, don't walk off the end.
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

		if base == 0 {
			continue
		}
		chosen := chooseHPLength(readLens, refLen, hasRef)
		for k := 0; k < chosen; k++ {
			out.WriteByte(base)
		}
	}
	return out.String()
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

// WriteClustal writes the alignment in CLUSTAL interleaved format.
//
// Blocks are 60 columns wide. A conservation line under each block marks
// columns where every displayed row has the same non-gap base with '*'
// (amino-acid similarity groups are not emitted — this is a DNA format).
// Row names are padded to a consistent width so the sequence columns stay
// aligned across blocks.
func (p *MSAAlignment) WriteClustal(w io.Writer) error {
	const blockWidth = 60

	gapped := p.GappedSequences()
	if len(gapped) == 0 {
		return nil
	}

	// Name-column width: max(10, longest name) + 6-space gutter to match
	// the classical CLUSTAL layout.
	nameWidth := 10
	for _, n := range p.Names {
		if len(n) > nameWidth {
			nameWidth = len(n)
		}
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
		for i, seq := range gapped {
			if _, err := fmt.Fprintf(w, "%-*s%s\n", nameWidth, p.Names[i], seq[start:end]); err != nil {
				return err
			}
		}
		// Conservation line: '*' where every row has the same non-gap base,
		// space otherwise. A single gap anywhere in the column suppresses
		// the mark.
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
