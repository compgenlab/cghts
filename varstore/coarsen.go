package varstore

import "fmt"

// Coarsening a gVCF's reference blocks into the spans a ref_blocks table holds.
//
// WHY COARSEN AT ALL. A GATK gVCF breaks a block at every change of GQ, which is
// about 90,000 blocks per megabase on a low-pass genome (docs/coverage-spans.md):
// the source's banding, not the coverage, sets the row count. A store gates on
// depth and never on GQ, so that resolution buys nothing a query can use. Two
// moves recover it:
//
//   - Band on depth instead. Adjacent blocks whose MIN_DP falls in the same
//     depth band merge, and the merged block's MinDP is the smallest of its
//     parts -- so it still bounds every base it covers. Bands are the store's
//     own DepthBands, which is what makes a query gate answer blocks and
//     catalog runs by the same rule.
//   - Tolerate short gaps. Fragmentation comes from brief excursions -- a
//     handful of bases below the gate, a block boundary -- not from real
//     coverage structure. Merging across an uncovered stretch of at most
//     MaxGap bases weakens the claim to "covered at MinDP, with no uncovered
//     run longer than MaxGap", which is why MaxGap is recorded beside the data.
//
// WHAT IS NEVER MERGED ACROSS: a barrier. A variant record -- any record with a
// real alternate, filtered or not -- is a position where the sample was NOT
// reference, or might not have been. Its gap is the whole reason coverage can be
// trusted not to claim reference where the sample carries something, so a
// tolerance that bridged it would do exactly that. The caller marks one with
// Barrier and the open span is closed there, whatever the gap.

// Coarsener merges a sample's reference blocks, streamed in coordinate order,
// into ref_blocks rows.
type Coarsener struct {
	// MinDP drops a block whose floor is below it. A dropped block is uncovered
	// territory at the gate, so it counts toward a gap like any other.
	MinDP int32
	// Bands are the depth boundaries a merge may not cross. Empty puts every
	// block in one band.
	Bands []int32
	// MaxGap is the longest uncovered stretch, in bases, a merge may span. Zero
	// still merges blocks that abut.
	MaxGap int32
	// Passthrough emits every kept block as the source wrote it: no banding, no
	// merging. For targeted and exome work, where blocks are few and fidelity is
	// cheap.
	Passthrough bool

	// Emit receives each finished block, in coordinate order.
	Emit func(RefBlock) error

	open     *RefBlock
	openBand int
}

// Block offers one reference block from the source. Start and End are 1-based
// and inclusive; MinDP is the block's MIN_DP.
func (c *Coarsener) Block(b RefBlock) error {
	if b.End < b.Start {
		return fmt.Errorf("reference block %s:%d-%d ends before it starts", b.Chrom, b.Start, b.End)
	}
	if c.open != nil && SameChrom(c.open.Chrom, b.Chrom) && b.Start <= c.open.End {
		return fmt.Errorf("reference block %s:%d-%d overlaps the one before it (ends %d); "+
			"a gVCF's blocks for one sample must be disjoint and in order",
			b.Chrom, b.Start, b.End, c.open.End)
	}
	if b.MinDP < c.MinDP {
		return nil // uncovered at the gate; a later block may bridge it, a barrier may not
	}
	band := DepthBand(c.Bands, b.MinDP)
	if c.open != nil && c.mergeable(b, band) {
		c.open.End = b.End
		if b.MinDP < c.open.MinDP {
			c.open.MinDP = b.MinDP
		}
		c.open.GQ = minKnown(c.open.GQ, b.GQ)
		return nil
	}
	if err := c.Flush(); err != nil {
		return err
	}
	bb := b
	c.open, c.openBand = &bb, band
	return nil
}

func (c *Coarsener) mergeable(b RefBlock, band int) bool {
	return !c.Passthrough &&
		SameChrom(c.open.Chrom, b.Chrom) &&
		band == c.openBand &&
		b.Start-c.open.End-1 <= c.MaxGap
}

// Barrier closes the open span at a position the sample may not be reference --
// a variant record, or a block that declined to call. Nothing merges across it.
func (c *Coarsener) Barrier() error { return c.Flush() }

// Flush emits the open span, if any. Call it at the end of the input.
func (c *Coarsener) Flush() error {
	if c.open == nil {
		return nil
	}
	b := *c.open
	c.open = nil
	return c.Emit(b)
}

// minKnown is the smaller of two GQ values, treating Missing as no claim rather
// than as the smallest number.
func minKnown(a, b int32) int32 {
	switch {
	case a == Missing:
		return b
	case b == Missing:
		return a
	case b < a:
		return b
	}
	return a
}
