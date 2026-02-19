package align

import (
	"fmt"
	"math"

	"github.com/compgen-io/cgltk/sequtils"
)

type pairwise struct {
	opts    *alignmentOptions
	isLocal bool
}

type swCellTrace uint8

const (
	swTraceStop swCellTrace = iota
	swTraceMatch
	swTraceIns
	swTraceDel
	swTraceClip
)

type swCell struct {
	scoreM float32
	scoreI float32
	scoreD float32
	traceM swCellTrace
	traceI swCellTrace
	traceD swCellTrace
}

/*
This is a simple Smith Waterman local aligner with affine gaps penalties.

It computes a full n*m matrix and backtracks to find an alignment. There are
many optimizations that are possible, but this is the old-school n*m matrix
DP approach.

There are two things that do make this particular method a bit more complex...

First, there are affine gaps. This means that the gap open/extension penalties
are not the same. Typically this means that opening a gap is more expensive
that extending a gap, but this is configurable.

Second, we will allow for discounts for homopolymer indels. If a gap is found
in a homopolymer (HP) region, we will (optionally) discount the penalty based
upon how long the HP is.
*/
func NewLocalAligner(opts *alignmentOptions) *pairwise {
	return &pairwise{opts: opts, isLocal: true}
}

/*
This is a global aligner, which means that the best alignment must include
the end of both the query and target. This is implemented by changing the
initialization and backtracking conditions.
*/
func NewGlobalAligner(opts *alignmentOptions) *pairwise {
	opts.ClippingDisable() // global alignment doesn't make sense with clipping
	return &pairwise{opts: opts, isLocal: false}
}

func rowColToIdx(row int, col int, colLen int) int {
	return row*colLen + col
}

func (sw *pairwise) Align(query string, target string) *PairwiseAlignment {
	if sw.opts.verbose {
		fmt.Println("Query: " + query)
		fmt.Println("Target: " + target)
	}
	qLen := len(query)
	tLen := len(target)
	var qRun []int
	var tRun []int
	if sw.opts.hpOpenScale > 0 || sw.opts.hpExtendScale > 0 {
		qRun = sequtils.HomopolymerRunLen(query)
		tRun = sequtils.HomopolymerRunLen(target)
	}

	rows := qLen + 1
	cols := tLen + 1

	// store cells in 'data', the row/col index is calculated
	// with the rowColToIdx function
	data := make([]swCell, rows*cols)

	negInf := float32(math.Inf(-1))

	if sw.opts.verbose {
		fmt.Printf("Setting up matrix (rows x cols => %d x %d)\n", rows, cols)
	}

	for i := range rows {
		idx := rowColToIdx(i, 0, cols)
		data[idx] = swCell{}
		data[idx].scoreM = 0
		data[idx].scoreI = negInf
		data[idx].scoreD = negInf
		data[idx].traceM = swTraceMatch
		data[idx].traceI = swTraceIns
		data[idx].traceD = swTraceDel
	}

	for j := range cols {
		idx := rowColToIdx(0, j, cols)
		data[idx] = swCell{}
		data[idx].scoreM = 0
		data[idx].scoreI = negInf
		data[idx].scoreD = negInf
		data[idx].traceM = swTraceMatch
		data[idx].traceI = swTraceIns
		data[idx].traceD = swTraceDel
	}

	bestScore := float32(0)
	bestRow, bestCol := 0, 0
	bestTrace := swTraceMatch

	// this is a simple dynamic programming loop that
	// fills out the entire data matrix.

	for i := 1; i < rows; i++ {
		// the baseline value is normally 0, but if we are clipping
		// from the left, then the baseline value should be the clipping
		// penalty.
		var leftClipBaseline float32 = 0
		var rightClipPenalty float32 = 0

		// what would be cigar penaltiy score be for the right side?
		if sw.opts.clippingOpenPenalty > 0 || sw.opts.clippingExtendPenalty > 0 {
			leftClipBaseline = -gapPenalty(i-1, sw.opts.clippingOpenPenalty, sw.opts.clippingExtendPenalty)
			rightClipPenalty = gapPenalty(rows-i-1, sw.opts.clippingOpenPenalty, sw.opts.clippingExtendPenalty)
		}

		for j := 1; j < cols; j++ {
			idx := rowColToIdx(i, j, cols)
			diag := rowColToIdx(i-1, j-1, cols)
			up := rowColToIdx(i-1, j, cols)
			left := rowColToIdx(i, j-1, cols)

			data[idx] = swCell{}

			// calculate gap penalties and hp discounts
			gio := sw.opts.gapOpenPenaltyIns
			gie := sw.opts.gapExtendPenaltyIns
			gdo := sw.opts.gapOpenPenaltyDel
			gde := sw.opts.gapExtendPenaltyDel

			if sw.opts.hpOpenScale > 0 || sw.opts.hpExtendScale > 0 {
				// choices here
				// * we could do the qRun for insertions, tRun for dels
				// * we could take the max
				//   runLen := max(qRun[i-1], tRun[j-1])
				// * we could take the average
				//   runLen := (qRun[i-1] + tRun[j-1])/2
				// * we could take a weighted average (a*tRun + (1-a)*qRun)
				//   runLen := ((1-a) * qRun[i-1] + (a * tRun[j-1])
				//
				// the most sensible option for *all cases* I think it the first.
				// remember the i,j coord is actually 1-based indexing

				if qRun[i-1] > 0 {
					// only process HP discounts if HP run is positive (i.e. not an N or other non-standard base)
					gio -= hpDiscount(qRun[i-1], sw.opts.hpOpenScale, sw.opts.hpOpenCap)
					gie -= hpDiscount(qRun[i-1], sw.opts.hpExtendScale, sw.opts.hpExtendCap)
				}
				if tRun[j-1] > 0 {
					// only process HP discounts if HP run is positive (i.e. not an N or other non-standard base)
					gdo -= hpDiscount(tRun[j-1], sw.opts.hpOpenScale, sw.opts.hpOpenCap)
					gde -= hpDiscount(tRun[j-1], sw.opts.hpExtendScale, sw.opts.hpExtendCap)
				}
			}

			// query and target are zero-based. the matrix has an extra row/col at start
			s := sw.opts.scoringMatrix.ScorePair(query[i-1], target[j-1])

			Mm := data[diag].scoreM + s
			Mi := data[diag].scoreI + s
			Md := data[diag].scoreD + s

			if Mm >= Mi && Mm >= Md {
				data[idx].scoreM = Mm
				data[idx].traceM = swTraceMatch
			} else if Mi >= Md {
				data[idx].scoreM = Mi
				data[idx].traceM = swTraceIns
			} else {
				data[idx].scoreM = Md
				data[idx].traceM = swTraceDel
			}

			// INS is from UP
			Im := data[up].scoreM - gio
			Ii := data[up].scoreI - gie

			if Im >= Ii {
				data[idx].scoreI = Im
				data[idx].traceI = swTraceMatch
			} else {
				data[idx].scoreI = Ii
				data[idx].traceI = swTraceIns
			}

			// DEL is from LEFT
			Dm := data[left].scoreM - gdo
			Dd := data[left].scoreD - gde

			if Dm >= Dd {
				data[idx].scoreD = Dm
				data[idx].traceD = swTraceMatch
			} else {
				data[idx].scoreD = Dd
				data[idx].traceD = swTraceDel
			}

			if sw.isLocal {
				// The baseline worst score is whatever the left clipping base is
				if data[idx].scoreM < leftClipBaseline {
					data[idx].scoreM = leftClipBaseline
				}
				if data[idx].scoreI < leftClipBaseline {
					data[idx].scoreI = leftClipBaseline
				}
				if data[idx].scoreD < leftClipBaseline {
					data[idx].scoreD = leftClipBaseline
				}
				// track global best endpoint + state (taking into account right clipping)
				if data[idx].scoreM-rightClipPenalty >= bestScore {
					bestScore, bestRow, bestCol, bestTrace = data[idx].scoreM-rightClipPenalty, i, j, swTraceMatch
				}
				if data[idx].scoreI-rightClipPenalty >= bestScore {
					bestScore, bestRow, bestCol, bestTrace = data[idx].scoreI-rightClipPenalty, i, j, swTraceIns
				}
				if data[idx].scoreD-rightClipPenalty >= bestScore {
					bestScore, bestRow, bestCol, bestTrace = data[idx].scoreD-rightClipPenalty, i, j, swTraceDel
				}
			}
		}
	}

	if !sw.isLocal {
		// for global, we want to start backtracking from the end of both sequences
		bestRow = rows - 1
		bestCol = cols - 1

		// determine which state is best at the end
		idx := rowColToIdx(bestRow, bestCol, cols)
		if data[idx].scoreM >= data[idx].scoreI && data[idx].scoreM >= data[idx].scoreD {
			bestScore = data[idx].scoreM
			bestTrace = swTraceMatch
		} else if data[idx].scoreI >= data[idx].scoreD {
			bestScore = data[idx].scoreI
			bestTrace = swTraceIns
		} else {
			bestScore = data[idx].scoreD
			bestTrace = swTraceDel
		}
	}

	// Now we can backtrack to get the alignment
	i, j := bestRow, bestCol
	curTrace := bestTrace
	cigar := ""
	var curScore float32
	switch curTrace {
	case swTraceMatch:
		curScore = data[rowColToIdx(i, j, cols)].scoreM
	case swTraceIns:
		curScore = data[rowColToIdx(i, j, cols)].scoreI
	case swTraceDel:
		curScore = data[rowColToIdx(i, j, cols)].scoreD
	}

	if sw.opts.verbose {
		for j := range cols {
			if j > 0 {
				fmt.Printf("%6s ", target[j-1:j])
			} else {
				fmt.Printf("      ")
			}
		}
		fmt.Println()
		for i := range rows {
			if i > 0 {
				fmt.Printf("%s ", query[i-1:i])
			} else {
				fmt.Printf("  ")
			}
			for j := range cols {
				fmt.Printf("%*.1f ", 6, data[rowColToIdx(i, j, cols)].scoreM)
			}
			fmt.Println()
		}

		fmt.Println("Backtracking")
	}
	limit := max(qLen * tLen)
	// for the backtrack, if we hit a cell that drops below zero,
	// we know we are in a left-clipping region
	for ((sw.isLocal && curScore > 0) || (!sw.isLocal && (i > 0 || j > 0))) && limit > 0 {
		limit--
		if sw.opts.verbose {
			fmt.Printf("(%d,%d) %v, %f %s\n", i, j, curTrace, curScore, cigar)
		}

		idx := rowColToIdx(i, j, cols)
		var nextTrace swCellTrace

		switch curTrace {
		case swTraceMatch:
			cigar = "M" + cigar
			nextTrace = data[idx].traceM
			i--
			j--
		case swTraceIns:
			cigar = "I" + cigar
			nextTrace = data[idx].traceI
			i--
		case swTraceDel:
			cigar = "D" + cigar
			nextTrace = data[idx].traceD
			j--
		}

		// This only happens when in global mode, but we'll add an if gate anyway.
		if !sw.isLocal {
			if i < 0 {
				i = 0
				cigar = "D" + cigar[1:]
			}
			if j < 0 {
				j = 0
				cigar = "I" + cigar[1:]
			}
		}

		curTrace = nextTrace
		idx2 := rowColToIdx(i, j, cols)

		// get the score for the next (now current) cell
		switch curTrace {
		case swTraceMatch:
			curScore = data[idx2].scoreM
		case swTraceIns:
			curScore = data[idx2].scoreI
		case swTraceDel:
			curScore = data[idx2].scoreD
		case swTraceStop:
			if sw.opts.verbose {
				fmt.Println("Hit stop trace, ending backtrack")
			}
			curScore = 0
		}
		if sw.opts.verbose {
			fmt.Printf("   => next_score: %f, next_trace: %v\n", curScore, curTrace)
		}
	}

	queryStart := i // we've alreday decremented i, but i is also offset by 1, so it's a wash
	targetStart := j
	queryEnd := bestRow

	if sw.opts.clippingOpenPenalty > 0 || sw.opts.clippingExtendPenalty > 0 {
		// left clipping
		for ; i > 0; i-- {
			cigar = "S" + cigar
		}

		// right clipping
		for i := bestRow; i < rows-1; i++ {
			cigar += "S"
		}
		queryStart = 0
		queryEnd = len(query)
	}

	return &PairwiseAlignment{
		Query:       query,
		Target:      target,
		Score:       bestScore,
		CIGAR:       CigarCondense(cigar),
		QueryStart:  queryStart,
		QueryEnd:    queryEnd,
		TargetStart: targetStart,
		TargetEnd:   bestCol,
	}
}
