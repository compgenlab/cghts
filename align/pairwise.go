package align

import (
	"fmt"
	"math"

	"github.com/compgen-io/cgltk/seqio"
	"github.com/compgen-io/cgltk/support/sequtils"
)

type pairwise struct {
	opts              *alignmentOptions
	isLocal           bool
	useHPDiscount     bool
	scoreTable        [256][256]float32
	hpDiscountOpen    []float32
	hpDiscountExt     []float32
	hpDiscountOpenMax int
	hpDiscountExtMax  int
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
	ret := &pairwise{opts: opts, isLocal: true}
	ret.precalc()
	return ret
}

/*
This is a global aligner, which means that the best alignment must include
the end of both the query and target. This is implemented by changing the
initialization and backtracking conditions.
*/
func NewGlobalAligner(opts *alignmentOptions) *pairwise {
	opts.ClippingDisable() // global alignment doesn't make sense with clipping
	ret := &pairwise{opts: opts, isLocal: false}
	ret.precalc()
	return ret
}

func (sw *pairwise) precalc() {
	sw.hpDiscountOpen = make([]float32, 0)
	sw.hpDiscountExt = make([]float32, 0)

	i := 1
	for {
		tmp := hpDiscount(i, sw.opts.hpOpenScale, sw.opts.hpOpenCap)
		sw.hpDiscountOpen = append(sw.hpDiscountOpen, tmp)
		if tmp >= sw.opts.hpOpenCap {
			break
		}
		i++
	}
	sw.hpDiscountOpenMax = i
	i = 1
	for {
		tmp := hpDiscount(i, sw.opts.hpExtendScale, sw.opts.hpExtendCap)
		sw.hpDiscountExt = append(sw.hpDiscountExt, tmp)
		if tmp >= sw.opts.hpExtendCap {
			break
		}
		i++
	}
	sw.hpDiscountExtMax = i
	sw.useHPDiscount = sw.opts.hpOpenScale > 0 || sw.opts.hpExtendScale > 0
	for i := range 256 {
		for j := range 256 {
			sw.scoreTable[i][j] = sw.opts.scoringMatrix.ScorePair(byte(i), byte(j))
		}
	}
}
func rowColToIdx(row int, col int, colLen int) int {
	return row*colLen + col
}

func (sw *pairwise) Align(query seqio.SeqQual, target seqio.SeqQual) *PairwiseAlignment {
	queryStr := query.Seq()
	targetStr := target.Seq()
	if sw.opts.verbose {
		fmt.Println("Query: " + queryStr)
		fmt.Println("Target: " + targetStr)
	}
	qLen := len(queryStr)
	tLen := len(targetStr)
	var qRun []int
	var tRun []int
	if sw.opts.hpOpenScale > 0 || sw.opts.hpExtendScale > 0 {
		qRun = sequtils.HomopolymerRunLen(queryStr)
		tRun = sequtils.HomopolymerRunLen(targetStr)
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

		rowStart := i * cols
		for j := 1; j < cols; j++ {
			idx := rowStart + j
			up := idx - cols
			diag := up - 1
			left := idx - 1

			// calculate gap penalties and hp discounts
			gio := sw.opts.gapOpenPenaltyIns
			gie := sw.opts.gapExtendPenaltyIns
			gdo := sw.opts.gapOpenPenaltyDel
			gde := sw.opts.gapExtendPenaltyDel

			if sw.useHPDiscount {
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
					if qRun[i-1] <= sw.hpDiscountOpenMax {
						gio -= sw.hpDiscountOpen[qRun[i-1]-1]
					} else {
						gio -= sw.opts.hpOpenCap
					}

					if qRun[i-1] <= sw.hpDiscountExtMax {
						gie -= sw.hpDiscountExt[qRun[i-1]-1]
					} else {
						gie -= sw.opts.hpExtendCap
					}

					// gio -= hpDiscount(qRun[i-1], sw.opts.hpOpenScale, sw.opts.hpOpenCap)
					// gie -= hpDiscount(qRun[i-1], sw.opts.hpExtendScale, sw.opts.hpExtendCap)
				}
				if tRun[j-1] > 0 {
					if tRun[j-1] <= sw.hpDiscountOpenMax {
						gdo -= sw.hpDiscountOpen[tRun[j-1]-1]
					} else {
						gdo -= sw.opts.hpOpenCap
					}

					if tRun[j-1] <= sw.hpDiscountExtMax {
						gde -= sw.hpDiscountExt[tRun[j-1]-1]
					} else {
						gde -= sw.opts.hpExtendCap
					}

					// only process HP discounts if HP run is positive (i.e. not an N or other non-standard base)
					// gdo -= hpDiscount(tRun[j-1], sw.opts.hpOpenScale, sw.opts.hpOpenCap)
					// gde -= hpDiscount(tRun[j-1], sw.opts.hpExtendScale, sw.opts.hpExtendCap)
				}
			}

			// query and target are zero-based. the matrix has an extra row/col at start
			s := sw.scoreTable[queryStr[i-1]][targetStr[j-1]]

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
				if data[idx].scoreI-rightClipPenalty > bestScore {
					bestScore, bestRow, bestCol, bestTrace = data[idx].scoreI-rightClipPenalty, i, j, swTraceIns
				}
				if data[idx].scoreD-rightClipPenalty > bestScore {
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
	cigarBuf := make([]byte, 0, qLen+tLen)
	var curScore float32
	idx := i*cols + j
	switch curTrace {
	case swTraceMatch:
		curScore = data[idx].scoreM
	case swTraceIns:
		curScore = data[idx].scoreI
	case swTraceDel:
		curScore = data[idx].scoreD
	}

	if sw.opts.verbose {
		for j := range cols {
			if j > 0 {
				fmt.Printf("%6s ", targetStr[j-1:j])
			} else {
				fmt.Printf("      ")
			}
		}
		fmt.Println()
		for i := range rows {
			if i > 0 {
				fmt.Printf("%s ", queryStr[i-1:i])
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

	// for the backtrack, if we are a local alignment and we hit a cell that drops below zero,
	// we know we are in a left-clipping region
	//
	// for global alignments, stop when we hit the upper left corner

	// failsafe to make sure we don't get stuck.
	limit := max(qLen * tLen)

	for ((sw.isLocal && curScore > 0) || (!sw.isLocal && (i > 0 || j > 0))) && limit > 0 {
		limit--
		if sw.opts.verbose {
			fmt.Printf("(%d,%d) %v, %f %s\n", i, j, curTrace, curScore, string(cigarBuf))
		}

		var nextTrace swCellTrace

		switch curTrace {
		case swTraceMatch:
			cigarBuf = append(cigarBuf, 'M')
			nextTrace = data[idx].traceM
			i--
			j--
			idx -= cols + 1
		case swTraceIns:
			cigarBuf = append(cigarBuf, 'I')
			nextTrace = data[idx].traceI
			i--
			idx -= cols
		case swTraceDel:
			cigarBuf = append(cigarBuf, 'D')
			nextTrace = data[idx].traceD
			j--
			idx--
		}

		// This only happens when in global mode, but we'll add an if gate anyway.
		if !sw.isLocal {
			if i < 0 {
				i = 0
				idx += cols
				cigarBuf[len(cigarBuf)-1] = 'D'
			}
			if j < 0 {
				j = 0
				idx++
				cigarBuf[len(cigarBuf)-1] = 'I'
			}
		}

		curTrace = nextTrace

		// get the score for the next (now current) cell
		switch curTrace {
		case swTraceMatch:
			curScore = data[idx].scoreM
		case swTraceIns:
			curScore = data[idx].scoreI
		case swTraceDel:
			curScore = data[idx].scoreD
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

	queryStart := i // we've already decremented i, but i is also offset by 1, so it's a wash
	targetStart := j
	queryEnd := bestRow

	if sw.opts.clippingOpenPenalty > 0 || sw.opts.clippingExtendPenalty > 0 {
		// left clipping: append before reversing so they end up at the front
		for ; i > 0; i-- {
			cigarBuf = append(cigarBuf, 'S')
		}
	}

	// reverse the buffer to get forward-order CIGAR
	for l, r := 0, len(cigarBuf)-1; l < r; l, r = l+1, r-1 {
		cigarBuf[l], cigarBuf[r] = cigarBuf[r], cigarBuf[l]
	}

	if sw.opts.clippingOpenPenalty > 0 || sw.opts.clippingExtendPenalty > 0 {
		// right clipping: append after reversing so they end up at the back
		for i := bestRow; i < rows-1; i++ {
			cigarBuf = append(cigarBuf, 'S')
		}
		queryStart = 0
		queryEnd = len(queryStr)
	}

	cigar := string(cigarBuf)

	return &PairwiseAlignment{
		Query:         query,
		Target:        target,
		QueryStart:    queryStart,
		QueryEnd:      queryEnd,
		TargetStart:   targetStart,
		TargetEnd:     bestCol,
		Score:         bestScore,
		CIGAR:         CigarCondense(cigar),
		cigarExpanded: cigar,
	}
}
