package align

import (
	"fmt"
	"math"
	"strconv"
	"strings"
	"sync"

	"github.com/compgen-io/cgltk/seqio"
	"github.com/compgen-io/cgltk/support/sequtils"
	"github.com/compgen-io/cgltk/support/utils"
)

type PairwiseAligner interface {
	Align(query seqio.SeqQual, target seqio.SeqQual) *PairwiseAlignment
}

type ScoringMatrix interface {
	ScorePair(one byte, two byte) float32
}

type PairwiseAlignment struct {
	Query        seqio.SeqQual
	Target       seqio.SeqQual
	QueryStart   int
	QueryEnd     int
	TargetStart  int
	TargetEnd    int
	Score        float32
	CIGAR        string
	cigarExpanded string
}

type matchMismatch struct {
	match    float32
	mismatch float32
}

type alignmentOptions struct {
	scoringMatrix         ScoringMatrix
	gapOpenPenaltyIns     float32
	gapExtendPenaltyIns   float32
	gapOpenPenaltyDel     float32
	gapExtendPenaltyDel   float32
	clippingOpenPenalty   float32
	clippingExtendPenalty float32
	hpOpenScale           float32
	hpExtendScale         float32
	hpOpenCap             float32
	hpExtendCap           float32
	verbose               bool
}

func DnaAlignmentDefaults() *alignmentOptions {
	// default short-read alignment scoring
	return &alignmentOptions{
		scoringMatrix:         MatchMismatchScoring(1, 2),
		gapOpenPenaltyIns:     6,
		gapExtendPenaltyIns:   1,
		gapOpenPenaltyDel:     6,
		gapExtendPenaltyDel:   1,
		clippingOpenPenalty:   5,
		clippingExtendPenalty: 1,

		// homopolymer errors aren't typical with illumina short reads
		hpOpenScale:   0,
		hpExtendScale: 0,

		hpOpenCap:   0,
		hpExtendCap: 0,
	}
}

func OntAlignmentDefaults() *alignmentOptions {
	// default short-read alignment scoring
	return &alignmentOptions{
		scoringMatrix:         MatchMismatchScoring(1, 1),
		gapOpenPenaltyIns:     2,
		gapExtendPenaltyIns:   1,
		gapOpenPenaltyDel:     3,
		gapExtendPenaltyDel:   1,
		clippingOpenPenalty:   5,
		clippingExtendPenalty: 1,

		// homopolymer errors are typical with oxford nanopore long reads
		hpOpenScale:   1,
		hpExtendScale: 0.4,

		hpOpenCap:   2,   // limit discount to at most make it a free indel (when hplen > 4, discount = gapOpenPenalty)
		hpExtendCap: 0.8, // going from 4->5 or 5->6 is cheap (0.2) -- not free
	}
}

// The scoring matrix (match/mismatch scoring) to use
func (a *alignmentOptions) ScoringMatrix(matrix ScoringMatrix) *alignmentOptions {
	a.scoringMatrix = matrix
	return a
}

// decay the gap extension length by the length of the homopolymer
// gap_penalty = (gap_open / hp_length) + n * (gap_extend / hp_length)
func (a *alignmentOptions) HomopolymerDiscount(openScale, openCap, extendScale, extendCap float32) *alignmentOptions {
	a.hpOpenScale = openScale
	a.hpOpenCap = openCap
	a.hpExtendScale = extendScale
	a.hpExtendCap = extendCap
	return a
}

// penalty for opening a gap (insertion)
// gap_penalty = gap_open + (n * gap_extend)
func (a *alignmentOptions) GapPenaltyIns(open, extend float32) *alignmentOptions {
	a.gapOpenPenaltyIns = float32(math.Abs(float64(open)))
	a.gapExtendPenaltyIns = float32(math.Abs(float64(extend)))
	return a
}

// penalty for opening a gap (deletions)
// gap_penalty = gap_open + (n * gap_extend)
func (a *alignmentOptions) GapPenaltyDel(open, extend float32) *alignmentOptions {
	a.gapOpenPenaltyDel = float32(math.Abs(float64(open)))
	a.gapExtendPenaltyDel = float32(math.Abs(float64(extend)))
	return a
}

// penalty for opening a 5' or 3' soft clipping gap
// gap_penalty = gap_open + (n * gap_extend)
func (a *alignmentOptions) ClippingPenalty(open, extend float32) *alignmentOptions {
	a.clippingOpenPenalty = float32(math.Abs(float64(open)))
	a.clippingExtendPenalty = float32(math.Abs(float64(extend)))
	return a
}

// penalty for opening a 5' or 3' soft clipping gap
// gap_penalty = gap_open + (n * gap_extend)
func (a *alignmentOptions) ClippingDisable() *alignmentOptions {
	a.clippingOpenPenalty = -1
	a.clippingExtendPenalty = -1
	return a
}

func (a *alignmentOptions) Verbose(verbose bool) *alignmentOptions {
	a.verbose = verbose
	return a
}

func MatchMismatchScoring(match int, mismatch int) *matchMismatch {
	return &matchMismatch{match: float32(match), mismatch: float32(math.Abs(float64(mismatch)))}
}

func (m *matchMismatch) ScorePair(one byte, two byte) float32 {
	if sequtils.DNAMatches(one, two) {
		return m.match
	}
	return -m.mismatch
}

// calculate a gap penalty for k bases with an open and extend penalty
// penalty = open + (k-1) * extend
func gapPenalty(k int, open, extend float32) float32 {
	if k <= 0 {
		return 0
	}
	if k == 1 {
		return open
	}
	return open + float32(k-1)*extend
}

// discounts to the gap penalties calculated for homopolymers
// hp discounts only occur if hpLen >= 2
// discount = min(cap, scale * log2(hpLen))
func hpDiscount(hpLen int, scale, cap float32) float32 {
	if hpLen < 2 {
		return 0
	}
	ret := scale * float32(math.Log2(float64(hpLen)))
	if ret > cap {
		ret = cap
	}
	return ret
}

// take an extended cigar string and convert
// it into a condensed string:
//
// IIMMMMMDMM => 2I5M1D2M
func CigarCondense(s string) string {
	if len(s) == 0 {
		return ""
	}
	var ret strings.Builder
	lastChar := s[0]
	count := 1
	for i := 1; i < len(s); i++ {
		if s[i] == lastChar {
			count++
		} else {
			fmt.Fprintf(&ret, "%d%c", count, lastChar)
			lastChar = s[i]
			count = 1
		}
	}
	fmt.Fprintf(&ret, "%d%c", count, lastChar)
	return ret.String()
}

// take an extended cigar string and convert
// it into a condensed string:
//
// 2I5M1D2M => IIMMMMMDMM
func CigarExpand(s string) (string, error) {
	countBuf := ""
	var ret strings.Builder
	for i := 0; i < len(s); i++ {
		if strings.ContainsAny(s[i:i+1], "0123456789") {
			countBuf += string(s[i])
		} else {
			op := string(s[i])
			if count, err := strconv.Atoi(countBuf); err != nil {
				return "", err
			} else {
				for range count {
					ret.WriteString(op)
				}
			}
			countBuf = ""
		}
	}
	if countBuf != "" {
		return "", fmt.Errorf("invalid CIGAR string")
	}
	return ret.String(), nil
}

func (a *PairwiseAlignment) String() string {
	// fmt.Printf("qPos: %d-%d, tPos: %d-%d\n", a.QueryStart, a.QueryEnd, a.TargetStart, a.TargetEnd)
	// fmt.Printf("CIGAR: %s\n", a.CIGAR)
	var qBuf, tBuf, alnBuf strings.Builder
	qPos := a.QueryStart
	tPos := a.TargetStart
	for i := 0; i < len(a.cigarExpanded); i++ {
		// fmt.Printf("qStr: %s\ntStr: %s\n-\n", qBuf.String(), tBuf.String())
		op := a.cigarExpanded[i]
		switch op {
		case 'M':
			qBuf.WriteByte(a.Query.Seq()[qPos])
			tBuf.WriteByte(a.Target.Seq()[tPos])
			if sequtils.DNAMatches(a.Query.Seq()[qPos], a.Target.Seq()[tPos]) {
				alnBuf.WriteByte('|')
			} else {
				alnBuf.WriteByte('.')
			}
			qPos++
			tPos++
		case 'D':
			qBuf.WriteByte('-')
			tBuf.WriteByte(a.Target.Seq()[tPos])
			alnBuf.WriteByte(' ')
			tPos++
		case 'I':
			qBuf.WriteByte(a.Query.Seq()[qPos])
			alnBuf.WriteByte(' ')
			tBuf.WriteByte('-')
			qPos++
		case 'S':
			qBuf.WriteByte(a.Query.Seq()[qPos])
			tBuf.WriteByte(' ')
			alnBuf.WriteByte('-')
			qPos++
		}
	}
	var qName, tName string
	if !a.Query.IsRevComp() {
		qName = fmt.Sprintf("%s (%d-%d)", a.Query.Name(), a.QueryStart+1, a.QueryEnd)
	} else {
		qName = fmt.Sprintf("%s (%d-%d)", a.Query.Name(), a.QueryEnd, a.QueryStart+1)
	}
	if !a.Target.IsRevComp() {
		tName = fmt.Sprintf("%s (%d-%d)", a.Target.Name(), a.TargetStart+1, a.TargetEnd)
	} else {
		tName = fmt.Sprintf("%s (%d-%d)", a.Target.Name(), a.TargetEnd, a.TargetStart+1)
	}

	maxNameLen := max(len(qName), len(tName))

	qName = fmt.Sprintf("%-*s", maxNameLen, qName)
	tName = fmt.Sprintf("%-*s", maxNameLen, tName)

	qStr := fmt.Sprintf("%s: %s", qName, qBuf.String())
	tStr := fmt.Sprintf("%s: %s", tName, tBuf.String())

	aName := fmt.Sprintf("%-*s", maxNameLen, " ")
	aStr := fmt.Sprintf("%s: %s", aName, alnBuf.String())

	ret := fmt.Sprintf(`%s
%s
%s
CIGAR: %s
Score: %s`, qStr, aStr, tStr, a.CIGAR, utils.TrimFloat(float64(a.Score), 2))
	return ret
}

func (a *PairwiseAlignment) Matches() int {
	qPos := a.QueryStart
	tPos := a.TargetStart
	matches := 0
	for i := 0; i < len(a.cigarExpanded); i++ {
		op := a.cigarExpanded[i]
		switch op {
		case 'M':
			if sequtils.DNAMatches(a.Query.Seq()[qPos], a.Target.Seq()[tPos]) {
				matches++
			}
			qPos++
			tPos++
		case 'D':
			tPos++
		case 'I':
			qPos++
		case 'S':
			qPos++
		}
	}
	return matches
}

func (a *PairwiseAlignment) TargetAlignedStr() string {
	var tBuf strings.Builder
	qPos := a.QueryStart
	tPos := a.TargetStart
	for i := 0; i < len(a.cigarExpanded); i++ {
		op := a.cigarExpanded[i]
		switch op {
		case 'M':
			tBuf.WriteByte(a.Target.Seq()[tPos])
			qPos++
			tPos++
		case 'D':
			tBuf.WriteByte(a.Target.Seq()[tPos])
			tPos++
		case 'I':
			tBuf.WriteByte('-')
			qPos++
		case 'S':
			tBuf.WriteByte(' ')
			qPos++
		}
	}
	return tBuf.String()
}

// Return the target string, relative to the plus strand
func (a *PairwiseAlignment) TargetStrPlus() string {
	if a.Target.IsRevComp() {
		return sequtils.ReverseCompliment(a.Target.Seq()[a.TargetStart:a.TargetEnd])
	}
	return a.Target.Seq()[a.TargetStart:a.TargetEnd]
}

func (a *PairwiseAlignment) TargetSub() seqio.SeqQual {
	return a.Target.Sub(a.TargetStart, a.TargetEnd)
}

// Return the target string
func (a *PairwiseAlignment) TargetStr() string {
	return a.Target.Seq()[a.TargetStart:a.TargetEnd]
}

func (a *PairwiseAlignment) QueryAlignedStr() string {
	var qBuf strings.Builder
	qPos := a.QueryStart
	tPos := a.TargetStart
	for i := 0; i < len(a.cigarExpanded); i++ {
		op := a.cigarExpanded[i]
		switch op {
		case 'M':
			qBuf.WriteByte(a.Query.Seq()[qPos])
			qPos++
			tPos++
		case 'D':
			qBuf.WriteByte('-')
			tPos++
		case 'I':
			qBuf.WriteByte(a.Query.Seq()[qPos])
			qPos++
		case 'S':
			qBuf.WriteByte(a.Query.Seq()[qPos])
			qPos++
		}
	}
	return qBuf.String()
}

// Return the query string, relative to the plus strand
func (a *PairwiseAlignment) QueryStrPlus() string {
	if a.Query.IsRevComp() {
		return sequtils.ReverseCompliment(a.Query.Seq()[a.QueryStart:a.QueryEnd])
	}
	return a.Query.Seq()[a.QueryStart:a.QueryEnd]
}

func (a *PairwiseAlignment) QuerySub() seqio.SeqQual {
	return a.Query.Sub(a.QueryStart, a.QueryEnd)
}

func (a *PairwiseAlignment) QueryStr() string {
	return a.Query.Seq()[a.QueryStart:a.QueryEnd]
}

type PairwiseAlignmentPromise struct {
	results []*PairwiseAlignment
	wg      sync.WaitGroup
}

func newPairwiseAlignmentPromise(resultCount int) *PairwiseAlignmentPromise {
	return &PairwiseAlignmentPromise{
		results: make([]*PairwiseAlignment, resultCount),
		wg:      sync.WaitGroup{},
	}
}

func (pap *PairwiseAlignmentPromise) add(delta int) {
	pap.wg.Add(delta)
}

func (pap *PairwiseAlignmentPromise) done() {
	pap.wg.Done()
}

func (pap *PairwiseAlignmentPromise) setResult(idx int, aln *PairwiseAlignment) {
	pap.results[idx] = aln
}

func (pap *PairwiseAlignmentPromise) Result() *PairwiseAlignment {
	pap.wg.Wait()

	bestAln := pap.results[0]
	// Find the best alignment
	for i, result := range pap.results {
		if i > 0 {
			if result.Score > bestAln.Score {
				bestAln = result
			}
		}
	}
	// fmt.Println(bestAln.String())
	return bestAln
}

func AlignBatch(aligner PairwiseAligner, sem utils.Semaphore, queries []seqio.SeqQual, targets []seqio.SeqQual) *PairwiseAlignmentPromise {
	// We will be doing all of the calls in parallel.
	// The semaphore will keep track of the number of concurrent jobs and
	// will cap it at maxWorkers

	promise := newPairwiseAlignmentPromise(len(queries) * len(targets))

	for i, query := range queries {
		query := query // capture loop var
		i := i         // capture loop var
		for j, target := range targets {
			target := target // capture loop var
			j := j           // capture loop var
			promise.add(1)

			sem.Acquire() // acquire
			go func() {
				defer promise.done()
				defer sem.Release() // release
				// run the alignment
				promise.setResult(i*len(targets)+j, aligner.Align(query, target))
			}()
		}
	}

	return promise
}
