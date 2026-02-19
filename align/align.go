package align

import (
	"fmt"
	"math"
	"strconv"
	"strings"

	"github.com/compgen-io/cgltk/sequtils"
	"github.com/compgen-io/cgltk/utils"
)

type PairwiseAligner interface {
	Align(query string, target string) *PairwiseAlignment
}

type ScoringMatrix interface {
	ScorePair(one byte, two byte) float32
}

type PairwiseAlignment struct {
	Query        string
	Target       string
	QueryStart   int
	QueryEnd     int
	TargetStart  int
	TargetEnd    int
	Score        float32
	CIGAR        string
	QueryRevComp bool
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
		gapOpenPenaltyDel:     5,
		gapExtendPenaltyDel:   2,
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
		gapOpenPenaltyDel:     2,
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
	last := ""
	var ret strings.Builder
	for i := 0; i < len(s); i++ {
		cur := string(s[i])
		if len(last) == 0 || cur[0] == last[0] {
			last += cur
		} else {
			fmt.Fprintf(&ret, "%d%s", len(last), string(last[0]))
			last = cur
		}
	}
	if len(last) > 0 {
		fmt.Fprintf(&ret, "%d%s", len(last), string(last[0]))
	}
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
	qStr := ""
	tStr := ""
	alnStr := ""
	qPos := a.QueryStart
	tPos := a.TargetStart
	cigarExpanded, err := CigarExpand(a.CIGAR)
	if err != nil {
		return fmt.Sprintf("Error expanding CIGAR: %v", err)
	}
	for i := 0; i < len(cigarExpanded); i++ {
		// fmt.Printf("qStr: %s\ntStr: %s\n-\n", qStr, tStr)
		op := cigarExpanded[i]
		switch op {
		case 'M':
			qStr += string(a.Query[qPos])
			tStr += string(a.Target[tPos])
			if sequtils.DNAMatches(a.Query[qPos], a.Target[tPos]) {
				alnStr += "|"
			} else {
				alnStr += "."
			}
			qPos++
			tPos++
		case 'D':
			qStr += "-"
			tStr += string(a.Target[tPos])
			alnStr += " "
			tPos++
		case 'I':
			qStr += string(a.Query[qPos])
			alnStr += " "
			tStr += "-"
			qPos++
		case 'S':
			qStr += string(a.Query[qPos])
			tStr += " "
			alnStr += "-"
			qPos++
		}
	}
	var qName string
	if !a.QueryRevComp {
		qName = fmt.Sprintf("Query (%d-%d)", a.QueryStart+1, a.QueryEnd)
	} else {
		qName = fmt.Sprintf("Query (%d-%d)", a.QueryEnd, a.QueryStart+1)
	}
	tName := fmt.Sprintf("Target (%d-%d)", a.TargetStart+1, a.TargetEnd)

	maxNameLen := max(len(qName), len(tName))

	qName = fmt.Sprintf("%-*s", maxNameLen, qName)
	tName = fmt.Sprintf("%-*s", maxNameLen, tName)

	qStr = fmt.Sprintf("%s: %s", qName, qStr)
	tStr = fmt.Sprintf("%s: %s", tName, tStr)

	aName := fmt.Sprintf("%-*s", maxNameLen, " ")
	aStr := fmt.Sprintf("%s: %s", aName, alnStr)

	ret := fmt.Sprintf(`%s
%s
%s
CIGAR: %s
Score: %s`, qStr, aStr, tStr, a.CIGAR, utils.TrimFloat(float64(a.Score), 2))
	return ret
}

func (a *PairwiseAlignment) TargetAlignedStr() string {
	tStr := ""
	qPos := a.QueryStart
	tPos := a.TargetStart
	cigarExpanded, err := CigarExpand(a.CIGAR)
	if err != nil {
		return fmt.Sprintf("Error expanding CIGAR: %v", err)
	}
	for i := 0; i < len(cigarExpanded); i++ {
		// fmt.Printf("qStr: %s\ntStr: %s\n-\n", qStr, tStr)
		op := cigarExpanded[i]
		switch op {
		case 'M':
			tStr += string(a.Target[tPos])
			qPos++
			tPos++
		case 'D':
			tStr += string(a.Target[tPos])
			tPos++
		case 'I':
			tStr += "-"
			qPos++
		case 'S':
			tStr += " "
			qPos++
		}
	}
	return tStr
}

func (a *PairwiseAlignment) QueryAlignedStr() string {
	qStr := ""
	qPos := a.QueryStart
	tPos := a.TargetStart
	cigarExpanded, err := CigarExpand(a.CIGAR)
	if err != nil {
		return fmt.Sprintf("Error expanding CIGAR: %v", err)
	}
	for i := 0; i < len(cigarExpanded); i++ {
		// fmt.Printf("qStr: %s\ntStr: %s\n-\n", qStr, tStr)
		op := cigarExpanded[i]
		switch op {
		case 'M':
			qStr += string(a.Query[qPos])
			qPos++
			tPos++
		case 'D':
			qStr += "-"
			tPos++
		case 'I':
			qStr += string(a.Query[qPos])
			qPos++
		case 'S':
			qStr += string(a.Query[qPos])
			qPos++
		}
	}
	return qStr
}
