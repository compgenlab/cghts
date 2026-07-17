package align

import (
	"fmt"
	"math"
	"strconv"
	"strings"
	"sync"

	"github.com/compgenlab/cghts/seqio"
	"github.com/compgenlab/cghts/support/sequtils"
	"github.com/compgenlab/cghts/support/utils"
)

// PairwiseAligner aligns a single query sequence against a single target
// sequence and returns the resulting [PairwiseAlignment]. It is implemented by
// the aligners returned from [NewLocalAligner], [NewGlobalAligner], and
// [NewSemiGlobalAligner].
type PairwiseAligner interface {
	Align(query seqio.SeqQual, target seqio.SeqQual) *PairwiseAlignment
}

// ScoringMatrix scores the alignment of one query base against one target base.
// ScorePair returns the (signed) score for aligning bytes one and two; positive
// values reward matches and negative values penalize mismatches.
type ScoringMatrix interface {
	ScorePair(one byte, two byte) float32
}

// PairwiseAlignment is the result of aligning a query against a target.
//
// QueryStart/QueryEnd and TargetStart/TargetEnd are 0-based, half-open
// coordinates (start inclusive, end exclusive) into the respective sequences,
// delimiting the aligned region. Score is the alignment score under the
// configured scoring matrix and gap penalties. CIGAR is the run-length-encoded
// alignment (for example "5M1D2M") using the operations M (match/mismatch),
// I (insertion in the query relative to the target), D (deletion relative to
// the target), and S (soft clip).
type PairwiseAlignment struct {
	Query         seqio.SeqQual
	Target        seqio.SeqQual
	QueryStart    int
	QueryEnd      int
	TargetStart   int
	TargetEnd     int
	Score         float32
	CIGAR         string
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

// DnaAlignmentDefaults returns alignment options tuned for Illumina short
// reads: match +1, mismatch -2, affine gap penalties of open 6 / extend 1 for
// both insertions and deletions, and soft-clipping penalties of open 5 /
// extend 1. Homopolymer discounts are disabled, since homopolymer-length errors
// are uncommon with short reads. The returned value can be customized via its
// fluent setter methods.
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

// OntAlignmentDefaults returns alignment options tuned for Oxford Nanopore
// long reads: match +1, mismatch -1, looser affine gap penalties (insertions
// open 2 / extend 1; deletions open 3 / extend 1), soft-clipping penalties of
// open 5 / extend 1, and enabled homopolymer discounts. The homopolymer
// discounts reduce the cost of indels inside homopolymer runs (scaled by the
// log of the run length and capped) so that a long run makes a single indel
// nearly free. The returned value can be customized via its fluent setter
// methods.
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

// ScoringMatrix sets the match/mismatch scoring matrix to use and returns the
// options for chaining.
func (a *alignmentOptions) ScoringMatrix(matrix ScoringMatrix) *alignmentOptions {
	a.scoringMatrix = matrix
	return a
}

// HomopolymerDiscount configures the discount applied to gap open and extend
// penalties for indels that fall inside a homopolymer run. The discount for a
// run of length L is min(cap, scale*log2(L)) and is subtracted from the gap
// penalty; openScale/openCap apply to gap opening and extendScale/extendCap to
// gap extension. Setting all values to 0 disables the discount. Returns the
// options for chaining.
func (a *alignmentOptions) HomopolymerDiscount(openScale, openCap, extendScale, extendCap float32) *alignmentOptions {
	a.hpOpenScale = openScale
	a.hpOpenCap = openCap
	a.hpExtendScale = extendScale
	a.hpExtendCap = extendCap
	return a
}

// GapPenaltyIns sets the affine gap penalties for insertions (query bases not
// present in the target). The penalty for a gap of length n is
// open + (n-1)*extend. The absolute values of open and extend are used (the
// penalties are subtracted internally). Returns the options for chaining.
func (a *alignmentOptions) GapPenaltyIns(open, extend float32) *alignmentOptions {
	a.gapOpenPenaltyIns = float32(math.Abs(float64(open)))
	a.gapExtendPenaltyIns = float32(math.Abs(float64(extend)))
	return a
}

// GapPenaltyDel sets the affine gap penalties for deletions (target bases not
// present in the query). The penalty for a gap of length n is
// open + (n-1)*extend. The absolute values of open and extend are used (the
// penalties are subtracted internally). Returns the options for chaining.
func (a *alignmentOptions) GapPenaltyDel(open, extend float32) *alignmentOptions {
	a.gapOpenPenaltyDel = float32(math.Abs(float64(open)))
	a.gapExtendPenaltyDel = float32(math.Abs(float64(extend)))
	return a
}

// ClippingPenalty sets the affine penalties for soft-clipping query bases at
// the 5' or 3' end (local alignment only). The penalty for clipping n bases is
// open + (n-1)*extend. The absolute values of open and extend are used. Returns
// the options for chaining.
func (a *alignmentOptions) ClippingPenalty(open, extend float32) *alignmentOptions {
	a.clippingOpenPenalty = float32(math.Abs(float64(open)))
	a.clippingExtendPenalty = float32(math.Abs(float64(extend)))
	return a
}

// ClippingDisable turns off soft clipping by setting the clipping penalties to
// a negative sentinel, forcing the query to be aligned end to end. It is
// applied automatically by [NewGlobalAligner] and [NewSemiGlobalAligner].
// Returns the options for chaining.
func (a *alignmentOptions) ClippingDisable() *alignmentOptions {
	a.clippingOpenPenalty = -1
	a.clippingExtendPenalty = -1
	return a
}

// Verbose enables or disables debug output (the DP matrix and backtrack trace)
// printed during alignment. Returns the options for chaining.
func (a *alignmentOptions) Verbose(verbose bool) *alignmentOptions {
	a.verbose = verbose
	return a
}

// MatchMismatchScoring returns a [ScoringMatrix] that awards match for matching
// bases and -|mismatch| for mismatches. Base comparison honors IUPAC ambiguity
// codes, so for example N matches any base.
func MatchMismatchScoring(match int, mismatch int) *matchMismatch {
	return &matchMismatch{match: float32(match), mismatch: float32(math.Abs(float64(mismatch)))}
}

// ScorePair returns the match reward when one and two match (per IUPAC
// matching) and the negated mismatch penalty otherwise.
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

// CigarCondense converts an expanded (per-base) CIGAR string into run-length
// encoded form, for example:
//
//	IIMMMMMDMM => 2I5M1D2M
//
// An empty input returns an empty string.
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

// CigarExpand converts a run-length encoded CIGAR string into expanded
// (per-base) form, the inverse of [CigarCondense], for example:
//
//	2I5M1D2M => IIMMMMMDMM
//
// It returns an error if the string is malformed (for example a trailing count
// with no operation, or a non-numeric count).
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

// String returns a human-readable, multi-line rendering of the alignment: the
// gapped query and target rows with a middle line marking matches ('|') and
// mismatches ('.'), followed by the CIGAR string and score. Coordinates in the
// row labels are 1-based and oriented to the strand of each sequence.
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

// Matches returns the number of aligned positions where the query and target
// bases match (per IUPAC matching). Insertions, deletions, and soft clips do
// not count.
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

// TargetAlignedStr returns the target sequence over the aligned region with
// gap characters ('-') inserted for insertions and a space for soft-clipped
// columns, so it lines up column-for-column with [PairwiseAlignment.QueryAlignedStr].
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

// TargetStrPlus returns the aligned target substring oriented to the plus
// strand: if the target is reverse-complemented, the reverse complement of the
// aligned region is returned.
func (a *PairwiseAlignment) TargetStrPlus() string {
	if a.Target.IsRevComp() {
		return sequtils.ReverseComplement(a.Target.Seq()[a.TargetStart:a.TargetEnd])
	}
	return a.Target.Seq()[a.TargetStart:a.TargetEnd]
}

// TargetSub returns the aligned region of the target as a [seqio.SeqQual]
// (the subsequence from TargetStart to TargetEnd).
func (a *PairwiseAlignment) TargetSub() seqio.SeqQual {
	return a.Target.Sub(a.TargetStart, a.TargetEnd)
}

// TargetStr returns the aligned target substring (TargetStart to TargetEnd) as
// stored, without reorienting to the plus strand.
func (a *PairwiseAlignment) TargetStr() string {
	return a.Target.Seq()[a.TargetStart:a.TargetEnd]
}

// QueryAlignedStr returns the query sequence over the aligned region with gap
// characters ('-') inserted for deletions, so it lines up column-for-column
// with [PairwiseAlignment.TargetAlignedStr]. Soft-clipped query bases are
// included.
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

// QueryStrPlus returns the aligned query substring oriented to the plus strand:
// if the query is reverse-complemented, the reverse complement of the aligned
// region is returned.
func (a *PairwiseAlignment) QueryStrPlus() string {
	if a.Query.IsRevComp() {
		return sequtils.ReverseComplement(a.Query.Seq()[a.QueryStart:a.QueryEnd])
	}
	return a.Query.Seq()[a.QueryStart:a.QueryEnd]
}

// QuerySub returns the aligned region of the query as a [seqio.SeqQual] (the
// subsequence from QueryStart to QueryEnd).
func (a *PairwiseAlignment) QuerySub() seqio.SeqQual {
	return a.Query.Sub(a.QueryStart, a.QueryEnd)
}

// QueryStr returns the aligned query substring (QueryStart to QueryEnd) as
// stored, without reorienting to the plus strand.
func (a *PairwiseAlignment) QueryStr() string {
	return a.Query.Seq()[a.QueryStart:a.QueryEnd]
}

// PairwiseAlignmentPromise holds the pending results of an [AlignBatch] call.
// Call its Result method to block until every alignment completes and obtain
// the single best-scoring alignment.
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

// Result blocks until all batched alignments have finished and returns the
// alignment with the highest score. It returns nil if the batch produced no
// alignments.
func (pap *PairwiseAlignmentPromise) Result() *PairwiseAlignment {
	pap.wg.Wait()

	if len(pap.results) == 0 {
		return nil
	}
	bestAln := pap.results[0]
	// Find the best alignment
	for _, result := range pap.results[1:] {
		if result.Score > bestAln.Score {
			bestAln = result
		}
	}
	return bestAln
}

// AlignBatch aligns every query against every target concurrently, using
// aligner for each pairwise alignment. The sem semaphore bounds the number of
// alignments running at once. AlignBatch returns immediately with a
// [PairwiseAlignmentPromise]; call its Result method to wait for completion and
// retrieve the highest-scoring alignment across all query/target pairs.
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
