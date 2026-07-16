package htsio

import (
	"fmt"
	"strconv"
)

// CigarOp is a single run-length-encoded CIGAR operation: an operation
// character from MIDNSHP=X and the number of bases it applies to.
type CigarOp struct {
	Len int
	Op  byte
}

// ConsumesRef reports whether the operation consumes reference bases
// (M, D, N, =, X).
func (o CigarOp) ConsumesRef() bool {
	switch o.Op {
	case 'M', 'D', 'N', '=', 'X':
		return true
	}
	return false
}

// ConsumesQuery reports whether the operation consumes query (SEQ) bases
// (M, I, S, =, X). Hard clips (H) do not: their bases are absent from SEQ.
func (o CigarOp) ConsumesQuery() bool {
	switch o.Op {
	case 'M', 'I', 'S', '=', 'X':
		return true
	}
	return false
}

// String returns the operation in SAM text form, e.g. "10M".
func (o CigarOp) String() string {
	return strconv.Itoa(o.Len) + string(o.Op)
}

func isCigarOp(c byte) bool {
	switch c {
	case 'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X':
		return true
	}
	return false
}

// ParseCigar parses a SAM CIGAR string into its operations. It returns nil, nil
// for an empty or unspecified ("*") CIGAR.
//
// An error is returned for a length with no operation, an operation with no
// length, or an operation character outside MIDNSHP=X. 'B' is rejected: it is
// not a SAM operation, though a malformed BAM can decode to it.
//
// Unlike CigarRefLen and CigarQueryLen, which silently skip anything they do not
// recognize, ParseCigar rejects malformed input, so callers that must round-trip
// a CIGAR can detect corruption.
func ParseCigar(cigar string) ([]CigarOp, error) {
	if cigar == "" || cigar == "*" {
		return nil, nil
	}
	var ops []CigarOp
	num := 0
	haveDigits := false
	for i := 0; i < len(cigar); i++ {
		c := cigar[i]
		if c >= '0' && c <= '9' {
			num = num*10 + int(c-'0')
			haveDigits = true
			continue
		}
		if !haveDigits {
			return nil, fmt.Errorf("cigar %q: operation %q at offset %d has no length", cigar, string(c), i)
		}
		if !isCigarOp(c) {
			return nil, fmt.Errorf("cigar %q: invalid operation %q at offset %d", cigar, string(c), i)
		}
		ops = append(ops, CigarOp{Len: num, Op: c})
		num = 0
		haveDigits = false
	}
	if haveDigits {
		return nil, fmt.Errorf("cigar %q: trailing length %d with no operation", cigar, num)
	}
	return ops, nil
}

// CigarAlignEnd returns the 0-based, half-open end of an alignment: the first
// reference position past the aligned bases, for a record with the given 1-based
// POS and CIGAR. It uses the same 0-based half-open convention as ParseRegion
// and SamReader.Query. An unspecified ("*") CIGAR consumes no reference, so the
// result is pos-1.
func CigarAlignEnd(pos int, cigar string) int {
	return pos - 1 + CigarRefLen(cigar)
}

// CigarSoftClips returns the lengths of the leading and trailing soft clips
// (S operations) in SEQ coordinates: SEQ[:leading] and SEQ[len(SEQ)-trailing:]
// are the clipped bases. Hard clips (H) are ignored, since their bases are
// absent from SEQ — "5H10S80M" reports a leading soft clip of 10. A CIGAR whose
// only operation is S reports it as leading, never as both. Returns 0, 0 for an
// empty, unspecified ("*"), or malformed CIGAR.
func CigarSoftClips(cigar string) (leading, trailing int) {
	ops, err := ParseCigar(cigar)
	if err != nil {
		return 0, 0
	}
	i := 0
	for i < len(ops) && ops[i].Op == 'H' {
		i++
	}
	if i < len(ops) && ops[i].Op == 'S' {
		leading = ops[i].Len
		i++
	}
	j := len(ops) - 1
	for j >= 0 && ops[j].Op == 'H' {
		j--
	}
	// j >= i keeps a lone leading S from also being counted as trailing.
	if j >= i && ops[j].Op == 'S' {
		trailing = ops[j].Len
	}
	return leading, trailing
}

// CigarQueryToRef maps a 0-based query index (an index into SEQ as stored) to
// its 0-based reference position, for a record with the given 1-based POS and
// CIGAR.
//
// ok is false when the query base has no reference position: it falls in a soft
// clip (S) or an insertion (I), the CIGAR is empty, unspecified ("*"), or
// malformed, or qIdx is outside the query length implied by the CIGAR. When ok
// is false the returned position is 0 and must not be used. Hard clips (H) are
// skipped, since those bases are absent from SEQ, so qIdx always indexes SEQ
// directly.
//
// CigarQueryToRef is the inverse of the reference walk implied by CigarRefLen,
// and is what tools with query-indexed data (soft-clip analysis, poly(A)
// calling, base modification tags) need to reach reference coordinates.
func CigarQueryToRef(pos int, cigar string, qIdx int) (refPos int, ok bool) {
	if qIdx < 0 {
		return 0, false
	}
	ops, err := ParseCigar(cigar)
	if err != nil {
		return 0, false
	}
	q := 0         // query bases consumed so far
	ref := pos - 1 // 0-based reference cursor
	for _, op := range ops {
		if op.ConsumesQuery() {
			if qIdx < q+op.Len {
				if op.ConsumesRef() {
					return ref + (qIdx - q), true
				}
				return 0, false // S or I: no reference position
			}
			q += op.Len
		}
		if op.ConsumesRef() {
			ref += op.Len
		}
	}
	return 0, false // qIdx past the end of the query
}

// CigarRefLen returns the number of reference bases consumed by a CIGAR string.
// Operations M, D, N, =, X consume reference; I, S, H, P do not.
//
// This deliberately duplicates the scanning loop rather than building on
// ParseCigar: it runs once per record on the BAM read path and would otherwise
// allocate a slice per record. TestParseCigarConsistency pins the two against
// each other.
func CigarRefLen(cigar string) int {
	if cigar == "*" {
		return 0
	}
	refLen := 0
	num := 0
	for i := 0; i < len(cigar); i++ {
		c := cigar[i]
		if c >= '0' && c <= '9' {
			num = num*10 + int(c-'0')
		} else {
			switch c {
			case 'M', 'D', 'N', '=', 'X':
				refLen += num
			}
			num = 0
		}
	}
	return refLen
}

// CigarQueryLen returns the number of query (read) bases consumed by a CIGAR
// string. Operations M, I, S, =, X consume query bases; D, N, H, P do not.
// It returns 0 for an empty or unspecified ("*") CIGAR.
func CigarQueryLen(cigar string) int {
	if cigar == "*" || cigar == "" {
		return 0
	}
	queryLen := 0
	num := 0
	for i := 0; i < len(cigar); i++ {
		c := cigar[i]
		if c >= '0' && c <= '9' {
			num = num*10 + int(c-'0')
		} else {
			switch c {
			case 'M', 'I', 'S', '=', 'X':
				queryLen += num
			}
			num = 0
		}
	}
	return queryLen
}

// ValidateCigarSeq checks that a CIGAR string and SEQ are mutually consistent.
// When both are present (neither is "*" or empty), the query length implied by
// the CIGAR must equal len(seq). A mismatch means the record is malformed:
// encoders that reconstruct SEQ from the CIGAR (e.g. the CRAM writer) would
// otherwise silently drop bases, so callers should reject such records.
func ValidateCigarSeq(cigar, seq string) error {
	if cigar == "*" || cigar == "" || seq == "*" || seq == "" {
		return nil
	}
	if ql := CigarQueryLen(cigar); ql != len(seq) {
		return fmt.Errorf("cigar query length %d does not match sequence length %d", ql, len(seq))
	}
	return nil
}
