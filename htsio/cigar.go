package htsio

import "fmt"

// CigarRefLen returns the number of reference bases consumed by a CIGAR string.
// Operations M, D, N, =, X consume reference; I, S, H, P do not.
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
