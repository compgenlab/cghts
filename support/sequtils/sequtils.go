package sequtils

import (
	"strings"
)

/*
Convert a DNA base to a number for easier ambiguity checking

bits:
3210
TGCA

	0x1 A : Adenine
	0x2 C : Cytosine
	0x4 G : Guanine
	0x8 T : Thymine
	0x5 R : puRine (A or G)
	0xA Y : pYrimidine (C or T)
	0x3 M : aMino (A or C)
	0xC K : Keto (G or T)
	0x6 S : Strong interaction (G or C)
	0x9 W : Weak interaction (A or T)
	0xE B : not A (C, G, or T)
	0xD D : not C (A, G, or T)
	0xB H : not G (A, C, or T)
	0x7 V : not T (A, C, or G)
	0xF N : Any nucleotide (A, C, G, or T)
*/
func ConvertDNATo4Bit(r byte) byte {
	switch r {
	case 'A':
		return 0x01
	case 'C':
		return 0x02
	case 'G':
		return 0x04
	case 'T':
		return 0x08
	case 'R':
		return 0x05
	case 'Y':
		return 0x0A
	case 'M':
		return 0x03
	case 'K':
		return 0x0C
	case 'S':
		return 0x06
	case 'W':
		return 0x09
	case 'B':
		return 0x0E
	case 'D':
		return 0x0D
	case 'H':
		return 0x0B
	case 'V':
		return 0x07
	case 'N':
		return 0x0F
	}
	return 0x0
}

/*
Return the complementary DNA base for a given char/rune

	Complementary Ambiguity Codes
	A-T / C-G
	R (A/G) is complementary to Y (C/T)
	M (A/C) is complementary to K (G/T)
	B (not A) is complementary to V (not T)
	D (not C) is complementary to H (not G)
	N is complementary to N
*/
func dnaComplement(r byte) byte {
	switch r {
	case 'A':
		return 'T'
	case 'C':
		return 'G'
	case 'G':
		return 'C'
	case 'T':
		return 'A'
	case 'R':
		return 'Y'
	case 'Y':
		return 'R'
	case 'M':
		return 'K'
	case 'K':
		return 'M'
	case 'S':
		return 'W'
	case 'W':
		return 'S'
	case 'B':
		return 'V'
	case 'D':
		return 'H'
	case 'H':
		return 'D'
	case 'V':
		return 'B'
	case 'N':
		return 'N'
	case 'a':
		return 't'
	case 'c':
		return 'g'
	case 'g':
		return 'c'
	case 't':
		return 'a'
	case 'r':
		return 'y'
	case 'y':
		return 'r'
	case 'm':
		return 'k'
	case 'k':
		return 'm'
	case 's':
		return 'w'
	case 'w':
		return 's'
	case 'b':
		return 'v'
	case 'd':
		return 'h'
	case 'h':
		return 'd'
	case 'v':
		return 'b'
	case 'n':
		return 'n'
	}
	return r
}

func ReverseComplement(seq string) string {
	b := []byte(seq)
	for i := 0; i < len(seq); i++ {
		b[i] = dnaComplement(seq[len(seq)-i-1])
	}
	return string(b)
}

func byteToUpperASCII(b byte) byte {
	if 'a' <= b && b <= 'z' {
		return b ^ 0x20 // equivalent to b - 'a' + 'A' (offset is 32, so this effectively subtracts 0x20, or 32d)
	}
	return b
}

func DNAMatches(one byte, two byte) bool {
	one = byteToUpperASCII(one)
	two = byteToUpperASCII(two)

	if one == two {
		return true
	}

	oneB := ConvertDNATo4Bit(one)
	twoB := ConvertDNATo4Bit(two)

	// AC  & A
	// 0x3 & 0x1 = 0x1

	if oneB&twoB > 0 {
		return true
	}
	return false
}

// HomopolymerRunLen returns, for each base position i, the total length
// of the homopolymer run that seq[i] belongs to.
// Example: "AAATCC" -> [3 3 3 1 2 2]
func HomopolymerRunLen(seq string) []int {
	n := len(seq)
	out := make([]int, n)
	if n == 0 {
		return out
	}

	start := 0
	for start < n {
		end := start + 1
		for end < n && seq[end] == seq[start] {
			end++
		}
		runLen := end - start
		for i := start; i < end; i++ {
			switch seq[start] {
			case 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't':
				out[i] = runLen
			default:
				// non-standard base, mark with negative run length
				out[i] = -runLen
			}
		}
		start = end
	}
	return out
}

// HomopolymerCompress returns a compressed version of the sequence
// where homopolymer runs are compressed to a single base.
// Example: "AAATCC" -> "ATC", []int{3, 1, 2}
func HomopolymerCompress(seq string) (string, []int) {
	var out strings.Builder
	n := len(seq)
	if n == 0 {
		return out.String(), []int{}
	}

	out.WriteByte(seq[0])
	runLens := []int{1}
	for i := 1; i < n; i++ {
		if seq[i] != seq[i-1] {
			out.WriteByte(seq[i])
			runLens = append(runLens, 1)
		} else {
			runLens[len(runLens)-1]++
		}
	}
	return out.String(), runLens
}
