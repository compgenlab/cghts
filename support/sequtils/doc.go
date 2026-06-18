// Package sequtils provides low-level DNA sequence utilities.
//
// The package works on plain ASCII byte sequences and supports the full set of
// IUPAC nucleotide ambiguity codes. It includes:
//
//   - [ConvertDNATo4Bit] for encoding a base as a 4-bit set of possible
//     nucleotides, the basis for ambiguity-aware comparisons.
//   - [DNAMatches] for testing whether two bases are compatible under IUPAC
//     rules, case-insensitively.
//   - [ReverseComplement] for computing the reverse complement of a sequence.
//   - [HomopolymerRunLen] and [HomopolymerCompress] for analyzing and
//     collapsing homopolymer runs.
package sequtils
