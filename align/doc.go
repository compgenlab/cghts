// Package align provides pairwise and multiple sequence alignment for DNA,
// with scoring tuned for both Illumina short reads and Oxford Nanopore (ONT)
// long reads.
//
// # Pairwise alignment
//
// The core engine is a Smith-Waterman style dynamic-programming aligner with
// affine gap penalties (separate open and extend costs, configured
// independently for insertions and deletions). Three modes are available, each
// constructed from a shared [alignmentOptions] value:
//
//   - [NewLocalAligner] — local alignment with soft clipping. The best-scoring
//     sub-alignment is returned and the unaligned ends of the query are
//     reported as soft clips (S) in the CIGAR.
//   - [NewGlobalAligner] — global (Needleman-Wunsch style) alignment in which
//     both the query and target are aligned end to end. Clipping is disabled.
//   - [NewSemiGlobalAligner] — the query is aligned end to end while the target
//     may have free end gaps on both sides, useful for placing a read inside a
//     longer reference or consensus.
//
// Each aligner satisfies the [PairwiseAligner] interface and returns a
// [PairwiseAlignment] describing the aligned coordinates (0-based, half-open on
// the end), the alignment score, and a CIGAR string.
//
// # Scoring and gap penalties
//
// Match/mismatch scoring is supplied through the [ScoringMatrix] interface;
// [MatchMismatchScoring] provides a simple match-reward / mismatch-penalty
// matrix that honors IUPAC ambiguity codes. [DnaAlignmentDefaults] and
// [OntAlignmentDefaults] return ready-made option sets, which can be further
// customized with the fluent setters on the returned value (for example
// GapPenaltyIns, GapPenaltyDel, and ClippingPenalty).
//
// # ONT homopolymer discounts
//
// Oxford Nanopore reads frequently miscall the length of homopolymer runs. To
// avoid over-penalizing these indels, the aligner can discount gap penalties
// inside homopolymer regions based on the run length. The discount grows with
// the log of the run length and is capped, so a long run makes a single indel
// nearly free. [OntAlignmentDefaults] enables this; HomopolymerDiscount tunes
// it; [DnaAlignmentDefaults] disables it (homopolymer errors are uncommon on
// short reads).
//
// # Parallel batch alignment
//
// [AlignBatch] aligns every query against every target concurrently, bounded
// by a [github.com/compgenlab/hts/support/utils.Semaphore]. It returns a
// [PairwiseAlignmentPromise]; calling its Result method blocks until all
// alignments finish and returns the single highest-scoring [PairwiseAlignment].
//
// # Multiple sequence alignment
//
// [MSA] builds a multiple sequence alignment from a set of sequences using an
// incremental, consensus-guided strategy: all pairwise alignments are computed,
// the best-scoring pair seeds the alignment, and remaining sequences are added
// one at a time by aligning each to the running consensus. The result is an
// [MSAAlignment] stored column-major; it supports majority-vote consensus,
// gapped-sequence and CLUSTAL/FASTA output, optional homopolymer
// compress/expand, and an optional designated reference row. Options are built
// with [NewMSAOptions].
//
// # CIGAR helpers
//
// CIGAR strings use the standard operations M (match/mismatch), I (insertion
// relative to the target), D (deletion relative to the target), and S (soft
// clip). [CigarCondense] converts a per-base (expanded) CIGAR such as
// "IIMMMMMDMM" into run-length form ("2I5M1D2M"); [CigarExpand] performs the
// reverse.
package align
