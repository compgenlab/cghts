// Package tabix reads and writes tabix-indexed, BGZF-compressed text files
// such as BED, VCF, and GFF/GTF, providing random access to records by
// genomic region.
//
// A tabix-indexed file consists of a BGZF-compressed, position-sorted text
// file plus a companion index that maps genomic regions to BGZF virtual
// offsets (an offset into the compressed stream paired with an offset into the
// decompressed block). Two index formats are supported:
//
//   - TBI: the classic tabix index, itself BGZF-compressed, using the standard
//     6-level binning scheme plus a 16 kb linear index. It cannot describe
//     reference sequences longer than 512 Mb.
//   - CSI: the coordinate-sorted index, also BGZF-compressed, using a
//     configurable min_shift and binning depth in place of the linear index,
//     allowing it to index sequences longer than 512 Mb.
//
// All region coordinates exposed by this package are 0-based, half-open
// [start, end), regardless of the on-disk coordinate convention recorded in
// the index metadata.
//
// The main exported types are:
//
//   - [Reader]: opens a BGZF file together with its .tbi or .csi index and
//     returns [Record] values overlapping a queried region.
//   - [Writer] and [WriterOpts]: sort, BGZF-compress, and (optionally) index
//     tabular input, with presets for [WriterOpts.BED], [WriterOpts.VCF], and
//     [WriterOpts.GFF].
//   - [BinIndex]: a parsed BAI or TBI index, loaded via [LoadBAI] or
//     [LoadTBI].
//   - [CSIIndex]: a parsed CSI index, loaded via [LoadCSI].
//   - [Chunk]: a contiguous range of BGZF virtual offsets returned by a query.
//
// Both index types answer region queries with Query, returning a sorted,
// merged list of [Chunk] ranges to scan in the data file.
package tabix
