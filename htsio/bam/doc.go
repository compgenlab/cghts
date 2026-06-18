// Package bam provides a native BAM reader and writer built directly on the
// BGZF (Blocked GNU Zip) layer, with no external dependencies (no samtools).
//
// BAM is the binary, BGZF-compressed encoding of SAM alignment data. This
// package decodes BAM records into and encodes them from the shared
// [github.com/compgenlab/hts/htsio.SamRecord] type, so callers work with the
// same record and header model used across the htsio packages.
//
// # Main types
//
//   - [Reader] reads a BAM stream: it parses the BAM magic, header text, and
//     reference dictionary, then exposes records via [Reader.Records]. When the
//     reader is file-backed and a companion .bai index is present,
//     [Reader.Query] returns only the records overlapping a region. Reader
//     implements the htsio.SamReader interface.
//   - [Writer] writes a BAM stream: it emits the header on the first write and
//     encodes records on a background goroutine, so [Writer.Write] is safe for
//     concurrent use. Writer implements the htsio.SamWriter interface.
//   - The sorted writer (see [NewSortedWriter] and [NewSortedWriterFromWriter])
//     buffers records in memory, spills sorted temporary BAM files to disk when
//     the buffer fills, and merge-sorts them into the final output on Close.
//     It supports both coordinate and read-name sort orders.
//
// # Coordinate conventions
//
// On the wire BAM positions are 0-based, while [github.com/compgenlab/hts/htsio.SamRecord]
// fields (Pos, PosNext) are 1-based to match SAM; the reader and writer convert
// between the two. [Reader.Query] region arguments are 0-based half-open
// [start, end), matching the rest of the htsio query API.
//
// # Compression
//
// All output is BGZF-compressed. The writer constructors that take a threads
// argument use parallel BGZF compression when threads > 1; otherwise
// compression is single-threaded.
package bam
