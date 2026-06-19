// Package bed provides streaming and indexed readers plus a writer for BED
// (Browser Extensible Data) interval files.
//
// # Records
//
// [BedRecord] is the core value type: a reference name, a 0-based half-open
// [start, end) span, and the optional BED6 fields name, score, and strand.
// Columns beyond the sixth are preserved verbatim as opaque strings in
// [BedRecord.Extras]; the package does not interpret BED12 block (exon)
// structure. The HasName flag distinguishes a bare BED3 record (where the
// name/score/strand fields are suppressed on output) from a full record.
//
// # Readers
//
// [BedReader] is a streaming, forward-only parser. Open one over an
// [io.Reader] with [NewBedReader] or over a file with [NewBedFile]; the
// file-based constructor detects gzip input by its magic bytes and
// transparently decompresses. Records are pulled one at a time with
// [BedReader.NextRecord], which returns [io.EOF] at end of input.
//
// [IndexedBedReader] provides random access to a tabix-indexed BED file
// (BGZF-compressed with a companion .tbi or .csi index). [IndexedBedReader.Query]
// returns the records overlapping a 0-based half-open region.
//
// # Writer
//
// [BedWriter] writes records one at a time. Open one over an [io.Writer] with
// [NewBedWriter] or over a file with [OpenBedWriter]. [BedWriterOpts] selects
// the output column layout (auto, strict BED3, or strict BED6), whether scores
// are coerced to integers, and whether output is written as a tabix-indexed
// BGZF file. When the filename ends in ".gz" (and tabix output is not
// requested) the output is whole-file gzip-compressed.
//
// All coordinates are 0-based half-open ([start, end)), matching both the BED
// format itself and the rest of this library.
package bed
