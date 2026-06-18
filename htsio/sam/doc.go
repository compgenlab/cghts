// Package sam provides plain-text SAM input and output, plus a samtools-backed
// reader for indexed and binary formats.
//
// All readers in this package implement the [htsio.SamReader] interface and
// produce [htsio.SamRecord] values, while [Writer] implements the
// [htsio.SamWriter] interface. Records are streamed lazily through Go 1.23
// range-over-func iterators (iter.Seq2[*htsio.SamRecord, error]).
//
// # Main types
//
//   - [TextReader] — reads plain SAM text from a file or io.ReadCloser.
//     Gzip-compressed input (.sam.gz) is auto-detected by peeking at the
//     first two bytes, so the same reader handles both compressed and
//     uncompressed streams. Because plain SAM has no index, region queries
//     are not supported.
//   - [SamtoolsReader] — reads SAM, BAM, or CRAM by shelling out to the
//     samtools view subprocess. It forwards thread, SAM-flag, and
//     mapping-quality options to samtools and, given an index, supports
//     region queries.
//   - [Writer] — writes records as SAM text to a file (or stdout when the
//     filename is "-"). The header is emitted lazily before the first
//     record, and the writer is safe for concurrent use.
//
// On construction, the package registers [TextReader] as the fallback
// [htsio] reader, so any input format not claimed by a more specific reader
// (such as BAM or CRAM) is read as SAM text.
//
// # Parsing and filtering
//
// Each tab-delimited data line must have at least the 11 mandatory SAM
// fields; the numeric fields (FLAG, POS, MAPQ, PNEXT, TLEN) are validated and
// a malformed value yields a parse error. Optional TAG:TYPE:VALUE fields are
// parsed into the record's tag map in order; fields that are not well-formed
// are skipped, matching samtools' tolerant behavior. Both reader types apply
// the filters carried by [htsio.SamReaderOpts] (required/excluded SAM flags
// and a minimum mapping quality), silently skipping records that do not pass.
//
// External processes are required only by [SamtoolsReader], which expects a
// samtools binary on the PATH; the text reader and writer have no external
// dependencies.
package sam
