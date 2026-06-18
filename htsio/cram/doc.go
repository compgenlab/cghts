// Package cram implements native reading and writing of CRAM files, the
// reference-based compressed alignment format used as a more compact
// alternative to BAM.
//
// The package speaks the CRAM container format directly rather than shelling
// out to samtools: it parses the file definition, containers, slices, blocks,
// the per-container compression header, and the bit-packed data series, and it
// reconstructs SEQ, CIGAR, quality strings, and optional tags from the encoded
// read-vs-reference differences. Writing performs the inverse, grouping records
// into slices and containers and emitting the appropriate block codecs.
//
// # Coordinate conventions
//
// All query coordinates exposed by [Reader.Query] are 0-based half-open
// ([start, end)), matching the rest of htsio. On-disk CRAM positions are 1-based
// and are converted internally. Records are decoded into htsio.SamRecord values
// so callers see the same record type as the BAM and SAM readers.
//
// # Versions
//
// CRAM v2.x and v3.x are supported. In v3 and later, container and block
// headers are protected by a trailing CRC32 that is verified on read, and the
// record counter is stored as LTF8 rather than ITF8. The writer emits [V31] by
// default; see [V2], [V3], and [V31].
//
// # Variable-length integers
//
// CRAM encodes integers with two self-describing variable-length schemes whose
// byte count is determined by leading set bits in the first byte:
//
//   - ITF8 encodes a 32-bit value in 1-5 bytes.
//   - LTF8 extends the scheme to a 64-bit value in 1-9 bytes.
//
// These appear throughout container, slice, and block headers as well as inside
// the data series.
//
// # Index support
//
// Random access uses a CRAI index. [Reader.Query] lazily loads the gzip-compressed
// "<file>.crai" sidecar, finds the slices overlapping the requested region, seeks
// to each container, and yields only the records that overlap the region.
// [Reader.Records] streams every record in the file in storage order without
// requiring an index, stopping cleanly at the CRAM EOF container.
//
// # Main types
//
//   - [Reader] reads CRAM files and implements htsio.SamReader. Construct one
//     with [NewReader] or [NewReaderFromStream].
//   - [Writer] writes CRAM files and implements htsio.SamWriter. Construct one
//     with [NewWriter] or [NewWriterFromWriter].
//   - [WriterOpts] (built via [NewWriterOpts]) configures the version, reference,
//     compression level, records per slice, and reference embedding.
//   - [Version] identifies the CRAM version to emit.
package cram
