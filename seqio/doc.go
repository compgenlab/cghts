// Package seqio provides streaming readers and writers for biological
// sequence files (FASTA and FASTQ) along with random-access readers for
// reference genomes.
//
// # Sequence values
//
// [SeqQual] is the core value type. It carries a sequence string, an
// optional quality string, a name, a position, and a strand orientation.
// SeqQual values are immutable: methods such as [SeqQual.RevComp] and
// [SeqQual.Sub] return new values rather than mutating the receiver.
//
// A [SeqRecord] represents a single FASTA/FASTQ entry. Records expose their
// name and comment and can yield their sequence either all at once with
// FullSeq or in bounded pieces with Chunks, which lets very large sequences
// (whole chromosomes) be processed without loading them entirely into memory.
//
// # Streaming readers and writers
//
// [FastaReader] and [FastqReader] are streaming, forward-only parsers. Open
// one over an [io.Reader] with [NewFastaReader]/[NewFastqReader], or over a
// file with [NewFastaFile]/[NewFastqFile]. The file-based constructors detect
// gzip input by its magic bytes and transparently decompress. Records are
// pulled one at a time with NextSeq; the readers are not seekable, so the
// Names and FetchRecord helpers consume (and then close) the underlying
// stream.
//
// [FastaWriter] and [FastqWriter] write records to an [io.Writer] or a file.
// The Open* constructors gzip-compress output when the filename ends in
// ".gz". [FastaWriter] supports optional line wrapping via [FastaWriterOpts].
//
// # Reference sequences
//
// The [ReferenceReader] interface provides random access to named reference
// sequences by 0-based half-open coordinate range. [OpenReference]
// auto-detects the source and returns an appropriate implementation:
//
//   - an HTTP/HTTPS URL backed by a remote .fai index and Range requests
//     ([RemoteFastaReader]),
//   - a local indexed FASTA (.fai index, plain or bgzip; [IndexedFastaReader]),
//   - or a plain unindexed FASTA loaded into memory as compressed chunks.
//
// Additional implementations resolve sequences by MD5 checksum: [RefgetReader]
// fetches from a GA4GH refget server, and [RefCacheReader] follows the htslib
// REF_PATH/REF_CACHE conventions. Indexed and remote readers load data in
// 10MB chunks behind an LRU cache (1GB max) and are safe for concurrent use.
//
// All reference coordinates are 0-based and half-open ([start, end)); returned
// bases are uppercased and out-of-range coordinates are clamped to the
// sequence boundaries.
package seqio
