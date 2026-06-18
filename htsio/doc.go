// Package htsio provides native Go readers and writers for high-throughput
// sequencing alignment files in the SAM, BAM, CRAM, and tabix formats.
//
// SAM (text) and BAM (BGZF-compressed binary) are read and written natively,
// without any external dependencies. CRAM support is delegated to the samtools
// executable, which must be available on the PATH; reading or writing CRAM
// therefore requires a working samtools installation.
//
// # Reading
//
// [NewSamReader] opens a file and auto-detects the format from its magic bytes
// (raw bytes for CRAM, post-gzip bytes for BAM, falling back to SAM text):
//
//	r, err := htsio.NewSamReader("alignments.bam")
//	if err != nil {
//		return err
//	}
//	defer r.Close()
//	for rec, err := range r.Records() {
//		if err != nil {
//			return err
//		}
//		fmt.Println(rec.ReadName, rec.RefName, rec.Pos)
//	}
//
// [NewSamReaderFromReader] performs the same detection on an arbitrary
// io.ReadCloser (for example, os.Stdin); readers built from a stream do not
// support [SamReader.Query]. Both constructors accept optional [SamReaderOpts]
// to apply flag, mapping-quality, and tag filters during iteration.
//
// The [SamReader] interface exposes records via Go 1.23 range-over-func
// iterators ([iter.Seq2]). A reader is single-pass and not safe for concurrent
// iteration: [SamReader.Records] and [SamReader.Query] advance the same
// underlying stream, so only one iterator may be active at a time. Iterators
// may be stopped early with a break.
//
// # Querying by region
//
// For indexed BAM files, [SamReader.Query] returns an iterator over the records
// overlapping a region. [ParseRegion] converts samtools-style region strings
// such as "chr1:1,000-2,000" into the coordinates Query expects.
//
// # Writing
//
// The [SamWriter] interface accepts [SamRecord] values via Write and finalizes
// the file (including the BGZF/CRAM EOF marker) on Close. Format-specific
// writer constructors live in the bam, sam, and cram subpackages.
//
// # Records, headers, and tags
//
// [SamRecord] holds the alignment fields plus parsed optional fields as
// [SamTag] values keyed by tag name. Flag accessor methods such as
// [SamRecord.IsUnmapped] and [SamRecord.IsReverse] decode the bitwise SAM flag.
// [SamHeader] holds the raw header lines and provides helpers to extract
// reference sequences ([SamHeader.References]), read groups, and MD5 checksums.
// [TagFilter] and the Tag* comparison operators express tag-based record
// filtering, which [SamReaderOpts] applies during iteration.
//
// # Coordinate convention
//
// All region query coordinates ([SamReader.Query], [ParseRegion]) are 0-based
// half-open: the interval [start, end) includes start and excludes end. Note
// that the [SamRecord.Pos] field, in contrast, is 1-based, matching the SAM
// text representation.
package htsio
