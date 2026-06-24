// Package vcf provides streaming and tabix-indexed readers, a writer, and a
// header/record model for VCF (Variant Call Format) files.
//
// # Lazy parsing
//
// Unlike a strict eager decoder, a [VcfRecord] parses only the cheap leading
// columns (CHROM, POS, REF) when it is read. The remaining columns — ID, ALT,
// QUAL, FILTER, INFO, FORMAT, and the per-sample genotype columns — are parsed
// on first access and cached. INFO is parsed only when [VcfRecord.Info] (or a
// convenience accessor) is called, and each sample column is parsed
// independently, so reading one sample's values never forces parsing the
// others. Many analyses need only CHROM/POS/REF/ALT (or only need to move the
// sample columns around), and for those the per-record cost stays small even on
// very wide, many-sample files.
//
// Because a record retains a reference to its raw line and may parse lazily
// after iteration has advanced, the reader allocates a fresh [VcfRecord] per
// call (it does not reuse a single record). A [VcfRecord] is not safe for
// concurrent use.
//
// # Coordinates
//
// VCF is intrinsically 1-based, and [VcfRecord.Pos] is reported 1-based to
// match the file and avoid off-by-one errors on round-trip. This is the one
// deliberate exception to the rest of the hts library, which uses 0-based
// half-open coordinates. Region queries ([IndexedVcfReader.Query]) take 0-based
// half-open arguments like every other indexed reader; only the record field is
// 1-based. Use [VcfRecord.ZeroBasedStart] when emitting BED-style output.
//
// # Readers and writers
//
// [NewVcfReader] wraps an io.Reader; [NewVcfFile] opens a file and transparently
// gunzips it. [NewIndexedVcfReader] provides random access over a
// BGZF-compressed, tabix-indexed VCF. [VcfWriter] writes records back out,
// preserving the original line verbatim when a record is unmodified.
package vcf
