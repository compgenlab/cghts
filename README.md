# cghts

[![Go Reference](https://pkg.go.dev/badge/github.com/compgenlab/cghts.svg)](https://pkg.go.dev/github.com/compgenlab/cghts)

A pure-Go library for computational genomics file formats: native readers and
writers for FASTA/FASTQ, SAM/BAM/CRAM, BGZF/tabix, BED, GTF, VCF, and
bigWig/bigBed, plus sequence alignment and analysis utilities. Minimal
dependencies and no cgo, so it builds and cross-compiles cleanly anywhere Go runs.

Full API documentation: <https://pkg.go.dev/github.com/compgenlab/cghts>

**Module:** `github.com/compgenlab/cghts`

This is the library half of the former `cgkit` project; the
[`cgkit`](https://github.com/compgenlab/cgkit) command-line toolkit is built on
it and exposes much of this functionality as CLI commands.

## Install

```bash
go get github.com/compgenlab/cghts
```

## Testing

```bash
make test     # GOCACHE=/tmp/go-build-cache go test ./...
```

## Packages

### seqio — FASTA/FASTQ I/O

Streaming readers and writers for FASTA and FASTQ files with transparent gzip support.

- `SeqReader` / `SeqRecord` interfaces for uniform access across formats
- `FastaReader` / `FastqReader` — lazy, streaming readers via `NextSeq()`; support indexed lookup by name
- `FastaWriter` / `FastqWriter` — writers with optional line wrapping (FASTA) and gzip output
- `SeqQual` — core type holding sequence, quality, name, strand, and position; supports `RevComp()` and `Sub()` extraction
- Memory-efficient chunked iteration via Go `iter.Seq`

### align — Pairwise and multiple sequence alignment

Smith-Waterman based alignment with affine gap penalties and Oxford Nanopore-aware homopolymer discounting.

- `NewLocalAligner()` — Smith-Waterman local alignment (soft clipping)
- `NewGlobalAligner()` — Needleman-Wunsch end-to-end alignment
- `NewSemiGlobalAligner()` — full query aligned, free target end gaps
- `DnaAlignmentDefaults()` / `OntAlignmentDefaults()` — preset scoring parameters
- Configurable scoring matrix, gap penalties, clipping, and homopolymer discount via builder pattern
- `AlignBatch()` — parallel alignment with semaphore-controlled goroutine pool
- `CigarCondense()` / `CigarExpand()` — convert between run-length and per-base CIGAR formats
- `MSA()` — incremental consensus multiple sequence alignment returning an `MSAAlignment` with optional homopolymer compression and reference sequence handling
- `MSAAlignment` — result type with `Consensus()`, `RehydratedConsensus()`, `WriteClustal()`, `WriteFasta()`, `GappedSequences()` for library-level output

### htsio — SAM/BAM/CRAM I/O

Native reading and writing of SAM, BAM, and tabix-indexed files. Samtools is only required for CRAM.

**Reading:**
- `SamReader` — interface with `Next()`, `Header()`, `Query()`, `Close()`
- `NewSamReader()` — auto-detects format: `.bam` → native BAM reader, `.sam`/`.sam.gz` → native text reader, `.cram` → samtools
- `Query(ref, start, end)` — returns `iter.Seq2[*SamRecord, error]` for indexed region queries (BAM via BAI, CRAM via samtools)
- Flag, MAPQ, and tag filtering via `SamReaderOpts`

**Writing:**
- `SamWriter` — interface with `Write()`, `Close()`
- `NewSamWriter()` — native BAM output (unsorted or coordinate/name sorted with merge sort), samtools for CRAM
- Sorted BAM writer buffers ~768MB, flushes to temp files, merge-sorts on Close

**Tabix:**
- `TabixReader` — query tabix-indexed BGZF files (BED, VCF, GFF) with TBI or CSI index auto-detection
- `TabixWriter` — sorted BGZF output with optional `.tbi` index generation; presets for BED, VCF, GFF
- Both use `iter.Seq2` for query results with 0-based half-open coordinates

**Index support:**
- BAI, TBI, CSI index parsers with shared `Query()` interface
- `ParseRegion()` — converts samtools-style region strings (`chr1:1000-2000`) to 0-based half-open

**Core types:**
- `SamRecord` — full SAM record with flag accessors (`IsUnmapped()`, `IsReverse()`, etc.) and typed tag access
- `SamHeader` — header manipulation including `@PG` line generation
- `TagFilter` — flexible tag-based filtering with comparison operators

### htsio/bgzf — BGZF compression

Low-level BGZF (Blocked GNU Zip Format) support used by BAM and tabix.

- `Reader` / `Writer` — streaming BGZF read/write with virtual offset tracking
- `IndexedReader` — random access with LRU block cache (default 64 blocks); supports virtual offset seeking and `.gzi` index for uncompressed offset seeking
- `NewBGZipFile()` — convenience constructor for file-backed BGZF output

### htsio/codec, htsio/bam, htsio/cram, htsio/sam, htsio/tabix

Format-specific subpackages backing the `htsio` facade — CRAM block codecs
(rANS, fqzcomp, arith), and the native BAM/SAM/CRAM/tabix reader and writer
implementations.

### htsio/bbi — bigWig / bigBed

Random-access reader for UCSC BBI files (bigWig and bigBed). BBI files are
self-indexed (chromosome B+ tree + spatial R-tree), so no sidecar index is
needed. The reader surface mirrors `htsio/tabix.Reader`.

- `Open()` → `Reader`; `Query(ref, start, end)` returns `iter.Seq2` of records
  over a 0-based half-open region (`Record.Value` for bigWig, `Record.Line` for bigBed)
- Base-resolution values only (zoom-level summaries ignored, so results are exact)
- Byte order auto-detected from the file magic; standard-library only

### bed — BED intervals

Streaming and tabix-indexed readers plus a writer for BED interval files.

- `BedRecord` — reference name, 0-based half-open `[start, end)`, optional BED6
  name/score/strand; columns past the sixth preserved verbatim in `Extras`
- `BedReader` — forward-only parser (`NextRecord`); file constructor auto-detects gzip
- `IndexedBedReader` — random access to a tabix-indexed BED via `Query()`
- `BedWriter` — sorted output with optional TBI/CSI index generation

### gtf — gene-model annotation

Parses a GTF file into an in-memory gene model (genes → transcripts →
exons/CDS/codons) with an interval index, and classifies genomic positions into
genic regions (coding exon, UTR, intron, junction, …). A port of ngsutilsj's
`GTFAnnotationSource`, reproducing its biotype derivation and region-code
precedence.

- `AnnotationSource` — in-memory model with position/region classifiers and gene iteration
- `IndexedAnnotationSource` — tabix-backed, per-position lookup with bounded memory
- Coordinates 0-based half-open (GTF's 1-based input converted on parse)

### vcf — Variant Call Format

Streaming and tabix-indexed readers, a writer, and a header/record model for VCF.

- **Lazy parsing** — a `VcfRecord` parses only CHROM/POS/REF up front; ID, ALT,
  QUAL, FILTER, INFO, FORMAT, and each sample column are parsed on first access
  and cached independently, so wide many-sample files stay cheap
- `VcfRecord.Pos` is 1-based (the one deliberate exception to the library's
  0-based half-open convention, matching the file for safe round-trips)
- **vcf/annotate** — composable framework writing INFO/FORMAT/ID fields onto
  records; locus annotators (`Indel`, `TsTv`, `AutoID`, `VariantDistance`, …) and
  sample annotators (`Dosage`, `VAF`, `FisherStrandBias`, …), plus bigWig/bigBed
  and GTF-backed annotators; run in order through a `Pipeline`
- **vcf/filter** — composable FILTER-stamping filters (comparison, list,
  flag-present/absent, zygosity, chrom, indel), chained via `Chain`; a port of
  ngsutilsj's vcf/filter framework

### varstore — sparse genotype store

Cohort genotypes in Parquet, plus a VCF-backed implementation of the same
interface. A store is three files sharing a base name, and all three are
required:

```
BASE.calls.parquet     one row per ALT-carrying genotype
BASE.sites.parquet     one row per interrogated site, with AC/AN
BASE.regions.parquet   runs of catalog sites at which a sample was called
```

Only the ALT-carrying genotypes are stored — a real cohort callset is
overwhelmingly reference, so this is a small fraction of the dense matrix. The
other two files are what make the absent rows interpretable: without them a
missing row cannot be told apart from a position nobody ever looked at, which is
the distinction the four states exist to draw.

One method queries either backend, and the two must agree — the same question
against a VCF and against a store converted from it:

```go
s, _ := varstore.OpenParquet("cohort")   // or OpenVcf, or a URL/s3:// prefix
defer s.Close()

calls, err := s.Calls(varstore.Query{
    Loci:    []varstore.Locus{{Chrom: "chr17", Pos: 43045703, Ref: "A", Alt: "G"}},
    Samples: []string{"HG00096"},        // empty means every sample
    Gate:    varstore.Gate{MinDP: 10},
})
for c, err := range calls { … }
```

Sites and samples are independent axes, and an empty selector restricts neither,
so the zero `Query` is every genotype in the store. A query is always **one
pass**: the Parquet side prunes row groups by the union of the selectors, so a
1000-variant panel costs about what one variant does.

- `IncludeRef` also emits `0/0` calls, which needs `sites` and `regions` —
  a store lacking them fails with `ErrNotClassifiable` rather than returning an
  ALT-only stream that would read as "nobody is reference here"
- **Run intervals are not coverage.** They compress a per-site record of "this
  sample was called at these catalog sites" and claim nothing about the bases
  between them, so a locus absent from `sites` is unassayed for every sample even
  where runs bracket it. A gVCF's reference blocks (`END`, `MIN_DP`) *are* span
  assertions, and only such a store could answer off-catalog — hence
  `SpanSemantics`
- **A span selects by overlap, not by start position.** A record covers more than
  its POS whenever REF is longer than one base, or a symbolic ALT or gVCF block
  declares a span, so a query for `chr1:3100-3150` returns a deletion beginning at
  3000. The `ref_end` column carries that extent; the writer derives it from
  `len(REF)` when a caller does not set it, and `(*vcf.VcfRecord).RefSpan` computes
  the full version from `INFO/END`, `SVLEN` and `FORMAT/LEN`. Stores written before
  the column existed select on position alone and report
  `ParquetStore.HasRefSpans() == false` — a different answer, not merely a slower
  one, so it is worth surfacing
- **gVCF reference blocks read as coverage.** A block asserts that a sample was
  called reference across a span at at least `MIN_DP` -- the one claim a plain VCF
  cannot make, and what lets a query answer for a position no variant record
  mentions. `OpenVcf` on a gVCF yields one reference row per block per sample, never
  one per base, with the block's `MIN_DP` as the depth bound and `Alt` `.`; the gap
  between blocks stays unanswerable, because nothing was reported there
- Stores keep the source's own contig spelling; compare through
  `CanonKey`/`SameChrom` so `chr17`, `17` and `NC_000017.11` all resolve

`cgkit`'s `vcf-toparquet` and `vcf-varquery` are the CLI over this package; see
its [`docs/vcf-toparquet.md`](https://github.com/compgenlab/cgkit/blob/main/docs/vcf-toparquet.md)
for the format walk-through and measurements.

### iosource — pluggable random-access I/O

Random-access byte sources for genomic files behind a concurrency-safe
`io.ReaderAt`, so index-driven readers fetch only the byte ranges they need.

- Local-file and **HTTP(S)-Range** implementations (standard-library only), so a
  BAM/CRAM/VCF/bigWig can be queried remotely without downloading the whole file
- `ByteSource` interface for other transports (SFTP, S3, …)
- Sibling-file resolution (`.bai`/`.tbi`/…) over both local and HTTP locators
- `ReadSeeker` adapts a `ByteSource` for readers that decompress sequentially
- **`iosource.Open(ctx, locator)`** dispatches on scheme: a plain path opens a
  file, `http(s)://` uses Range requests, and other schemes come from registered
  transports

**S3** lives in `iosource/s3`. Importing it registers the scheme, so a blank
import is all it takes:

```go
import _ "github.com/compgenlab/cghts/iosource/s3"

src, _ := iosource.Open(ctx, "s3://bucket/clinvar.vcf.gz")
```

Credentials come from the standard AWS chain — environment, `~/.aws/credentials`
and `~/.aws/config` (honouring `AWS_PROFILE`), then a container or EC2 instance
role — so SSO profiles, `role_arn`+`source_profile` chains and
`credential_process` all work. `s3.WithProfile`, `WithRegion` and `WithEndpoint`
override explicitly; `AWS_ENDPOINT_URL` targets an S3-compatible gateway.

**Wired into tabix and bbi.** `tabix.NewReaderFromSource(data, index, …)` takes a
seekable data stream and a raw `.tbi`/`.csi` stream, detecting which index format
it was given from the magic; `bbi.NewReaderFromSource(readerAt, …)` does the same
for bigWig/bigBed, which are pure random access. `LoadTBIFrom` / `LoadCSIFrom`
parse an index from any reader. The path-based constructors are unchanged.

```go
src, _ := iosource.NewHTTPRange(url)          // or any ByteSource
rs, _ := iosource.ReadSeeker(src)
idx, _, _ := iosource.ResolveSibling(url, []string{".tbi", ".csi"}, iosource.HTTPSibling)
r, _ := tabix.NewReaderFromSource(rs, idx, tabix.WithCloser(src))
defer r.Close()
```

A narrow query over a 4.75 MB indexed VCF fetches ~1.8% of the file.

**BAM and CRAM read remotely too.** `htsio.OpenSamReader(ctx, locator)` opens a
SAM/BAM/CRAM from a path, a URL or an `s3://` locator, resolving the `.bai` or
`.crai` through the same transport, so indexed region queries work without
downloading the file. For CRAM the reference is independent of where the data
lives — `SetRefPath` takes any locator and `SetRefReader` an already-open one,
so a remote CRAM with a local reference (or the reverse) is fine.

**varstore reads remotely too.** `OpenParquetContext(ctx, locator)` opens a
Parquet store from a path, a URL or an `s3://` prefix. Parquet suits this
unusually well: the footer carries per-row-group statistics and the package
already prunes groups by them, so a locus query skips pruned groups without
transferring them at all. Store members are held open for the store's lifetime
rather than reopened per query — which also removes a repeated footer parse on
local reads.

### support packages

- **support/sequtils** — IUPAC ambiguity matching, reverse complement, homopolymer run analysis, 4-bit DNA encoding
- **support/stats** — 2×2 Fisher exact test, Phred/log2 conversions
- **support/utils** — `Semaphore` for concurrency control, `PositionTrackingReader`, float formatting
- **support/stringutils** — string helpers
- **analysis/seq** — GC content calculation
