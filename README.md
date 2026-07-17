# cghts

[![Go Reference](https://pkg.go.dev/badge/github.com/compgenlab/cghts.svg)](https://pkg.go.dev/github.com/compgenlab/cghts)

A Go library for computational genomics: sequence I/O, alignment, and native
SAM/BAM/CRAM/tabix handling. Particular focus on Oxford Nanopore (long-read)
sequencing workflows.

Full API documentation: <https://pkg.go.dev/github.com/compgenlab/cghts>

**Module:** `github.com/compgenlab/cghts`

This is the library half of the former `cgkit` project. The CLIs built on it live
in separate repos: [`cgkit`](https://github.com/compgenlab/cgkit), the general
genomics CLI, and `nupa`, a focused Oxford Nanopore UMI and poly(A) toolkit.

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

### iosource — pluggable random-access I/O

Random-access byte sources for genomic files behind a concurrency-safe
`io.ReaderAt`, so index-driven readers fetch only the byte ranges they need.

- Local-file and **HTTP(S)-Range** implementations (standard-library only), so a
  BAM/CRAM/VCF/bigWig can be queried remotely without downloading the whole file
- `ByteSource` interface for other transports (SFTP, S3, …)
- Sibling-file resolution (`.bai`/`.tbi`/…) over both local and HTTP locators

### support packages

- **support/sequtils** — IUPAC ambiguity matching, reverse complement, homopolymer run analysis, 4-bit DNA encoding
- **support/stats** — 2×2 Fisher exact test, Phred/log2 conversions
- **support/utils** — `Semaphore` for concurrency control, `PositionTrackingReader`, float formatting
- **support/stringutils** — string helpers
- **analysis/seq** — GC content calculation
