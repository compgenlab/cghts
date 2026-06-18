# hts

A Go library for computational genomics: sequence I/O, alignment, and native
SAM/BAM/CRAM/tabix handling. Particular focus on Oxford Nanopore (long-read)
sequencing workflows.

**Module:** `github.com/compgenlab/hts`

This is the library half of the former `cgkit` project; the CLI lives in
[`cgio`](https://github.com/compgenlab/cgio).

## Install

```bash
go get github.com/compgenlab/hts
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

### support packages

- **support/sequtils** — IUPAC ambiguity matching, reverse complement, homopolymer run analysis, 4-bit DNA encoding
- **support/utils** — `Semaphore` for concurrency control, `PositionTrackingReader`, float formatting
- **support/stringutils** — string helpers
- **analysis/seq** — GC content calculation
