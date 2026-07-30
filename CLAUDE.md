# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`cghts` is a pure-Go library for computational genomics file formats: native
readers and writers for FASTA/FASTQ, SAM/BAM/CRAM, BGZF/tabix, BED, GTF, VCF, and
bigWig/bigBed, plus sequence alignment and analysis utilities. No cgo; the only
third-party dependency is `ulikunitz/xz`. It is the library half of the former
`cgkit` project.

The CLIs that consume it live in separate repos: `cgkit`
(`github.com/compgenlab/cgkit`), the general-purpose genomics CLI, and `nupa`
(`github.com/compgenlab/nupa`, private), a focused Oxford Nanopore UMI and poly(A)
toolkit. The ONT-specific commands (`trim`, `umi-cluster`, `umi-dedup`,
`polya-site`) live in nupa, so this library stays domain-general — any
Nanopore-specific handling here (e.g. the aligner's homopolymer discounts) is a
reusable feature, not the library's focus.

**Module:** `github.com/compgenlab/cghts`
**Go version:** 1.24 (raised from 1.23 by `parquet-go`, which `varstore/` requires)

## Commands

```bash
# Run all tests
make test
# equivalent to:
GOCACHE=/tmp/go-build-cache go test ./...

# Run a single test
go test ./align/... -run TestCigarCondense
```

When developing alongside `cgkit` and `nupa`, the modules are joined by a
`go.work` workspace in the parent directory so those CLIs resolve this checkout
directly. This Makefile exports `GOWORK=off`, so the library always builds
standalone regardless of the ambient workspace.

## Architecture

### Package Layout

- **`seqio/`** — FASTA/FASTQ readers with gzip support. Core type is `SeqQual`, which holds sequence, quality scores, name, and strand. Readers are streaming via `NextSeq()`.
- **`align/`** — Smith-Waterman local alignment with affine gap penalties. Includes special handling for Oxford Nanopore homopolymer error profiles, plus MSA via incremental consensus.
- **`htsio/`** — SAM/BAM/CRAM reading and writing. Native BAM and SAM readers/writers; samtools only for CRAM. Includes BAI/TBI/CSI index parsers, tabix reader/writer, sorted BAM writer with merge sort. Subpackages: `bam`, `bgzf`, `cram`, `codec`, `sam`, `tabix`, `bbi`.
- **`htsio/bgzf/`** — BGZF (Blocked GNU Zip Format) reader, writer, and indexed reader with LRU block cache. Used by BAM and tabix layers.
- **`htsio/bbi/`** — Random-access reader for UCSC BBI files (bigWig/bigBed); self-indexed, standard-library only. Reader surface mirrors `htsio/tabix`.
- **`bed/`** — Streaming and tabix-indexed BED readers plus a writer. Core type `BedRecord` (0-based half-open, BED6 fields + verbatim `Extras`).
- **`gtf/`** — GTF gene-model parsing (genes → transcripts → exons/CDS) with an interval index and genic-region classification; a port of ngsutilsj's `GTFAnnotationSource`. `AnnotationSource` (in-memory) and `IndexedAnnotationSource` (tabix-backed).
- **`vcf/`** — Streaming and tabix-indexed VCF reader/writer with a lazy record model (`VcfRecord` parses columns on first access; `Pos` is 1-based). Subpackages: `vcf/annotate` (composable INFO/FORMAT annotators + `Pipeline`) and `vcf/filter` (composable FILTER-stamping filters).
- **`iosource/`** — Pluggable random-access byte sources behind a concurrency-safe `io.ReaderAt`: local-file and HTTP(S)-Range implementations, plus sibling-index resolution. Lets index-driven readers fetch byte ranges from remote files. Wired into `htsio/tabix` (`NewReaderFromSource`), `htsio/bbi` (`NewReaderFromSource`) and `varstore` (`OpenParquetContext`); use `iosource.ReadSeeker` to adapt a `ByteSource` for BGZF. `iosource.Open(ctx, locator)` dispatches on scheme; **`iosource/s3`** registers `s3://` (blank import) and is the reason cghts depends on the AWS SDK.
- **`varstore/`** — Sparse alt-only genotype store in Parquet (`calls`/`sites`/`regions` members) plus a VCF-backed store behind the same `Store` interface. The point of the interface is that both backends answer identically: a query against a VCF and against a store converted from it must agree, which is what the consuming CLI's cross-backend tests assert. **One method does the querying** — `Calls(Query) (iter.Seq2[Call, error], error)`. `Query` selects sites (`Loci`, `Spans`) on one axis and `Samples` on the other, and an **empty selector means no restriction** on that axis, so the zero `Query` is every genotype in the store. `IncludeRef` additionally emits reference (`0/0`) calls. Implementations pick their own access strategy from the query's shape — pruned per-locus lookups for a few loci, one ordered pass for many — so a caller cannot turn a variant panel into one row-group decode per variant. `CollectCalls` buffers into a slice for callers that do not stream.
  Rules callers must respect: rows arrive in the **store's own order** (contig order as written, *not* lexicographic — sorting `chrom` as a string puts `chr10` before `chr2` and yields an unindexable VCF); `Missing == -1` (columns are required, absence is in-band, so test before use); `SpansSites` (a locus absent from `sites.parquet` is unassayed for *every* sample, and the region runs must not be consulted); and `CanonKey`/`SameChrom` (stores keep the source's own contig spelling, so `chr17`/`17`/`NC_000017.11` must be compared canonically). A reference call is emitted only where the genotype was **all**-reference: at a multiallelic record a 0/2 sample carries a different alternate, so writing `0/0` for it would fabricate a genotype — the exclusion is keyed on the source *record* (`Locus.Record()`), not the split locus. `Calls` with `IncludeRef` and `Classify` both return `ErrNotClassifiable` rather than degrade when the store lacks the evidence (no `sites`/`regions`, or built `--no-callable`).
- **`support/sequtils/`** — DNA utilities: IUPAC ambiguity code matching, reverse complement, homopolymer run analysis, 4-bit DNA encoding.
- **`support/stats/`** — 2×2 Fisher exact test, Phred/log2 conversions.
- **`support/utils/`** — General utilities: semaphore for concurrency, float formatting, position-tracking reader.
- **`support/stringutils/`** — String helpers.
- **`analysis/seq/`** — Sequence analysis (GC content); package `seqanalysis`.

### HTS I/O System

The `htsio/` package provides native SAM/BAM I/O without external dependencies (samtools only for CRAM):

- `SamReader` interface: `Next()`, `Header()`, `Query()`, `Close()`
- `Query()` returns `iter.Seq2[*SamRecord, error]` — uses Go 1.23 range-over-func
- `NewSamReader()` auto-detects: `.bam` → `BamReader`, `.sam`/`.sam.gz` → `SamTextReader`, `.cram` → `SamtoolsSamReader`
- `NewSamWriter()` auto-selects: BAM (sorted/unsorted) → native, CRAM → samtools
- All query coordinates are 0-based half-open
- `ParseRegion()` converts samtools-style strings to 0-based half-open
- `IterReader()` bridges `iter.Seq2` back to `SamReader` for legacy callers
- `TabixReader`/`TabixWriter` handle tabix-indexed text files (BED, VCF, GFF) with TBI or CSI indexes
- `bgzf.IndexedReader` has an LRU block cache shared by BAI and tabix query paths

### Alignment System

The aligner (`align/`) is the most complex component:

- `NewLocalAligner()` — Smith-Waterman with soft clipping (for partial matches)
- `NewGlobalAligner()` — Full-sequence alignment
- `DnaAlignmentDefaults()` — Presets for Illumina short reads
- `OntAlignmentDefaults()` — Presets for Oxford Nanopore (looser gap penalties, homopolymer discounts)
- `AlignBatch()` — Parallel alignment using a semaphore-controlled goroutine pool
- Homopolymer discounts are precalculated and cached for performance

CIGAR strings use standard ops: M (match), I (insertion), D (deletion), S (soft clip). Helper functions `CigarCondense`/`CigarExpand` convert between run-length encoded and per-base forms.

## Note

This library carries **no CLI dependencies** (no cobra/pflag) — CLI concerns belong
in `cgkit` or `nupa`. Keep it that way.

Third-party dependencies are deliberately few, and each one must be a *file-format
or codec* need that the standard library cannot meet:

- `github.com/ulikunitz/xz` — CRAM LZMA.
- `github.com/parquet-go/parquet-go` — the `varstore/` Parquet genotype store.
  This is the heaviest dependency here: it pulls ten modules transitively,
  including `google.golang.org/protobuf` and `twpayne/go-geom`, and it is why the
  `go` directive is 1.24. It was added knowingly, because `varstore` is consumed by
  both `cgkit` and Cohort Studio and its four-state classification semantics must
  not be reimplemented per consumer. Go links only imported packages, so binaries
  that never touch `varstore` are unaffected in size — the cost is `go.sum` weight
  and supply-chain surface for every consumer.

Do not add anything further without a comparable justification.
