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
  The tabix index records its **preset** in the header `format` field, and the VCF preset derives each record's end **from the record** rather than from a column: `col_end = 0` is not a claim that VCF records are one base wide, which is a separate spec rule applying when `col_beg == col_end`. Writer and reader must both apply it -- the reader's end feeds the overlap test, so an understated record is dropped even when the index pointed straight at it. Before this the preset was hardcoded to 0, so every VCF `.tbi` written here declared itself generic and nothing could apply VCF rules on read; indexes written earlier still read that way, so a gVCF or SV VCF must be re-indexed to benefit.
- **`htsio/bgzf/`** — BGZF (Blocked GNU Zip Format) reader, writer, and indexed reader with LRU block cache. Used by BAM and tabix layers.
- **`htsio/bbi/`** — Random-access reader for UCSC BBI files (bigWig/bigBed); self-indexed, standard-library only. Reader surface mirrors `htsio/tabix`.
- **`bed/`** — Streaming and tabix-indexed BED readers plus a writer. Core type `BedRecord` (0-based half-open, BED6 fields + verbatim `Extras`).
- **`gtf/`** — GTF gene-model parsing (genes → transcripts → exons/CDS) with an interval index and genic-region classification; a port of ngsutilsj's `GTFAnnotationSource`. `AnnotationSource` (in-memory) and `IndexedAnnotationSource` (tabix-backed).
- **`vcf/`** — Streaming and tabix-indexed VCF reader/writer with a lazy record model (`VcfRecord` parses columns on first access; `Pos` is 1-based). Subpackages: `vcf/annotate` (composable INFO/FORMAT annotators + `Pipeline`) and `vcf/filter` (composable FILTER-stamping filters).
  `RefSpan`/`RefSpanEnd` give a record's reference extent -- the widest of `len(REF)`, `INFO/END`, `INFO/SVLEN` for symbolic alternates measured on the reference, and `FORMAT/LEN` for reference blocks -- and `IsRefBlock`/`IsRefBlockAlt` tell a gVCF block from a variant. Do **not** reach for `AltPositions` for this: it resolves an SV's *partner breakpoint*, which may be on another chromosome and may precede POS, and using it as an extent gave `vcf-tobed` a one-base overhang on every deletion. The rules live once in `internal/vcfspan`, shared with `htsio/tabix` (which cannot import `vcf`, since `vcf` imports it) so an index and a reader cannot disagree about how far a record reaches. `IsRefBlock` means *no* alternate is real, because GATK writes the block allele beside a genuine one as `G,<NON_REF>` and treating its presence as sufficient drops every real call in a gVCF. `BlockDepth` reads `FORMAT/MIN_DP` (falling back to `DP`, which is a **weaker** claim -- the depth at POS, not across the block) and `BlockRGQ` reads `RGQ`, which GATK writes in place of `GQ` on a block. Absent, empty and `.` are all unknown, never zero: a zero depth is a real claim and unknown depth is not, and conflating them makes a gate admit everything.
- **`iosource/`** — Pluggable random-access byte sources behind a concurrency-safe `io.ReaderAt`: local-file and HTTP(S)-Range implementations, plus sibling-index resolution. Lets index-driven readers fetch byte ranges from remote files. Wired into `htsio/tabix` (`NewReaderFromSource`), `htsio/bbi` (`NewReaderFromSource`) and `varstore` (`OpenParquetContext`); use `iosource.ReadSeeker` to adapt a `ByteSource` for BGZF. `iosource.Open(ctx, locator)` dispatches on scheme; **`iosource/s3`** registers `s3://` (blank import) and is the reason cghts depends on the AWS SDK.
- **`varstore/`** — Sparse alt-only genotype store in Parquet plus a VCF-backed store behind the same `Store` interface. `OpenStore(ctx, locator, kind)`/`StoreKind` own backend resolution (a VCF suffix, else a manifest found at the locator) so the CLI does not re-derive it; `StoreKind` reports an unregistered scheme as a *transport* problem rather than an unrecognized format, since the fix is importing a package, not checking the path. `VcfStore` holds its ctx (`scan` opens lazily on every call and the interface has no ctx-taking method), probes for a tabix index **once at open** rather than per locus — `SiteKnown` runs once per queried variant, so remotely the old per-locus `os.Stat(path+".tbi")` was one index re-fetch per lookup, and could not work for a URL at all — and exposes `Indexed()` so a caller can warn that an unindexed remote VCF streams the whole object per query. `Sites(fn)` is on the interface and implemented by both: for Parquet it reads the catalog, for a VCF it is a full pass computing AC/AN/NCarriers on the way, leaving `NCalled`/`NLowDP` zero because both are defined against a `--min-dp` a plain VCF has no record of. **A store is a directory** holding `calls.parquet`/`sites.parquet`/`regions.parquet`; `cohort`, `cohort/` and `cohort/calls.parquet` are all accepted and `TrimStoreSuffix` normalizes them without touching the filesystem (it cannot, for a remote locator). The former filename-prefix form (`cohort.calls.parquet`) is gone: keeping both meant every path decision first had to guess which was meant, and `MemberPath` uses `joinStore` rather than `filepath.Join` because Join would clean `s3://bucket/x` into `s3:/bucket/x`. `NewWriter` creates the directory itself, since creating a store *is* creating its directory.
  **A store is readable only with its `manifest.json.gz`**, written last by `Writer.Finish` (gzipped JSON, temp+rename, `Sync` before the rename so the marker cannot become visible ahead of what it certifies). `NewWriter` deletes any existing manifest *first*, so a store under construction never carries a completion marker — otherwise a `--force` retry that died would leave the previous run's manifest vouching for this run's half-written members. Open verifies `complete`, refuses a `format_version` from the future, and checks each member's recorded row count against that member's own parquet footer; that last check is the only thing that catches a well-formed member belonging to a *different* conversion, since sites and regions carry no metadata of their own. The per-chromosome census is accumulated by the writer in `WriteCall`/`WriteSite` rather than supplied by the caller, so it records rows rather than intent — which matters because `MetaSource`/`MetaContigs` are stamped at construction and a three-chromosome store from a 22-input run names all 22 of both. A member may be absent only if the manifest recorded zero rows in it (the `--no-callable` shape). No escape hatch: an unfinished store answers "not assayed" for everything it never reached, which is indistinguishable from the honest answer.
  **A conversion that cannot finish must not leave something that looks finished.** A parquet footer is written only by the writer's `Close`, so a valid footer means a member was completed — which is why `Discard` unlinks *without* closing the writers (closing them first wrote complete footers into files about to be removed, so a kill mid-`Discard`, or a failed unlink, left a well-formed partial store), `abort` reports the removal error instead of swallowing it, and `Close` stops at the first failure rather than finalizing the rest. Three well-formed members that disagree is the one outcome a reader cannot detect. `sites`/`regions` footers are also parsed at open, so a truncated member is an error while a genuinely absent one (`--no-callable`) is not. The point of the interface is that both backends answer identically: a query against a VCF and against a store converted from it must agree, which is what the consuming CLI's cross-backend tests assert. **One method does the querying** — `Calls(Query) (iter.Seq2[Call, error], error)`. `Query` selects sites (`Loci`, `Spans`) on one axis and `Samples` on the other, and an **empty selector means no restriction** on that axis, so the zero `Query` is every genotype in the store. `IncludeRef` additionally emits reference (`0/0`) calls. `CollectCalls` buffers into a slice for callers that do not stream.
  **A query is always one pass, however many selectors it names**; its shape changes only how much of that pass is skipped — the Parquet side prunes row groups against the *union* of the selectors, the VCF side takes a tabix seek when there is exactly one selector and otherwise scans and filters per row. So a caller cannot turn an N-variant panel into N lookups, which is the whole reason the four earlier query methods collapsed into this one: **bulk is flat where a per-locus loop is linear** (~19 ms vs 6.6 s at 1000 targets, 349×, on the machine that ran it). Treat the ratio as the durable claim and the milliseconds as hardware-bound; `cgkit`'s `internal/cmd/vcfcmd/bench_test.go` regenerates the corpus from a seed and re-measures. Do not add a per-locus fast path without re-running it.
  Rules callers must respect: rows arrive in the **store's own order** (contig order as written, *not* lexicographic — sorting `chrom` as a string puts `chr10` before `chr2` and yields an unindexable VCF); `Missing == -1` (columns are required, absence is in-band, so test before use); `SpansSites` (a locus absent from `sites.parquet` is unassayed for *every* sample, and the region runs must not be consulted); and `CanonKey`/`SameChrom` (stores keep the source's own contig spelling, so `chr17`/`17`/`NC_000017.11` must be compared canonically). A reference call is emitted only where the genotype was **all**-reference: at a multiallelic record a 0/2 sample carries a different alternate, so writing `0/0` for it would fabricate a genotype — the exclusion is keyed on the source *record* (`Locus.Record()`), not the split locus. `Calls` with `IncludeRef` and `Classify` both return `ErrNotClassifiable` rather than degrade when the store lacks the evidence (no `sites`/`regions`, or built `--no-callable`).
  One asymmetry the interface does *not* hide: the backends promise identical **answers**, not identical memory. `IncludeRef` on the Parquet side has to interleave rows derived from three files in one order, so it materializes before it yields — bounded by the query, not by the store, which is fine for a panel or a named sample but means an *unrestricted* `IncludeRef` query holds the whole dense matrix. (The VCF side streams, since a record already carries every sample.) A streaming merge of the sites and calls scans is the fix if that case becomes real. **A span selects by overlap** (`Span.Overlaps`), not by start position: a record reaching into a region from before it must match, or a region query silently loses every long deletion, SV and gVCF block. The extent lives in the `ref_end` column on `Call` and `Site`; the `Writer` derives it from `len(REF)` when a caller leaves it zero, so no writer has to opt in and the column is never a field of zeros that would make `HasRefSpans` lie. Two consequences to keep in mind when touching this: `wantsSite` and `wantsRecord` must share one span test, since `wantsSite ⊆ wantsRecord` is a structural invariant and widening only the former fabricates `0/0` for a sample carrying a sibling ALT; and the pruning lower bound cannot come from `pos`, because a matching record may start below a row group's minimum — it comes from `max(ref_end)`, which is unsorted and so bounds only that side, the same asymmetry `spanRunFilter` documents for the regions file. **`VcfStore` reads gVCF reference blocks as coverage.** A block yields one reference row per sample -- never one per base, however wide -- with `Alt` `.`, `RefEnd` carrying the extent, `DP` left `Missing` because a block measures no single base, and `MinDP` carrying the floor. The recorded genotype survives, so a `./.` block is skipped rather than reported `0/0`. `Classify`/`SiteKnown` treat block coverage as an observation, which an explicit record at the locus outranks. `Call.MinDP` is the tightest lower bound the source vouches for and is tagged `parquet:"-"`: derived per query today, a real column when a blocks store exists. **`Gate.Admits` prefers it over `DP`, and that is load-bearing** -- `Missing` passes every gate by design, so a block-derived call would otherwise be admitted unconditionally, the gate silently doing nothing on exactly the state gVCF makes trustworthy. `HasRefSpans()` reads the schema rather than a metadata key, since the column either exists or it does not where a key could disagree with the data. `MetaContigs`/`ParquetStore.Contigs()` carry the source's `##contig` lines verbatim, because a store is expected to be exported back to VCF and nothing else in the store records which reference it was called against — lengths especially cannot be recovered from the loci a query happens to touch.
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
