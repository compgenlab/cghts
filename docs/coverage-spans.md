# Coverage spans — answering off the catalog

## The problem

A varstore is built from a joint-called cohort VCF, and its `regions` table records
runs of **catalog sites** a sample was called at. `SpanSemantics` says so outright:
the interval marks catalog sites, and *nothing is claimed about the bases in between*.

That is the honest reading of a pVCF, and it has two consequences that only appear
later.

**Across batches.** Two callsets of 10,000 samples each, joint-called separately,
have different catalogs. A variant carried by three people in batch 1 is simply
absent from batch 2's catalog, because nobody there carried it. Ask "who carries
X?" and batch 1 answers with 3 carriers and 9,997 reference — while all 10,000 of
batch 2 come back `NotAssayed`, because their pVCF never proved they were reference
there. The carrier rate then divides by 10,000 instead of 20,000.

Nothing is *wrong*: `Roll` drops them (`rollup.go`, `Called == 0`), the column has
no value for them, and Cohort Studio reports them as not eligible. It is correct
and it is severely pessimistic, and it compounds — with five batches, a variant
seen in one leaves 80% of the cohort unevaluable.

**Across time**, which is worse. The catalog is frozen at conversion. A variant
discovered next year — a new batch, a new ClinVar release, a clinical question
about a site nobody in the first cohort happened to carry — is permanently
unanswerable for every sample already converted. Re-converting recovers nothing,
because the pVCF never held the information. **The only moment it can be captured
is at conversion, from a source that still has it.**

## What already exists

More than it first appears. `CalledSiteRun` is `{sample, chrom, start, end,
n_sites, min_dp}` — an interval carrying a depth bound — and the converter already
breaks a run when depth leaves its band, which is structurally what a gVCF does
when it breaks a block on GQ. The rows would look the same.

What differs is the **claim**, not the shape:

| | claims |
|---|---|
| `SpansSites` | this sample was called at the catalog sites inside this interval |
| `SpansBlocks` | every base in this interval was covered at depth ≥ D |

That is why `SpanSemantics` is a manifest field rather than a schema difference,
and why `SpansBlocks` has always been defined and never produced.

## The design

**A new optional table, `coverage.parquet`**, holding genomic block spans, beside
the existing `regions.parquet` holding catalog-site runs.

Two tables rather than one table with a kind column, for three reasons:

- **They arrive at different times from different sources.** Catalog runs come from
  the pVCF during conversion. Blocks come from per-sample gVCFs, which may be
  processed separately or later. A separate table can be attached to a volume
  without rewriting anything.
- **The common query must not pay for the rare one.** A catalog-site lookup is the
  hot path; making it scan or skip block rows to find site rows taxes every query
  for a capability most of them never use.
- **Absent must mean unknown.** A volume converted before this exists has no
  coverage table, which is *"nobody said"* — not *"covered nowhere"*. An optional
  table gives that for free, the way `openOptionalShardSet` already does.

### Schema

```go
// CoverageBlock is a positive claim about every base in a span.
type CoverageBlock struct {
    SampleID string `parquet:"sample_id,dict"`
    Chrom    string `parquet:"chrom,dict"`
    Start    int32  `parquet:"start"`
    End      int32  `parquet:"end"`

    // MinDP is the floor the source vouches for across the whole block --
    // a gVCF's MIN_DP. This is the field a Ref inferred here rests on.
    MinDP int32 `parquet:"min_dp"`

    // GQ as the source stated it for this block. Recorded because it cannot be
    // recovered later and costs nothing here, NOT because anything gates on it:
    // GQ is not comparable across callers, saturates at 99, and would make a
    // gate mean different things in two parts of one release. See "GQ vs DP".
    GQ int32 `parquet:"gq"`
}
```

No `n_sites`: this table is not about catalog sites.

### Do not re-band

Take the gVCF's own block boundaries as given.

GATK already bands reference blocks by GQ, with the ranges declared in the
`##GVCFBlock` header, so the boundaries are chosen upstream and a converter can
only ever coarsen them. Passing them through means `MinDP` and `GQ` are both
exactly what the source said, no banding parameter has to be defended, and the
conversion is a translation rather than a re-interpretation.

If the row count proves unaffordable, `--coverage-bands` can merge adjacent blocks
within a depth class afterwards — the same dial `--depth-bands` already is, applied
as an option rather than a default. **Do not pick that default before measuring**;
the depth-band fixture swung an order of magnitude on binning choice alone.

### Precedence — the load-bearing rule

> **At a catalog site, `regions` is authoritative. `coverage` answers only off
> catalog.**

Get this backwards and every declined genotype silently becomes a reference call.
In a joint-called pVCF every sample has a genotype at every record, so *"in the
catalog but not in a run"* means the caller **declined** — allele-balance
rejection, a local assembly failure, a mapping problem. A coverage block spanning
that position must not overrule it. Those disagreements cluster exactly where
calling is hard, which is where a fabricated Ref does the most damage.

This is also what `SpansBlocks` has always said: *"Only such a store may answer for
positions absent from the catalog."*

### What becomes answerable

`Classify` at a locus outside the sites catalog. Today it is `NotAssayed` for
everyone, unconditionally. With coverage present: a sample whose block covers the
position at `MinDP ≥ gate` is **Ref**; otherwise still `NotAssayed`.

That is the whole point — batch 2's samples become provable reference at batch 1's
variant, and the denominator is the cohort again.

### A block-derived Ref is weaker, and must say so

*"Reads were here at depth 30"* is not *"the caller would have called this person
reference."* A Ref inferred from a block is a weaker claim than one inferred from a
catalog run, and the difference matters precisely where the two would disagree.

So the state must carry its provenance, and Cohort Studio's `wins()` must prefer a
catalog-run Ref over a block-derived one when two parts observe the same person.
Not a fifth state — the four states stay as they are — but a field on the returned
observation. Without it, "positive evidence dominates" runs in reverse: a
block-derived Ref from a clean region beats a real no-call from a hard one.

### Guards

- **Rolling up a part with coverage and one without is legal and must be stated.**
  The part without contributes `NotAssayed` off catalog, which is correct; the
  coverage note has to say that a source could not answer rather than letting the
  asymmetry pass silently. Same shape as the existing assay and depth-band guards.
- **A volume's coverage table is verified against its own roster at open**, exactly
  as calls and regions are. A coverage table from another conversion is the same
  class of error `TableInfo.Rows` already exists to catch.

## Open, and needs a measurement

**Row count per sample is the one number that could change the design**, and it is
unmeasured. Estimates so far, all the same order of magnitude:

| | intervals per sample |
|---|---|
| today's banded catalog runs | ~1.4M *(measured, synthetic fixture)* |
| GIAB high-confidence BED, NA12878 | 649k *(measured)* |
| GATK gVCF blocks | ~1–5M *(literature, unmeasured)* |

At the ~8 bytes/row measured for the regions table, 1–5M blocks is 8–40 MB per
sample, so 80–400 GB across 10,000. The bottom of that range is comfortable and
the top is not, and the difference is entirely how finely the source banded.

**Measure a real gVCF before fixing any default.** One chromosome of one sample
answers it.

## Not in scope

- Reconstructing a gVCF from a varstore. Coverage spans make reference calls
  provable; they do not make the store lossless, and `gvcf_info` stays discarded.
- Merging batches on the sample axis. Still the open architectural question, and
  still unaffected by this.
- Per-base depth. A block carries one floor, as a gVCF block does.
