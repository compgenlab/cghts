// Package varstore holds the on-disk schema for a sparse genotype store and a
// uniform way to query genotypes regardless of which format backs them.
//
// A Parquet store is a directory. All three tables are required: the calls file
// alone cannot distinguish a confidently-called reference genotype from a
// position that was never assayed, which is the distinction Classify exists to
// make.
//
//	cohort/
//	  calls.parquet     one row per ALT-carrying genotype
//	  sites.parquet     one row per interrogated site, sample-independent
//	  regions.parquet   runs of catalog sites at which a sample was called
//
// The tables are only meaningful together, so the store is one thing to copy,
// move and delete. A trailing separator on the base is optional, and a table
// path is accepted anywhere a store name is: all of "cohort", "cohort/" and
// "cohort/calls.parquet" name the same store.
//
// # What a store built from a plain VCF can answer
//
// Only the variants the VCF actually contains. A plain VCF asserts nothing
// about the positions between its records: an absent base was not observed to
// be reference, it simply was not reported. So the sites catalog is the exact
// boundary of what is knowable, and a query for any locus outside it yields
// StateNotAssayed for every sample -- not StateNonCarrier.
//
// This is why the run intervals in the regions file must never be read as
// coverage. They compress a per-sample, per-site record of "this sample was
// called at these catalog sites"; the gaps between those sites are not part of
// the claim. A gVCF is different, because its reference blocks (END, MIN_DP)
// are positive statements about spans, and only such a store could answer
// off-catalog positions.
//
// The sites file cannot be reconstructed from the calls file. Taking the
// distinct loci out of the calls recovers every site only when the store holds
// the entire joint callset; over a subset of samples the sites where nobody in
// that subset carries an ALT vanish silently, and a query would then report
// "never interrogated" for a position that was in fact interrogated.
package varstore

import (
	"fmt"
	"os"
	"strings"
)

// Missing marks an absent integer field (DP, GQ, AD) in a Parquet row. VCFs
// routinely omit these -- a GT-only phased panel has no DP at all -- and the
// columns are kept non-optional so reads stay a flat scan, so the absence has
// to be in-band. Callers must test against Missing before using a value; a
// naive comparison would treat it as an extremely low quality score.
const Missing int32 = -1

// Call is one ALT-carrying genotype for one sample at one biallelic site.
// Records are split so that every row carries exactly one ALT allele.
type Call struct {
	SampleID string `parquet:"sample_id,dict"`
	Chrom    string `parquet:"chrom,dict"`
	Pos      int32  `parquet:"pos"`
	Ref      string `parquet:"ref,dict"`
	Alt      string `parquet:"alt,dict"`
	// RefEnd is the 0-based exclusive end of the reference bases the source
	// record covered, which for a long deletion, a symbolic ALT or a gVCF block
	// is past Pos. Zero in stores written before this column existed; treat that
	// as "unknown" and fall back to Pos, never as a zero-length record.
	RefEnd int32  `parquet:"ref_end"`
	GT     string `parquet:"gt,dict"`
	DP     int32  `parquet:"dp"`
	ADRef  int32  `parquet:"ad_ref"`
	ADAlt  int32  `parquet:"ad_alt"`
	GQ     int32  `parquet:"gq"`

	// MinDP is the tightest lower bound on depth the source vouches for here, where
	// DP is the depth actually recorded at this position. They come apart for a gVCF
	// reference block: MIN_DP is the floor across the whole block and there is no
	// per-base depth at all, so DP stays Missing while MinDP carries the claim.
	// Zero means unknown.
	//
	// Deliberately not a Parquet column. Today it is derived per query by the VCF
	// backend and a store has nowhere to have recorded it; when a gVCF-derived blocks
	// store exists this becomes a real column. Note that an absent column reads back
	// as 0, which is why unknown is 0 here rather than Missing -- the two conventions
	// would otherwise disagree the moment the column appears.
	MinDP int32 `parquet:"-"`
}

// Locus returns the site this call belongs to.
func (c Call) Locus() Locus {
	return Locus{Chrom: c.Chrom, Pos: c.Pos, Ref: c.Ref, Alt: c.Alt}
}

// HomRefGT is the genotype a reconstructed reference call reports.
const HomRefGT = "0/0"

// HomRefCall synthesizes the row for a reference genotype recovered from a
// Parquet store.
//
// A store keeps only ALT-carrying genotypes, so the reference call itself was
// never written down; what survives conversion is the *fact* that the sample was
// called at this catalog site, held in the regions file. So the genotype string
// is synthesized and every quality field is Missing. The row asserts exactly
// "observed here, and observed to be reference" and nothing more -- notably not
// the ploidy or phasing of the original genotype, which are unrecoverable.
//
// A VCF-backed store does not go through this: it still has the real genotype
// and its DP/AD/GQ, and reports them.
func HomRefCall(sample string, l Locus) Call {
	return Call{
		SampleID: sample,
		Chrom:    l.Chrom,
		Pos:      l.Pos,
		Ref:      l.Ref,
		Alt:      l.Alt,
		GT:       HomRefGT,
		DP:       Missing,
		ADRef:    Missing,
		ADAlt:    Missing,
		GQ:       Missing,
	}
}

// Site is one interrogated variant site, independent of any sample. Counts are
// taken across every sample present in the source, so a site with AC == 0 still
// records that the position was examined.
//
// Allele counts and sample counts are deliberately both present, because they
// answer different questions and neither can be derived from the other:
//
//   - AC / AN are ALLELE counts, per the VCF convention. A 1/1 genotype
//     contributes 2 to each. They are computed from GT alone and are not
//     depth-gated, so AF is exactly AC/AN.
//   - NCarriers / NCalled / NLowDP are SAMPLE counts. A 1/1 genotype is one
//     carrier, and NCalled/NLowDP additionally reflect the --min-dp threshold
//     used at conversion.
//
// So AC >= NCarriers whenever any homozygote is present, and AN is unrelated to
// NCalled both in unit and in gating. Counts are over the samples in this
// store, not copied from the source's INFO fields, which would be wrong after
// splitting multiallelics or converting a subset of a cohort.
type Site struct {
	Chrom     string `parquet:"chrom,dict"`
	Pos       int32  `parquet:"pos"`
	Ref       string `parquet:"ref,dict"`
	Alt       string `parquet:"alt,dict"`
	RefEnd    int32  `parquet:"ref_end"`    // 0-based exclusive; see Call.RefEnd
	AC        int32  `parquet:"ac"`         // alt alleles observed, this ALT
	AN        int32  `parquet:"an"`         // called alleles at the site
	NCarriers int32  `parquet:"n_carriers"` // samples with >=1 copy of this ALT
	NLowDP    int32  `parquet:"n_lowdp"`    // samples failing --min-dp here
	NCalled   int32  `parquet:"n_called"`   // samples called and passing --min-dp
}

// AF returns the alternate allele frequency, or 0 when nothing was called.
func (s Site) AF() float64 {
	if s.AN == 0 {
		return 0
	}
	return float64(s.AC) / float64(s.AN)
}

// Locus returns this site's identity.
func (s Site) Locus() Locus {
	return Locus{Chrom: s.Chrom, Pos: s.Pos, Ref: s.Ref, Alt: s.Alt}
}

// CalledSiteRun records that one sample was successfully called, at adequate
// depth, at every catalog site in [Start, End]. Start and End are the first and
// last such site positions, and NSites is how many catalog sites the run covers.
//
// IT SAYS NOTHING ABOUT THE GAPS BETWEEN THOSE SITES. A plain VCF records
// nothing between its variant records, so there is no basis for calling an
// intervening base reference -- the caller may simply never have looked. The
// interval form is a COMPRESSION of a per-sample, per-site fact ("called here,
// and here, and here"), not a statement about genomic territory. Reading it as
// a coverage interval would silently manufacture reference observations for
// positions that were never interrogated, which is the precise error the
// four-way classification exists to prevent.
//
// Consequently a run is only meaningful at positions that appear in the sites
// catalog, and Classify checks the catalog first for exactly that reason. A
// gVCF, whose reference blocks carry END and MIN_DP, is what would license
// answering off-catalog positions; see SpanSemantics.
type CalledSiteRun struct {
	SampleID string `parquet:"sample_id,dict"`
	Chrom    string `parquet:"chrom,dict"`
	Start    int32  `parquet:"start"`
	End      int32  `parquet:"end"`
	NSites   int32  `parquet:"n_sites"`

	// MinDP is the lowest depth seen anywhere inside the run: the tightest
	// bound the source vouches for across the whole span.
	//
	// WITHOUT IT A REFERENCE CALL CARRIES NO CONFIDENCE. `0/0` is never stored,
	// so a Ref is inferred from being inside a run -- and a run alone says only
	// "at or above the conversion gate", which makes a reference at 60x and one
	// at exactly 10x indistinguishable forever. That is fine until two sources
	// disagree about a person, at which point the reference side of the
	// argument has nothing to say.
	//
	// Zero in stores written before this column existed. Treat that as unknown
	// and fall back to the manifest's MinDP, never as a depth of zero -- the
	// same rule the schema already carries for RefEnd.
	MinDP int32 `parquet:"min_dp"`
}

// CoverageBlock is a positive claim about every base in a span.
//
// WHAT IT IS FOR. regions records runs of CATALOG SITES, which is all a
// joint-called pVCF supports: it says a sample was called at the sites inside
// an interval and claims nothing about the bases between them. That is honest
// and it is pessimistic in two ways that only surface later.
//
// Across batches, two callsets joint-called separately have different catalogs,
// so a variant carried by three people in one is simply absent from the other's
// catalog -- and every sample in that second batch reads as NotAssayed rather
// than reference, halving a denominator. Across fifteen batches it leaves
// fourteen fifteenths of a cohort unevaluable at any batch-specific variant.
//
// Across time it is worse. The catalog is frozen at conversion, so a variant
// discovered next year is permanently unanswerable for everything already
// converted, and reconverting recovers nothing because the pVCF never held it.
// A coverage block is the only thing that can answer it, and the only moment to
// capture one is at conversion from a source that still has it -- a gVCF, or a
// per-sample callable BED.
//
// THE CLAIM IS BOUNDED BY MaxGap. A block built by merging across small
// uncovered stretches does not mean "every base"; it means "covered at MinDP or
// above, with no uncovered run longer than the manifest's MaxGap". That
// tolerance is what makes the table affordable -- merging a real per-sample mask
// across 10bp collapsed it 8.4x for 0.35% over-claim -- and it is recorded in
// ManifestParams for the same reason MinDP is: two stores with different
// tolerances do not mean the same thing by "covered".
type CoverageBlock struct {
	SampleID string `parquet:"sample_id,dict"`
	Chrom    string `parquet:"chrom,dict"`

	// Start and End bound the covered span, inclusive, in the source's own
	// coordinates -- the same convention CalledSiteRun uses.
	Start int32 `parquet:"start"`
	End   int32 `parquet:"end"`

	// MinDP is the floor the source vouches for across the whole span, a gVCF's
	// MIN_DP. A reference call inferred from this block rests on it, so it is
	// the field the gate compares against, exactly as for a callable run.
	MinDP int32 `parquet:"min_dp"`

	// GQ as the source stated it, or Missing.
	//
	// RECORDED, NOT GATED ON. It cannot be recovered later and costs nothing
	// here, which is the same reasoning that captures INFO fields. But nothing
	// should gate on it: GQ is not comparable across callers, saturates at 99,
	// and a gate resting on it would mean different things in two parts of one
	// release -- which is the failure the depth gate is arranged to avoid.
	GQ int32 `parquet:"gq"`
}

// CoveragePath returns the file holding a store's coverage blocks.
func CoveragePath(base string) string { return TablePath(base, CoverageTable) }

// The tables of a store. A store IS a directory, and they sit inside it under
// fixed names:
//
//	cohort/
//	  calls.parquet
//	  sites.parquet
//	  regions.parquet
//
// A trailing separator on the base is optional and means nothing;
// "cohort" and "cohort/" name the same store.
//
// There used to be a second, filename-prefix form ("cohort.calls.parquet"), and
// dropping it removed more than it cost. The tables are only meaningful
// together, so a store wants to be one thing to copy, move and delete; a fixed
// directory is also how every table format built on Parquet ships, so external
// readers can be pointed straight at a table. And path resolution collapses:
// there is no longer a base that might be a prefix or might be a directory, no
// separator-sensitive TablePath, and no need to guess which was meant.
const (
	CallsTable   = "calls"
	SitesTable   = "sites"
	RegionsTable = "regions"

	// CoverageTable is OPTIONAL, and its absence means "nobody said" rather
	// than "covered nowhere". See CoverageBlock.
	CoverageTable = "coverage"
)

// TableFile returns the file name one table of a store is stored under,
// without a directory. What a Sink is addressed by: where the store lives is
// the sink's business, and the table is the same name wherever it goes.
func TableFile(table string) string { return table + ".parquet" }

// TablePath returns the file holding one table of the store at base.
func TablePath(base, table string) string {
	return joinStore(base, table+".parquet")
}

// joinStore appends a name to a store directory.
//
// It is deliberately not filepath.Join, which cleans its result and would turn
// "s3://bucket/cohort" into "s3:/bucket/cohort" by collapsing the scheme's
// double slash. A forward slash is used as the separator on every platform,
// since Go accepts it in filesystem paths on Windows too and it is the only
// thing a locator can use.
func joinStore(base, name string) string {
	switch {
	case base == "":
		return name
	case strings.HasSuffix(base, "/"), strings.HasSuffix(base, string(os.PathSeparator)):
		return base + name
	default:
		return base + "/" + name
	}
}

// trimStoreDir removes a trailing separator from a store directory, so that the
// base is spelled one way internally whatever the caller typed. The filesystem
// root is left alone: "/" is a directory name, not a decorated empty one.
func trimStoreDir(base string) string {
	for len(base) > 1 {
		last := base[len(base)-1]
		if last != '/' && last != os.PathSeparator {
			break
		}
		base = base[:len(base)-1]
	}
	return base
}

// CallsPath returns the calls file for a store base name.
func CallsPath(base string) string { return TablePath(base, CallsTable) }

// SitesPath returns the sites file for a store base name.
func SitesPath(base string) string { return TablePath(base, SitesTable) }

// RegionsPath returns the callable-regions file for a store base name.
func RegionsPath(base string) string { return TablePath(base, RegionsTable) }

// EnsureStoreDir creates the directory a store lives in.
//
// LOCAL PATHS ONLY, and a remote base is an error rather than a no-op: on an
// s3:// locator this used to reach MkdirAll and create a local directory
// literally named "s3:". Callers that may be handed either should go through a
// Sink, which creates whatever its destination needs when a table is created
// -- there is nothing to make in advance on an object store.
//
// Note that it is not undone if the conversion later fails, so an abandoned run
// can leave an empty directory behind. That is harmless -- ExistingTables finds
// nothing in it, so the retry is not blocked -- but it does mean the presence of
// a directory says nothing about whether a store is there.
func EnsureStoreDir(base string) error {
	if scheme := schemeOf(base); scheme != "" {
		return fmt.Errorf(
			"varstore: EnsureStoreDir is for local paths; %s is a %s locator and its sink creates what it needs",
			base, scheme)
	}
	return os.MkdirAll(trimStoreDir(base), 0o755)
}

// TableFiles are the tables a store is made of, manifest included, in the
// order a conversion writes them.
//
// The manifest counts as a table here because a conversion writes one and
// removes one: leaving it out of an overwrite check would let a re-run orphan a
// marker vouching for tables it replaced.
func TableFiles() []string {
	return []string{
		TableFile(CallsTable), TableFile(SitesTable),
		TableFile(RegionsTable), VolumeManifestFile,
	}
}

// ExistingTables lists which table files already exist at base, the manifest
// included: it is written by a conversion and removed by one, so leaving it out
// would let a re-run silently orphan a marker vouching for tables it replaced.
//
// LOCAL PATHS ONLY. Use ExistingTablesIn for a base that may be remote -- this
// one answers "none" for every locator it cannot stat, which is the wrong
// answer in the one direction that matters.
func ExistingTables(base string) []string {
	var found []string
	for _, m := range []string{CallsTable, SitesTable, RegionsTable} {
		if p := TablePath(base, m); fileExists(p) {
			found = append(found, p)
		}
	}
	if p := VolumeManifestPath(base); fileExists(p) {
		found = append(found, p)
	}
	return found
}

// ExistingTablesIn lists the tables already present in a sink, by locator.
func ExistingTablesIn(sink Sink) ([]string, error) {
	var found []string
	for _, name := range TableFiles() {
		size, ok, err := sink.Stat(name)
		_ = size
		if err != nil {
			return nil, fmt.Errorf("checking %s in %s: %w", name, sink.Describe(), err)
		}
		if ok {
			found = append(found, joinStore(sink.Describe(), name))
		}
	}
	return found, nil
}

func CheckStoreTarget(base string, force bool) error {
	if force {
		return nil
	}
	// THROUGH A SINK, so the check means the same thing wherever the store is
	// going. It used to stat local paths, which against a remote base found
	// nothing and waved through exactly the overwrite it exists to prevent --
	// silently, because the tables it looks for cannot be seen with os.Stat.
	sink, err := OpenSink(base)
	if err != nil {
		return err
	}
	return CheckStoreTargetIn(sink, force)
}

// CheckStoreTargetIn is CheckStoreTarget against a sink the caller already has.
//
// Worth using when the caller goes on to write through that same sink: opening
// the destination twice is how the check and the writer come to disagree about
// where the store is.
func CheckStoreTargetIn(sink Sink, force bool) error {
	if force {
		return nil
	}
	existing, err := ExistingTablesIn(sink)
	if err != nil {
		return err
	}
	if len(existing) > 0 {
		return fmt.Errorf("refusing to overwrite an existing store: %s; pass --force to replace it",
			strings.Join(existing, ", "))
	}
	return nil
}

// ShardFile returns the file name one shard of a split table is stored under,
// relative to the store: "calls/00007.parquet".
//
// A DIRECTORY PER MEMBER rather than a directory per shard, because the tables
// are read independently -- a locus lookup touches calls and regions and may
// never open sites -- and because it keeps the unsplit layout recognisable:
// `calls.parquet` becomes `calls/`, which reads as the same thing having grown.
//
// Zero-padded so a plain listing sorts into coordinate order, which is how
// anybody debugging one will first look at it.
func ShardFile(table string, i int) string {
	return fmt.Sprintf("%s/%05d.parquet", table, i)
}
