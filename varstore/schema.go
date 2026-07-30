// Package varstore holds the on-disk schema for a sparse genotype store and a
// uniform way to query genotypes regardless of which format backs them.
//
// A Parquet store is a set of three files derived from one base name. All three
// are required: the calls file alone cannot distinguish a confidently-called
// reference genotype from a position that was never assayed, which is the
// distinction Classify exists to make.
//
//	BASE.calls.parquet     one row per ALT-carrying genotype
//	BASE.sites.parquet     one row per interrogated site, sample-independent
//	BASE.regions.parquet   runs of catalog sites at which a sample was called
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
	"path/filepath"
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
	GT       string `parquet:"gt,dict"`
	DP       int32  `parquet:"dp"`
	ADRef    int32  `parquet:"ad_ref"`
	ADAlt    int32  `parquet:"ad_alt"`
	GQ       int32  `parquet:"gq"`
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
}

// The three members of a store, and the two ways a base name can address them.
//
// A base ending in a path separator names a DIRECTORY, and the members sit
// inside it under their bare names:
//
//	--out cohort/   ->  cohort/calls.parquet, cohort/sites.parquet, ...
//
// Any other base is a filename PREFIX, and the member name is appended:
//
//	--out cohort    ->  cohort.calls.parquet, cohort.sites.parquet, ...
//
// The directory form keeps the set as one thing to copy, move or delete, which
// matters because the three files are only meaningful together.
const (
	CallsMember   = "calls"
	SitesMember   = "sites"
	RegionsMember = "regions"

	CallsSuffix   = "." + CallsMember + ".parquet"
	SitesSuffix   = "." + SitesMember + ".parquet"
	RegionsSuffix = "." + RegionsMember + ".parquet"
)

// IsDirBase reports whether a base names a directory rather than a filename
// prefix, i.e. whether it ends in a path separator.
func IsDirBase(base string) bool {
	return strings.HasSuffix(base, "/") || strings.HasSuffix(base, string(os.PathSeparator))
}

// MemberPath returns the file holding one member of the store at base.
func MemberPath(base, member string) string {
	if IsDirBase(base) {
		return filepath.Join(base, member+".parquet")
	}
	return base + "." + member + ".parquet"
}

// CallsPath returns the calls file for a store base name.
func CallsPath(base string) string { return MemberPath(base, CallsMember) }

// SitesPath returns the sites file for a store base name.
func SitesPath(base string) string { return MemberPath(base, SitesMember) }

// RegionsPath returns the callable-regions file for a store base name.
func RegionsPath(base string) string { return MemberPath(base, RegionsMember) }

// EnsureStoreDir creates the containing directory for a directory-form base.
// It is a no-op for the prefix form, whose parent is the caller's business.
func EnsureStoreDir(base string) error {
	if !IsDirBase(base) {
		return nil
	}
	return os.MkdirAll(base, 0o755)
}

// ExistingMembers lists which of the three member files already exist at base.
func ExistingMembers(base string) []string {
	var found []string
	for _, m := range []string{CallsMember, SitesMember, RegionsMember} {
		if p := MemberPath(base, m); fileExists(p) {
			found = append(found, p)
		}
	}
	return found
}

// CheckStoreTarget refuses to write over an existing store unless force is set.
//
// Conversion truncates all three members, so an accidental re-run against a
// populated base destroys the previous store outright -- and because the three
// files are only meaningful together, a half-replaced set is worse than either
// outcome. Any single surviving member is therefore enough to stop.
//
// A prefix-form base that is itself a directory is refused separately: writing
// "cohort.calls.parquet" beside an existing "cohort/" directory is legal but
// almost certainly means the trailing slash was forgotten.
func CheckStoreTarget(base string, force bool) error {
	if force {
		return nil
	}
	if !IsDirBase(base) {
		if st, err := os.Stat(base); err == nil && st.IsDir() {
			return fmt.Errorf("%s is a directory; write inside it with --out %s%c, "+
				"or pass --force to use the name as a filename prefix",
				base, base, os.PathSeparator)
		}
	}
	if existing := ExistingMembers(base); len(existing) > 0 {
		return fmt.Errorf("refusing to overwrite an existing store: %s; pass --force to replace it",
			strings.Join(existing, ", "))
	}
	return nil
}
