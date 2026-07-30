package varstore

import (
	"errors"
	"fmt"
	"iter"
	"strconv"
	"strings"

	"github.com/compgenlab/cghts/vcf"
)

// State is how one sample stands at one locus.
//
// The four states exist because "no ALT call" is ambiguous on its own. A sample
// with no ALT may have been confidently called reference, or may never have been
// assayed there at all; collapsing those two into "not a carrier" silently
// inflates the denominator of every cohort query.
type State string

const (
	StateCarrier    State = "carrier"     // an ALT call passing the gate
	StateUncertain  State = "uncertain"   // an ALT call below the gate
	StateNonCarrier State = "non_carrier" // observed, and observed to be reference
	StateNotAssayed State = "not_assayed" // no observation to draw on
)

// ErrNotClassifiable is returned by Classify when the backing store cannot
// separate StateNonCarrier from StateNotAssayed. Callers must surface it rather
// than degrading to a carrier/non-carrier answer: reporting every unobserved
// sample as a non-carrier is the exact error the four states exist to prevent.
var ErrNotClassifiable = errors.New("store cannot distinguish non-carrier from not-assayed")

// Locus is a single biallelic variant.
type Locus struct {
	Chrom string
	Pos   int32
	Ref   string
	Alt   string
}

func (l Locus) String() string {
	return fmt.Sprintf("%s:%d:%s:%s", l.Chrom, l.Pos, l.Ref, l.Alt)
}

// RecordKey identifies the source record a locus was split out of. One VCF
// record holds one genotype per sample, so every alternate allele split out of
// it shares a key -- which is how a sample carrying a *different* alternate of
// the same record is told apart from one that is genuinely all-reference.
// Chromosome naming is canonicalized so spellings compare equal.
type RecordKey struct {
	Chrom string
	Pos   int32
	Ref   string
}

// Record returns the key of the record this locus came from.
func (l Locus) Record() RecordKey {
	return RecordKey{Chrom: CanonKey(l.Chrom), Pos: l.Pos, Ref: l.Ref}
}

// Span is a genomic interval in **0-based half-open** coordinates, matching
// htsio.ParseRegion and the tabix query API it feeds. VCF record positions are
// 1-based, so use Contains rather than comparing against Start/End directly.
type Span struct {
	Chrom string
	Start int32
	End   int32
}

// Contains reports whether a 1-based position falls in the span.
func (s Span) Contains(chrom string, pos int32) bool {
	if !SameChrom(s.Chrom, chrom) {
		return false
	}
	p := pos - 1 // to 0-based
	return p >= s.Start && p < s.End
}

// Gate is a per-genotype quality threshold. A zero field imposes no constraint.
//
// A gate is only as good as the fields present in the data: a store built from
// a VCF carrying no GQ has Missing in every gq column, and MinGQ then admits
// every call rather than rejecting it. Fail-open is deliberate -- absent quality
// is not evidence of poor quality -- but it means a gate can silently do nothing.
type Gate struct {
	MinDP int32

	// MinGQ gates ALT calls only, and is NOT exposed by vcf-varquery for that
	// reason. Conversion builds its callable runs from depth alone, so a Parquet
	// store retains no GQ for a genotype it never recorded -- HomRefs cannot honor
	// this field, while a VCF-backed HomRefs does, and the two backends then
	// disagree about which samples are reference. Set it only where every row
	// carries a recorded GQ.
	MinGQ int32
}

// Admits reports whether a call passes the gate. Missing values pass.
func (g Gate) Admits(c Call) bool {
	if g.MinDP > 0 && c.DP != Missing && c.DP < g.MinDP {
		return false
	}
	if g.MinGQ > 0 && c.GQ != Missing && c.GQ < g.MinGQ {
		return false
	}
	return true
}

// IsZero reports whether the gate constrains nothing.
func (g Gate) IsZero() bool { return g.MinDP <= 0 && g.MinGQ <= 0 }

// Query selects sites on one axis and samples on the other.
//
// The two axes are independent, and an empty selector imposes no restriction on
// its axis -- the way a nil Span has always read. So the zero Query asks for every
// genotype in the store, naming loci narrows the sites, naming samples narrows the
// columns, and doing both is the variants-by-samples question that used to need
// two separate calls.
//
// It is a struct rather than a set of arguments so that it can grow -- richer site
// selectors, requested INFO fields -- without changing the interface.
type Query struct {
	// Site axis, unioned. All empty means every site the store holds.
	Loci  []Locus
	Spans []Span

	// Sample axis. Empty means every sample in the store.
	//
	// Note this defaults to *all*, matching the site axis. A caller that must not
	// accidentally ask for a whole cohort should require its own selection first;
	// the library will not guess that an unset filter meant "none".
	Samples []string

	// Gate drops ALT calls below a depth threshold. The zero Gate admits everything.
	Gate Gate

	// IncludeRef also emits the reference (0/0) calls, turning "which variants does
	// this subject carry" into "every site interrogated for this subject". It needs
	// the sites and regions members, so a store built with --no-callable refuses
	// with ErrNotClassifiable rather than reporting unobserved samples as reference.
	IncludeRef bool
}

// plan indexes a Query's selectors for repeated matching during a scan.
type plan struct {
	q       Query
	loci    map[Locus]bool     // canonicalized; nil when no locus was named
	records map[RecordKey]bool // the records those loci were split out of
	samples map[string]bool    // nil when every sample is wanted
	anySite bool               // neither loci nor spans were named
}

// plan prepares the query for matching. Loci are canonicalized once here so a
// panel spelled "22" matches a store spelled "chr22" without re-canonicalizing
// per row.
func (q Query) plan() *plan {
	p := &plan{q: q, anySite: len(q.Loci) == 0 && len(q.Spans) == 0}
	if len(q.Loci) > 0 {
		p.loci = make(map[Locus]bool, len(q.Loci))
		p.records = make(map[RecordKey]bool, len(q.Loci))
		for _, l := range q.Loci {
			p.loci[canonLocus(l)] = true
			p.records[l.Record()] = true
		}
	}
	if len(q.Samples) > 0 {
		p.samples = make(map[string]bool, len(q.Samples))
		for _, s := range q.Samples {
			p.samples[s] = true
		}
	}
	return p
}

// canonLocus normalizes a locus's chromosome so differing conventions compare
// and hash equal.
func canonLocus(l Locus) Locus {
	l.Chrom = CanonKey(l.Chrom)
	return l
}

// wantsSample reports whether the sample axis admits a sample.
func (p *plan) wantsSample(name string) bool {
	return p.samples == nil || p.samples[name]
}

// wantsSite reports whether the site axis admits a locus.
func (p *plan) wantsSite(l Locus) bool {
	if p.anySite {
		return true
	}
	if p.loci != nil && p.loci[canonLocus(l)] {
		return true
	}
	for i := range p.q.Spans {
		if p.q.Spans[i].Contains(l.Chrom, l.Pos) {
			return true
		}
	}
	return false
}

// wantsRecord reports whether a locus belongs to a source record the query
// touches, which is broader than wantsSite: naming chr1:200:C:G also reaches
// chr1:200:C:T, because both were split out of one record and so share one
// genotype per sample.
//
// Reference-call reconstruction needs this. A sample carrying T is not reference
// at the G locus, and the only evidence of that is its call at the *sibling*
// locus -- which a site-level test discards, silently turning a 0/2 sample into a
// fabricated 0/0.
func (p *plan) wantsRecord(l Locus) bool {
	if p.anySite {
		return true
	}
	if p.records != nil && p.records[l.Record()] {
		return true
	}
	for i := range p.q.Spans {
		if p.q.Spans[i].Contains(l.Chrom, l.Pos) {
			return true
		}
	}
	return false
}

// wantsCall applies both axes and the gate to one ALT call.
func (p *plan) wantsCall(c Call) bool {
	return p.wantsSample(c.SampleID) && p.wantsSite(c.Locus()) && p.q.Gate.Admits(c)
}

// CollectCalls buffers a query into a slice, for callers that do not need to
// stream. It is a convenience over Calls, never a different answer.
func CollectCalls(s Store, q Query) ([]Call, error) {
	seq, err := s.Calls(q)
	if err != nil {
		return nil, err
	}
	var out []Call
	for c, err := range seq {
		if err != nil {
			return nil, err
		}
		out = append(out, c)
	}
	return out, nil
}

// SampleState pairs a sample with its state at a locus. Call is non-nil only
// for StateCarrier and StateUncertain.
type SampleState struct {
	SampleID string
	State    State
	Call     *Call
}

// Store answers carrier questions over a genotype backend.
//
// Implementations differ in how they derive a state, never in what the states
// mean -- that equivalence is the point of the interface, and the same query
// against a VCF and against a Parquet store converted from it must agree.
type Store interface {
	// Samples returns every sample the store knows about, in file order.
	Samples() ([]string, error)

	// Calls streams the genotypes a query selects.
	//
	// Rows arrive in the store's own order -- (chrom, pos) as written, then ALT
	// order, then sample roster order. That is contig order rather than
	// lexicographic, which is what an indexable VCF requires and what a
	// lexicographic sort gets wrong the moment a store holds chr10.
	//
	// The returned error reports setup failure, so a caller learns about an
	// unusable query before iterating; per-row failures arrive through the
	// iterator. With Query.IncludeRef a store lacking sites or regions fails here
	// with ErrNotClassifiable rather than silently omitting reference calls.
	//
	// Implementations choose their own access strategy from the query's shape --
	// pruned per-locus lookups for a handful of loci, one ordered pass for many --
	// so a caller cannot accidentally turn a panel into one lookup per variant.
	//
	// A reference call is emitted only where the genotype was ALL-reference. At a
	// multiallelic record a 0/2 sample is not a carrier of allele 1 but is not
	// reference either, so it appears as neither: writing 0/0 for it would
	// fabricate a genotype the source never contained.
	Calls(q Query) (iter.Seq2[Call, error], error)

	// Classify resolves every sample to a state at the locus, or returns
	// ErrNotClassifiable if the backend lacks the evidence to do so.
	Classify(l Locus, g Gate) ([]SampleState, error)

	// SiteKnown reports whether the source actually reported this locus. For a
	// plain VCF, and for a Parquet store derived from one, this is the limit of
	// what is answerable: nothing is claimed about positions the source did not
	// report, so an unknown locus yields StateNotAssayed throughout rather than
	// a set of reference calls.
	SiteKnown(l Locus) (bool, error)

	// Close releases any open files.
	Close() error
}

// ParseLocus parses "chrom:pos:ref:alt".
func ParseLocus(s string) (Locus, error) {
	f := strings.Split(s, ":")
	if len(f) != 4 {
		return Locus{}, fmt.Errorf("invalid variant %q (want chrom:pos:ref:alt)", s)
	}
	pos, err := strconv.Atoi(f[1])
	if err != nil || pos < 1 {
		return Locus{}, fmt.Errorf("invalid position in variant %q", s)
	}
	if f[2] == "" || f[3] == "" {
		return Locus{}, fmt.Errorf("invalid variant %q (ref and alt are required)", s)
	}
	return Locus{Chrom: strings.TrimSpace(f[0]), Pos: int32(pos), Ref: f[2], Alt: f[3]}, nil
}

// CanonKey reduces a chromosome name to a canonical identity so that UCSC
// ("chr17"), Ensembl ("17") and NCBI RefSeq ("NC_000017.11") spellings of the
// same chromosome compare equal. Stores record whatever the source called it;
// a query should not have to know which convention that was.
//
// Names with no canonical identity -- unplaced scaffolds, alt loci, non-human
// contigs -- fall back to exact matching, which is what cghts recommends.
func CanonKey(c string) string {
	c = strings.TrimSpace(c)
	if key, ok := vcf.CanonicalContig(c); ok {
		return key
	}
	return c
}

// SameChrom reports whether two chromosome names denote the same contig.
func SameChrom(a, b string) bool { return CanonKey(a) == CanonKey(b) }

// SameLocus compares two loci, tolerating differing chromosome conventions.
func SameLocus(a, b Locus) bool {
	return SameChrom(a.Chrom, b.Chrom) &&
		a.Pos == b.Pos && a.Ref == b.Ref && a.Alt == b.Alt
}
