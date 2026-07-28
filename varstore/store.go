package varstore

import (
	"errors"
	"fmt"
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

	// Carriers returns the ALT-carrying calls at a locus that pass the gate.
	Carriers(l Locus, g Gate) ([]Call, error)

	// Variants returns the ALT-carrying calls for one sample, optionally
	// restricted to a span. A nil span means the whole store.
	Variants(sample string, span *Span, g Gate) ([]Call, error)

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
