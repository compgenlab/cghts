// Package filter provides a composable set of VCF filters. Each [Filter]
// examines a record and, when the record fails the filter's test, stamps the
// filter's ID onto the record's FILTER column (via [vcf.VcfRecord.AddFilter]).
// A record with no stamped codes is PASS. Filters are reusable from any Go code,
// not just the cgkit CLI. It ports ngsutilsj's vcf/filter framework.
package filter

import (
	"strconv"
	"strings"

	"github.com/compgenlab/cghts/vcf"
)

// Filter tests a record and marks it (with its ID) when the record fails.
type Filter interface {
	ID() string
	SetupHeader(h *vcf.VcfHeader) error
	Filter(rec *vcf.VcfRecord) error
	Close() error
}

// filterFunc is the common implementation: SetupHeader adds the ##FILTER def and
// Filter stamps the ID when test reports the record should be filtered.
type filterFunc struct {
	id   string
	desc string
	test func(rec *vcf.VcfRecord) (bool, error)
}

func newFilter(id, desc string, test func(rec *vcf.VcfRecord) (bool, error)) Filter {
	return &filterFunc{id: id, desc: desc, test: test}
}

func (f *filterFunc) ID() string { return f.id }

func (f *filterFunc) SetupHeader(h *vcf.VcfHeader) error {
	h.AddFilter(&vcf.FilterDef{ID: f.id, Description: f.desc})
	return nil
}

func (f *filterFunc) Filter(rec *vcf.VcfRecord) error {
	fail, err := f.test(rec)
	if err != nil {
		return err
	}
	if fail {
		rec.AddFilter(f.id)
	}
	return nil
}

func (f *filterFunc) Close() error { return nil }

// Chain applies a set of filters to each record in order. SetupHeaders adds
// every filter's ##FILTER def (in order); Apply runs each filter on a record.
type Chain struct {
	filters []Filter
}

// NewChain returns an empty Chain.
func NewChain() *Chain { return &Chain{} }

// Add appends a filter.
func (c *Chain) Add(f Filter) { c.filters = append(c.filters, f) }

// Len reports how many filters are in the chain.
func (c *Chain) Len() int { return len(c.filters) }

// SetupHeaders runs each filter's header step in order.
func (c *Chain) SetupHeaders(h *vcf.VcfHeader) error {
	for _, f := range c.filters {
		if err := f.SetupHeader(h); err != nil {
			return err
		}
	}
	return nil
}

// Apply runs each filter on the record (stamping FILTER codes as they fail).
func (c *Chain) Apply(rec *vcf.VcfRecord) error {
	for _, f := range c.filters {
		if err := f.Filter(rec); err != nil {
			return err
		}
	}
	return nil
}

// Close releases every filter.
func (c *Chain) Close() error {
	var firstErr error
	for _, f := range c.filters {
		if err := f.Close(); err != nil && firstErr == nil {
			firstErr = err
		}
	}
	return firstErr
}

// --- self-contained filters ---

// NewChromPass flags variants NOT on any of the given chromosomes.
func NewChromPass(chroms []string) Filter {
	set := toSet(chroms)
	id := "CHROM_PASS_" + strings.Join(chroms, "_")
	desc := "Chromosome is not " + strings.Join(chroms, ",")
	return newFilter(id, desc, func(rec *vcf.VcfRecord) (bool, error) {
		return !set[rec.Chrom], nil
	})
}

// NewChromFail flags variants on any of the given chromosomes.
func NewChromFail(chroms []string) Filter {
	set := toSet(chroms)
	id := "CHROM_FAIL_" + strings.Join(chroms, "_")
	desc := "Chromosome is " + strings.Join(chroms, ",")
	return newFilter(id, desc, func(rec *vcf.VcfRecord) (bool, error) {
		return set[rec.Chrom], nil
	})
}

// NewSNV flags variants that are SNVs (REF and every ALT length 1).
func NewSNV() Filter {
	return newFilter("SNV", "Variant is an SNV", func(rec *vcf.VcfRecord) (bool, error) {
		if len(rec.Ref) != 1 {
			return false, nil
		}
		for _, a := range rec.Alt() {
			if len(a) != 1 {
				return false, nil
			}
		}
		return true, nil
	})
}

// NewIndel flags variants that are indels (REF or any ALT length != 1).
func NewIndel() Filter {
	return newFilter("INDEL", "Variant is an indel", func(rec *vcf.VcfRecord) (bool, error) {
		if len(rec.Ref) != 1 {
			return true, nil
		}
		for _, a := range rec.Alt() {
			if len(a) != 1 {
				return true, nil
			}
		}
		return false, nil
	})
}

// NewMaxIns flags variants with an insertion longer than thres.
func NewMaxIns(thres int) Filter {
	return newFilter("INS_max_"+strconv.Itoa(thres), "Insertion longer than "+strconv.Itoa(thres),
		func(rec *vcf.VcfRecord) (bool, error) {
			for _, a := range rec.Alt() {
				if len(a)-1 > thres {
					return true, nil
				}
			}
			return false, nil
		})
}

// NewMaxDel flags variants with a deletion longer than thres.
func NewMaxDel(thres int) Filter {
	return newFilter("DEL_max_"+strconv.Itoa(thres), "Deletion longer than "+strconv.Itoa(thres),
		func(rec *vcf.VcfRecord) (bool, error) {
			return len(rec.Ref)-1 > thres, nil
		})
}

// NewQual flags variants whose QUAL is less than thres.
func NewQual(thres float64) Filter {
	t := javaDouble(thres)
	return newFilter("QUAL_lt_"+t, "Quality score less than "+t, func(rec *vcf.VcfRecord) (bool, error) {
		return rec.Qual() < thres, nil
	})
}

// NewHomozygous flags variants that are homozygous in any sample (GT a/a).
func NewHomozygous() Filter {
	return newFilter("homozygous", "Homozygous variant (in any sample)", gtTest(func(a, b string) bool { return a == b }))
}

// NewHeterozygous flags variants that are heterozygous in any sample (GT a/b).
func NewHeterozygous() Filter {
	// ngsutilsj's description carries the original "Heterzygous" spelling.
	return newFilter("heterozygous", "Heterzygous variant (in any sample)", gtTest(func(a, b string) bool { return a != b }))
}

// gtTest builds a sample-GT test: it reports true when any sample's diploid GT
// satisfies match(allele0, allele1).
func gtTest(match func(a, b string) bool) func(rec *vcf.VcfRecord) (bool, error) {
	return func(rec *vcf.VcfRecord) (bool, error) {
		for i := 0; i < rec.NumSamples(); i++ {
			s, err := rec.Sample(i)
			if err != nil {
				return false, err
			}
			gt, ok := s.Get("GT")
			if !ok {
				continue
			}
			v := gt.String()
			var parts []string
			if strings.Contains(v, "/") {
				parts = strings.Split(v, "/")
			} else if strings.Contains(v, "|") {
				parts = strings.Split(v, "|")
			} else {
				continue
			}
			if len(parts) == 2 && match(parts[0], parts[1]) {
				return true, nil
			}
		}
		return false, nil
	}
}

func toSet(vals []string) map[string]bool {
	m := make(map[string]bool, len(vals))
	for _, v := range vals {
		m[v] = true
	}
	return m
}

// javaDouble formats a float the way Java's Double.toString does for ordinary
// magnitudes (an integral value keeps ".0"), to match ngsutilsj's generated IDs.
func javaDouble(f float64) string {
	s := strconv.FormatFloat(f, 'f', -1, 64)
	if !strings.Contains(s, ".") {
		s += ".0"
	}
	return s
}
