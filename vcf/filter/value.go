package filter

import (
	"fmt"
	"math"
	"strconv"
	"strings"

	"github.com/compgenlab/hts/vcf"
)

// valueFilter ports ngsutilsj's VCFAbstractFilter value/flag filters: it
// addresses a single INFO or FORMAT key (optionally a specific sample, or every
// sample) and a per-allele selector, then runs a test on the resolved value.
//
// sampleName == ""     -> test every sample's value, any match flags the record
// sampleName == "INFO" -> test the INFO field
// otherwise            -> test that named sample (resolved in SetupHeader)
type valueFilter struct {
	id         string
	desc       string
	key        string
	sampleName string
	sampleIdx  int
	test       func(v vcf.AttrValue, present bool) (bool, error)
}

func newValueFilter(id, desc, key, sampleName string, test func(v vcf.AttrValue, present bool) (bool, error)) *valueFilter {
	return &valueFilter{id: id, desc: desc, key: key, sampleName: sampleName, sampleIdx: -2, test: test}
}

func (f *valueFilter) ID() string { return f.id }

func (f *valueFilter) SetupHeader(h *vcf.VcfHeader) error {
	h.AddFilter(&vcf.FilterDef{ID: f.id, Description: f.desc})
	if f.sampleName != "" && f.sampleName != "INFO" {
		f.sampleIdx = h.SampleIndex(f.sampleName)
		if f.sampleIdx < 0 {
			return fmt.Errorf("filter: unable to find sample: %s", f.sampleName)
		}
	}
	return nil
}

func (f *valueFilter) Filter(rec *vcf.VcfRecord) error {
	fail, err := f.evaluate(rec)
	if err != nil {
		return err
	}
	if fail {
		rec.AddFilter(f.id)
	}
	return nil
}

func (f *valueFilter) evaluate(rec *vcf.VcfRecord) (bool, error) {
	if f.sampleName == "INFO" {
		v, ok := rec.Info().Get(f.key)
		return f.test(v, ok)
	}
	if f.sampleIdx < 0 {
		// not set, or set as "" -> any sample matches.
		for i := 0; i < rec.NumSamples(); i++ {
			s, err := rec.Sample(i)
			if err != nil {
				return false, err
			}
			v, ok := s.Get(f.key)
			hit, err := f.test(v, ok)
			if err != nil {
				return false, err
			}
			if hit {
				return true, nil
			}
		}
		return false, nil
	}
	s, err := rec.Sample(f.sampleIdx)
	if err != nil {
		return false, err
	}
	v, ok := s.Get(f.key)
	return f.test(v, ok)
}

func (f *valueFilter) Close() error { return nil }

// --- flag filters (INFO-only) ---

// NewFlagPresent flags variants whose INFO contains the given flag key.
func NewFlagPresent(key string) Filter {
	return newFilter(key+"_present", "Contains info flag "+key, func(rec *vcf.VcfRecord) (bool, error) {
		return rec.Info().Contains(key), nil
	})
}

// NewFlagAbsent flags variants whose INFO is missing the given flag key.
func NewFlagAbsent(key string) Filter {
	return newFilter(key+"_absent", "Missing flag "+key, func(rec *vcf.VcfRecord) (bool, error) {
		return !rec.Info().Contains(key), nil
	})
}

// NewValueMissing flags variants whose key is absent or missing (".") in the
// addressed sample(s). For GT it also treats "./." and ".|." as missing. An
// empty sampleName tests every sample; "INFO" tests the INFO field.
func NewValueMissing(key, sampleName string) Filter {
	id := key + "_missing"
	loc := "all samples"
	if sampleName != "" {
		id += "_" + sampleName
		loc = "sample: " + sampleName
	}
	desc := "Missing value: " + key + " (in " + loc + ")"
	return newValueFilter(id, desc, key, sampleName, func(v vcf.AttrValue, present bool) (bool, error) {
		if !present || v.IsMissing() {
			return true, nil
		}
		if key == "GT" {
			s := v.String()
			if s == "./." || s == ".|." {
				return true, nil
			}
		}
		return false, nil
	})
}

// --- value comparison filters ---

// NewEquals flags variants where key equals val in the addressed value.
func NewEquals(key, val, sampleName, alleleName string) Filter {
	id := sanitizeID(key + "_eq_" + val + idSuffix(sampleName, alleleName))
	desc := valueDesc(key+" equals "+val, sampleName, alleleName)
	return newValueFilter(id, desc, key, sampleName, func(v vcf.AttrValue, present bool) (bool, error) {
		if alleleName == "sum" && present {
			return floatEq(v, val)
		}
		if !present {
			return val == "", nil
		}
		s, err := v.StringFor(alleleName)
		if err != nil {
			return false, err
		}
		return val == s, nil
	})
}

// NewNotEquals flags variants where key does not equal val.
func NewNotEquals(key, val, sampleName, alleleName string) Filter {
	id := sanitizeID(key + "_neq_" + val + idSuffix(sampleName, alleleName))
	desc := valueDesc(key+" not equals "+val, sampleName, alleleName)
	return newValueFilter(id, desc, key, sampleName, func(v vcf.AttrValue, present bool) (bool, error) {
		if alleleName == "sum" && present {
			eq, err := floatEq(v, val)
			return !eq, err
		}
		if !present {
			return val != "", nil
		}
		s, err := v.StringFor(alleleName)
		if err != nil {
			return false, err
		}
		return val != s, nil
	})
}

// NewContains flags variants where the addressed value contains val.
func NewContains(key, val, sampleName, alleleName string) Filter {
	id := sanitizeID(key + "_contains_" + val + idSuffix(sampleName, alleleName))
	desc := valueDesc(key+" contains "+val, sampleName, alleleName)
	return newValueFilter(id, desc, key, sampleName, func(v vcf.AttrValue, present bool) (bool, error) {
		if !present {
			return false, nil
		}
		s, err := v.StringFor(alleleName)
		if err != nil {
			return false, err
		}
		return strings.Contains(s, val), nil
	})
}

// NewNotContains flags variants where the addressed value does not contain val.
func NewNotContains(key, val, sampleName, alleleName string) Filter {
	id := sanitizeID(key + "_notcontains_" + val + idSuffix(sampleName, alleleName))
	desc := valueDesc(key+" does not contain "+val, sampleName, alleleName)
	return newValueFilter(id, desc, key, sampleName, func(v vcf.AttrValue, present bool) (bool, error) {
		if !present {
			return false, nil
		}
		s, err := v.StringFor(alleleName)
		if err != nil {
			return false, err
		}
		return !strings.Contains(s, val), nil
	})
}

// NewInList flags variants where the addressed value equals any of vals.
func NewInList(key string, vals []string, sampleName, alleleName string) Filter {
	id := sanitizeID(key + "_in_" + strings.Join(vals, "_") + idSuffix(sampleName, alleleName))
	desc := valueDesc(key+" in "+strings.Join(vals, ", "), sampleName, alleleName)
	return newValueFilter(id, desc, key, sampleName, func(v vcf.AttrValue, present bool) (bool, error) {
		if !present {
			return false, nil
		}
		q, err := v.StringFor(alleleName)
		if err != nil {
			return false, err
		}
		for _, x := range vals {
			if x == q {
				return true, nil
			}
		}
		return false, nil
	})
}

// NewNotInList flags variants where the addressed value is not any of vals (an
// absent value counts as not-in-list, matching ngsutilsj).
func NewNotInList(key string, vals []string, sampleName, alleleName string) Filter {
	id := sanitizeID(key + "_notin_" + strings.Join(vals, "_") + idSuffix(sampleName, alleleName))
	desc := valueDesc(key+" not in "+strings.Join(vals, ", "), sampleName, alleleName)
	return newValueFilter(id, desc, key, sampleName, func(v vcf.AttrValue, present bool) (bool, error) {
		if present {
			q, err := v.StringFor(alleleName)
			if err != nil {
				return false, err
			}
			for _, x := range vals {
				if x == q {
					return false, nil
				}
			}
		}
		return true, nil
	})
}

// --- numeric comparison filters ---

// NewLessThan flags variants where the addressed numeric value is < thres.
func NewLessThan(key string, thres float64, sampleName, alleleName string) Filter {
	return newMathFilter("lt", "less than", key, thres, sampleName, alleleName, func(d float64) bool { return d < thres })
}

// NewLessThanEqual flags variants where the addressed numeric value is <= thres.
func NewLessThanEqual(key string, thres float64, sampleName, alleleName string) Filter {
	return newMathFilter("lte", "less than or equal", key, thres, sampleName, alleleName, func(d float64) bool { return d <= thres })
}

// NewGreaterThan flags variants where the addressed numeric value is > thres.
func NewGreaterThan(key string, thres float64, sampleName, alleleName string) Filter {
	return newMathFilter("gt", "greater than", key, thres, sampleName, alleleName, func(d float64) bool { return d > thres })
}

// NewGreaterThanEqual flags variants where the addressed numeric value is >= thres.
func NewGreaterThanEqual(key string, thres float64, sampleName, alleleName string) Filter {
	return newMathFilter("gte", "greater than or equal", key, thres, sampleName, alleleName, func(d float64) bool { return d >= thres })
}

func newMathFilter(idTag, descTag, key string, thres float64, sampleName, alleleName string, op func(d float64) bool) Filter {
	t := javaDouble(thres)
	// Math-filter IDs are NOT sanitized in ngsutilsj.
	id := key + "_" + idTag + "_" + t + idSuffix(sampleName, alleleName)
	desc := valueDesc(key+" "+descTag+" "+t, sampleName, alleleName)
	return newValueFilter(id, desc, key, sampleName, func(v vcf.AttrValue, present bool) (bool, error) {
		if !present || v.IsMissing() {
			return false, nil
		}
		d, err := v.FloatFor(alleleName)
		if err != nil {
			return false, err
		}
		if math.IsNaN(d) {
			return true, nil
		}
		return op(d), nil
	})
}

// --- shared helpers ---

func floatEq(v vcf.AttrValue, val string) (bool, error) {
	d, err := v.FloatFor("sum")
	if err != nil {
		return false, err
	}
	t, err := strconv.ParseFloat(val, 64)
	if err != nil {
		return false, err
	}
	return d == t, nil
}

// valueDesc builds a "PREFIX (in all samples[, allele: X])" description matching
// ngsutilsj's generated FILTER descriptions.
func valueDesc(prefix, sampleName, alleleName string) string {
	s := prefix + " (in "
	if sampleName == "" {
		s += "all samples"
	} else {
		s += "sample: " + sampleName
	}
	if alleleName == "" {
		s += ")"
	} else {
		s += ", allele: " + alleleName + ")"
	}
	return s
}

// idSuffix appends "_SAMPLE" and "_ALLELE" when present, as ngsutilsj does.
func idSuffix(sampleName, alleleName string) string {
	s := ""
	if sampleName != "" {
		s += "_" + sampleName
	}
	if alleleName != "" {
		s += "_" + alleleName
	}
	return s
}

// sanitizeID replaces characters ngsutilsj strips from generated FILTER IDs
// (",;<> \t\r\n") with "_" and collapses runs of "_".
func sanitizeID(s string) string {
	repl := strings.NewReplacer(
		",", "_", ";", "_", ">", "_", "<", "_",
		" ", "_", "\t", "_", "\r", "_", "\n", "_",
	)
	s = repl.Replace(s)
	for strings.Contains(s, "__") {
		s = strings.ReplaceAll(s, "__", "_")
	}
	return s
}
