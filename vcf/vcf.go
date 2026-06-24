package vcf

import (
	"fmt"
	"math"
	"strconv"
	"strings"
)

const (
	// pass is the FILTER value for a variant that passed all filters.
	pass = "PASS"
	// missing is the VCF missing-value marker.
	missing = "."
)

// AttrValue is a single INFO or FORMAT value. It carries the raw string and
// defers all type interpretation to the accessor methods, mirroring ngsutilsj's
// VCFAttributeValue. A missing value is "." and a bare flag is "".
type AttrValue struct {
	raw string
}

// String returns the raw value.
func (v AttrValue) String() string { return v.raw }

// IsMissing reports whether the value is the missing marker ".".
func (v AttrValue) IsMissing() bool { return v.raw == missing }

// IsEmpty reports whether the value is empty (a bare INFO flag).
func (v AttrValue) IsEmpty() bool { return v.raw == "" }

// Int parses the value as an integer.
func (v AttrValue) Int() (int, error) { return strconv.Atoi(v.raw) }

// Float parses the value as a float. An empty or missing value yields NaN,
// matching VCFAttributeValue.asDouble(null).
func (v AttrValue) Float() (float64, error) {
	if v.raw == "" || v.raw == missing {
		return math.NaN(), nil
	}
	return strconv.ParseFloat(v.raw, 64)
}

// StringFor extracts a string for a multi-allele selector. The selector is one
// of: "" (whole value), "ref" (first), "alt1" (second), an integer index, or an
// aggregate ("sum", "nref", "min", "max"). It ports VCFAttributeValue.asString.
func (v AttrValue) StringFor(sel string) (string, error) {
	if sel == "" {
		return v.raw, nil
	}
	parts := strings.Split(v.raw, ",")
	switch sel {
	case "ref":
		return parts[0], nil
	case "alt1":
		if len(parts) < 2 {
			return "", fmt.Errorf("vcf: no alt1 allele in %q", v.raw)
		}
		return parts[1], nil
	default:
		if i, err := strconv.Atoi(sel); err == nil {
			if i < 0 || i >= len(parts) {
				return "", fmt.Errorf("vcf: allele index %d out of range in %q", i, v.raw)
			}
			return parts[i], nil
		}
		d, err := v.FloatFor(sel)
		if err != nil {
			return "", err
		}
		if math.IsNaN(d) {
			return "", fmt.Errorf("vcf: unable to find allele: %s", sel)
		}
		return formatFloat(d), nil
	}
}

// FloatFor extracts a float for a multi-allele selector (see [AttrValue.StringFor]
// for the selector forms). It ports VCFAttributeValue.asDouble.
func (v AttrValue) FloatFor(sel string) (float64, error) {
	switch sel {
	case "":
		return v.Float()
	case "sum":
		return v.aggregate(0, sumAgg)
	case "nref":
		return v.aggregate(1, sumAgg)
	case "min":
		return v.aggregate(0, minAgg)
	case "max":
		return v.aggregate(0, maxAgg)
	default:
		s, err := v.StringFor(sel)
		if err != nil {
			return 0, err
		}
		if s == "" || s == missing {
			return math.NaN(), nil
		}
		f, err := strconv.ParseFloat(s, 64)
		if err != nil {
			return 0, fmt.Errorf("vcf: invalid value %q, expected a number", s)
		}
		return f, nil
	}
}

type aggKind int

const (
	sumAgg aggKind = iota
	minAgg
	maxAgg
)

func (v AttrValue) aggregate(from int, kind aggKind) (float64, error) {
	parts := strings.Split(v.raw, ",")
	acc := 0.0
	ext := math.NaN()
	for i := from; i < len(parts); i++ {
		if parts[i] == "" {
			continue
		}
		d, err := strconv.ParseFloat(parts[i], 64)
		if err != nil {
			return 0, fmt.Errorf("vcf: invalid value %q, expected a number", parts[i])
		}
		switch kind {
		case sumAgg:
			acc += d
		case minAgg:
			if math.IsNaN(ext) || d < ext {
				ext = d
			}
		case maxAgg:
			if math.IsNaN(ext) || d > ext {
				ext = d
			}
		}
	}
	if kind == sumAgg {
		return acc, nil
	}
	return ext, nil
}

// Attributes is an ordered collection of INFO or FORMAT key/value pairs.
type Attributes struct {
	keys []string
	vals map[string]AttrValue
}

func newAttributes() *Attributes {
	return &Attributes{vals: map[string]AttrValue{}}
}

func (a *Attributes) put(key string, v AttrValue) {
	if _, ok := a.vals[key]; !ok {
		a.keys = append(a.keys, key)
	}
	a.vals[key] = v
}

// Get returns the value for key. The boolean is false when the key is absent
// (distinct from a present-but-missing "." value).
func (a *Attributes) Get(key string) (AttrValue, bool) {
	v, ok := a.vals[key]
	return v, ok
}

// Contains reports whether key is present.
func (a *Attributes) Contains(key string) bool {
	_, ok := a.vals[key]
	return ok
}

// Keys returns the keys in insertion order.
func (a *Attributes) Keys() []string { return a.keys }

// FindKeys returns the keys matching the given glob (supporting * and ?).
func (a *Attributes) FindKeys(glob string) []string {
	var out []string
	for _, k := range a.keys {
		if globMatch(k, glob) {
			out = append(out, k)
		}
	}
	return out
}

// VcfRecord is a single VCF data line. The leading CHROM/POS/REF columns are
// parsed eagerly; everything else is parsed on first access and cached. See the
// package documentation for the lazy-parsing contract.
type VcfRecord struct {
	line   string
	header *VcfHeader

	tabs  [8]int // byte offsets of the first up-to-8 tab characters
	ntabs int

	Chrom string
	Pos   int // 1-based
	Ref   string

	altDone bool
	alt     []string

	qualDone bool
	qual     float64

	filtDone bool
	filters  []string

	infoDone bool
	info     *Attributes

	fmtDone   bool
	fmtKeys   []string
	formatRaw string
	sampleRaw []string
	samples   []*Attributes
}

// fixedCol returns fixed column k (0..7) and whether it is present.
func (r *VcfRecord) fixedCol(k int) (string, bool) {
	if k < 0 || k > 7 {
		return "", false
	}
	var start int
	if k == 0 {
		start = 0
	} else {
		if k-1 >= r.ntabs {
			return "", false
		}
		start = r.tabs[k-1] + 1
	}
	var end int
	switch {
	case k < r.ntabs:
		end = r.tabs[k]
	case k == r.ntabs && r.ntabs < 8:
		end = len(r.line)
	default:
		return "", false
	}
	return r.line[start:end], true
}

// afterInfo returns the raw FORMAT-and-samples portion of the line (everything
// after the INFO column), or "" when the record has no sample data.
func (r *VcfRecord) afterInfo() string {
	if r.ntabs < 8 {
		return ""
	}
	return r.line[r.tabs[7]+1:]
}

// prefixThroughInfo returns the raw line up to and including the INFO column.
func (r *VcfRecord) prefixThroughInfo() string {
	if r.ntabs < 8 {
		return r.line
	}
	return r.line[:r.tabs[7]]
}

func newRecord(line string, header *VcfHeader) (*VcfRecord, error) {
	line = strings.TrimRight(line, "\r\n")

	r := &VcfRecord{line: line, header: header}
	for i := 0; i < len(line) && r.ntabs < 8; i++ {
		if line[i] == '\t' {
			r.tabs[r.ntabs] = i
			r.ntabs++
		}
	}

	// Need at least CHROM, POS, ID, REF, ALT (4 tabs).
	if r.ntabs < 4 {
		return nil, fmt.Errorf("vcf: too few columns: %q", line)
	}

	chrom, _ := r.fixedCol(0)
	posStr, _ := r.fixedCol(1)
	ref, _ := r.fixedCol(3)
	pos, err := strconv.Atoi(posStr)
	if err != nil {
		return nil, fmt.Errorf("vcf: invalid POS %q: %w", posStr, err)
	}
	r.Chrom = chrom
	r.Pos = pos
	r.Ref = ref
	return r, nil
}

// Line returns the raw source line (without a trailing newline).
func (r *VcfRecord) Line() string { return r.line }

// Header returns the header this record was parsed against (may be nil).
func (r *VcfRecord) Header() *VcfHeader { return r.header }

// ID returns the ID column, or "" when it is the missing marker ".".
func (r *VcfRecord) ID() string {
	s, _ := r.fixedCol(2)
	if s == missing {
		return ""
	}
	return s
}

// AltOrig returns the raw ALT column verbatim.
func (r *VcfRecord) AltOrig() string {
	s, _ := r.fixedCol(4)
	return s
}

// Alt returns the alternate alleles, dropping "." entries. It is nil when the
// ALT column is ".".
func (r *VcfRecord) Alt() []string {
	if !r.altDone {
		raw, _ := r.fixedCol(4)
		for _, a := range strings.Split(raw, ",") {
			if a != missing {
				r.alt = append(r.alt, a)
			}
		}
		r.altDone = true
	}
	return r.alt
}

// Qual returns the QUAL value, or -1 when it is missing.
func (r *VcfRecord) Qual() float64 {
	if !r.qualDone {
		r.qual = -1
		if raw, ok := r.fixedCol(5); ok && raw != missing {
			if f, err := strconv.ParseFloat(raw, 64); err == nil {
				r.qual = f
			}
		}
		r.qualDone = true
	}
	return r.qual
}

// Filters returns the FILTER codes. It is nil when the record passed (PASS) and
// an empty non-nil slice when the column was ".".
func (r *VcfRecord) Filters() []string {
	if !r.filtDone {
		if raw, ok := r.fixedCol(6); ok && raw != pass {
			r.filters = []string{}
			for _, f := range strings.Split(raw, ";") {
				if f != missing {
					if r.header == nil || r.header.filterAllowed(f) {
						r.filters = append(r.filters, f)
					}
				}
			}
		}
		r.filtDone = true
	}
	return r.filters
}

// IsFiltered reports whether the record carries any (non-PASS) filter.
func (r *VcfRecord) IsFiltered() bool {
	return len(r.Filters()) > 0
}

// Info returns the parsed INFO attributes, parsing the INFO column on first
// access.
func (r *VcfRecord) Info() *Attributes {
	if !r.infoDone {
		r.info = newAttributes()
		if raw, ok := r.fixedCol(7); ok && raw != missing {
			for _, el := range strings.Split(raw, ";") {
				if eq := strings.IndexByte(el, '='); eq < 0 {
					r.info.put(el, AttrValue{raw: ""})
				} else {
					r.info.put(el[:eq], AttrValue{raw: el[eq+1:]})
				}
			}
		}
		r.infoDone = true
	}
	return r.info
}

// InfoValue is a shortcut for Info().Get(key).
func (r *VcfRecord) InfoValue(key string) (AttrValue, bool) {
	return r.Info().Get(key)
}

func (r *VcfRecord) ensureFormat() {
	if r.fmtDone {
		return
	}
	r.fmtDone = true
	rest := r.afterInfo()
	if rest == "" {
		return
	}
	cols := strings.Split(rest, "\t")
	r.formatRaw = cols[0]
	r.fmtKeys = strings.Split(cols[0], ":")
	if len(cols) > 1 {
		r.sampleRaw = cols[1:]
		r.samples = make([]*Attributes, len(r.sampleRaw))
	}
}

// ReorderSamplesLine returns the record's raw line with its sample columns
// permuted by order (each entry is an original 0-based sample index). FORMAT
// values are not parsed; the raw sample columns are moved verbatim. An
// out-of-range index emits ".".
func (r *VcfRecord) ReorderSamplesLine(order []int) string {
	r.ensureFormat()
	var b strings.Builder
	b.WriteString(r.prefixThroughInfo())
	if len(r.sampleRaw) == 0 {
		return b.String()
	}
	b.WriteString("\t")
	b.WriteString(r.formatRaw)
	for _, idx := range order {
		b.WriteByte('\t')
		if idx >= 0 && idx < len(r.sampleRaw) {
			b.WriteString(r.sampleRaw[idx])
		} else {
			b.WriteString(missing)
		}
	}
	return b.String()
}

// FormatKeys returns the FORMAT keys (the colon-separated keys in column 9).
func (r *VcfRecord) FormatKeys() []string {
	r.ensureFormat()
	return r.fmtKeys
}

// NumSamples returns the number of sample columns.
func (r *VcfRecord) NumSamples() int {
	r.ensureFormat()
	return len(r.sampleRaw)
}

// Sample returns the parsed FORMAT attributes for sample i, parsing that
// sample's column on first access. Other samples are left unparsed.
func (r *VcfRecord) Sample(i int) (*Attributes, error) {
	r.ensureFormat()
	if i < 0 || i >= len(r.sampleRaw) {
		return nil, fmt.Errorf("vcf: sample index %d out of range (%d samples)", i, len(r.sampleRaw))
	}
	if r.samples[i] == nil {
		attrs := newAttributes()
		vals := strings.Split(r.sampleRaw[i], ":")
		for k, key := range r.fmtKeys {
			if k < len(vals) {
				attrs.put(key, AttrValue{raw: vals[k]})
			} else {
				attrs.put(key, AttrValue{raw: missing})
			}
		}
		r.samples[i] = attrs
	}
	return r.samples[i], nil
}

// SampleByName returns the parsed FORMAT attributes for the named sample.
func (r *VcfRecord) SampleByName(name string) (*Attributes, error) {
	if r.header == nil {
		return nil, fmt.Errorf("vcf: missing header, cannot resolve sample %q", name)
	}
	i := r.header.SampleIndex(name)
	if i < 0 {
		return nil, fmt.Errorf("vcf: sample not found: %s", name)
	}
	return r.Sample(i)
}

// ZeroBasedStart returns Pos-1, the 0-based start used for BED-style output.
func (r *VcfRecord) ZeroBasedStart() int { return r.Pos - 1 }

// IsIndel reports whether the REF or any ALT allele is longer than one base.
func (r *VcfRecord) IsIndel() bool {
	if len(r.Ref) != 1 {
		return true
	}
	for _, a := range r.Alt() {
		if len(a) != 1 {
			return true
		}
	}
	return false
}

// CalcTsTv classifies a SNV as a transition (-1) or transversion (1). It
// returns 0 for anything that is not a single-base biallelic substitution
// (indels, multiallelic sites, MNVs). It ports VCFRecord.calcTsTv.
func (r *VcfRecord) CalcTsTv() int {
	if len(r.Ref) != 1 {
		return 0
	}
	alt := r.Alt()
	if len(alt) != 1 || len(alt[0]) != 1 {
		return 0
	}
	ref := strings.ToUpper(r.Ref)
	a := strings.ToUpper(alt[0])
	if ref == a {
		return 0
	}
	switch ref {
	case "A":
		if a == "G" {
			return -1
		}
		return 1
	case "G":
		if a == "A" {
			return -1
		}
		return 1
	case "C":
		if a == "T" {
			return -1
		}
		return 1
	case "T":
		if a == "C" {
			return -1
		}
		return 1
	}
	return 0
}

// formatFloat renders a float as a plain decimal, trimming a trailing ".0".
func formatFloat(f float64) string {
	s := strconv.FormatFloat(f, 'f', -1, 64)
	return strings.TrimSuffix(s, ".0")
}

// globMatch reports whether s matches a glob pattern with * (any run) and ?
// (any single character).
func globMatch(s, pattern string) bool {
	return globHelper(s, pattern)
}

func globHelper(s, p string) bool {
	for len(p) > 0 {
		switch p[0] {
		case '*':
			// Collapse consecutive stars.
			for len(p) > 1 && p[1] == '*' {
				p = p[1:]
			}
			if len(p) == 1 {
				return true
			}
			for i := 0; i <= len(s); i++ {
				if globHelper(s[i:], p[1:]) {
					return true
				}
			}
			return false
		case '?':
			if len(s) == 0 {
				return false
			}
			s = s[1:]
			p = p[1:]
		default:
			if len(s) == 0 || s[0] != p[0] {
				return false
			}
			s = s[1:]
			p = p[1:]
		}
	}
	return len(s) == 0
}
