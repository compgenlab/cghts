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

// Set adds or replaces a key with a raw string value.
func (a *Attributes) Set(key, value string) { a.put(key, AttrValue{raw: value}) }

// SetFlag adds a key as a bare flag (empty value).
func (a *Attributes) SetFlag(key string) { a.put(key, AttrValue{raw: ""}) }

// SetValue adds or replaces a key with an [AttrValue].
func (a *Attributes) SetValue(key string, v AttrValue) { a.put(key, v) }

// Remove deletes a key.
func (a *Attributes) Remove(key string) {
	if _, ok := a.vals[key]; !ok {
		return
	}
	delete(a.vals, key)
	for i, k := range a.keys {
		if k == key {
			a.keys = append(a.keys[:i], a.keys[i+1:]...)
			break
		}
	}
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

// infoString renders the attributes as a VCF INFO field: bare key for a flag,
// "key=value" otherwise, in insertion order; "." when empty. Ports
// VCFAttributes.toString.
func (a *Attributes) infoString() string {
	if len(a.keys) == 0 {
		return missing
	}
	parts := make([]string, 0, len(a.keys))
	for _, k := range a.keys {
		v := a.vals[k]
		if v.IsEmpty() {
			parts = append(parts, k)
		} else {
			parts = append(parts, k+"="+v.raw)
		}
	}
	return strings.Join(parts, ";")
}

// formatString renders the attributes as a per-sample FORMAT value for the given
// ordered format keys: missing ("." ) for absent keys, GT kept if present, and
// trailing missing values trimmed. Ports VCFAttributes.toString(format).
func (a *Attributes) formatString(format []string) string {
	vals := make([]string, len(format))
	for i, k := range format {
		if v, ok := a.vals[k]; ok {
			vals[i] = v.raw
		} else {
			vals[i] = missing
		}
	}
	limit := 0
	if len(format) > 0 && format[0] == "GT" {
		limit = 1
	}
	end := len(vals)
	for end > limit && vals[end-1] == missing {
		end--
	}
	return strings.Join(vals[:end], ":")
}

// VcfRecord is a single VCF data line. The leading CHROM/POS/REF columns are
// parsed eagerly; everything else is parsed on first access and cached. See the
// package documentation for the lazy-parsing contract.
type VcfRecord struct {
	line   string
	header *VcfHeader

	tabs  [8]int // byte offsets of the first up-to-8 tab characters
	ntabs int

	// dirty is set by any mutation; a dirty record is reconstructed from the
	// parsed model on write rather than emitted verbatim.
	dirty bool

	Chrom string
	Pos   int // 1-based
	Ref   string

	idDone bool
	id     string

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

// NewRecord builds a VcfRecord from a bare CHROM/POS/REF/ALT tuple, with no ID,
// QUAL, FILTER (PASS), INFO, or samples. The record is "dirty" — it has no
// backing line, so it is always serialized from its model on write. Annotators
// can mutate it like any parsed record. This lets annotation code run on plain
// variant tuples, not just parsed VCF lines.
func NewRecord(chrom string, pos int, ref string, alt []string) *VcfRecord {
	return &VcfRecord{
		Chrom:    chrom,
		Pos:      pos,
		Ref:      ref,
		dirty:    true,
		idDone:   true,
		altDone:  true,
		alt:      append([]string(nil), alt...),
		qualDone: true,
		qual:     -1,
		filtDone: true, // nil filters => PASS
		infoDone: true,
		info:     newAttributes(),
		fmtDone:  true, // no samples
	}
}

// Line returns the raw source line (without a trailing newline).
func (r *VcfRecord) Line() string { return r.line }

// Header returns the header this record was parsed against (may be nil).
func (r *VcfRecord) Header() *VcfHeader { return r.header }

// ID returns the ID column, or "" when it is the missing marker ".".
func (r *VcfRecord) ID() string {
	if !r.idDone {
		s, _ := r.fixedCol(2)
		if s == missing {
			s = ""
		}
		r.id = s
		r.idDone = true
	}
	return r.id
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

// Dirty reports whether the record has been modified since it was read (and so
// must be reconstructed on write rather than emitted verbatim).
func (r *VcfRecord) Dirty() bool { return r.dirty }

// MarkDirty flags the record as modified. Call this after mutating Info() or a
// Sample() in place; the record-level mutators below do it for you.
func (r *VcfRecord) MarkDirty() { r.dirty = true }

// SetChrom replaces the CHROM column.
func (r *VcfRecord) SetChrom(chrom string) {
	r.Chrom = chrom
	r.dirty = true
}

// SubsetSamples rebuilds the record's sample columns to the given 0-based column
// indices (in order). An empty slice drops all sample (and FORMAT) columns. An
// out-of-range index emits ".". It ports the sample-column projection used by
// vcf-reorder/vcf-strip; serialize() then emits the kept samples.
func (r *VcfRecord) SubsetSamples(indices []int) {
	r.ensureFormat()
	raw := make([]string, 0, len(indices))
	parsed := make([]*Attributes, 0, len(indices))
	for _, idx := range indices {
		if idx >= 0 && idx < len(r.sampleRaw) {
			raw = append(raw, r.sampleRaw[idx])
			parsed = append(parsed, r.samples[idx])
		} else {
			raw = append(raw, missing)
			parsed = append(parsed, nil)
		}
	}
	r.sampleRaw = raw
	r.samples = parsed
	r.dirty = true
}

// GenotypeBases resolves the GT of sample i to ref/alt bases (e.g. "A/G"),
// using REF and the parsed ALT alleles. It returns ok=false when GT is absent,
// not diploid, has a missing allele, or references an unknown allele index. It
// ports vcf-sample-export's --gt conversion (output is always "/"-joined).
func (r *VcfRecord) GenotypeBases(i int) (string, bool) {
	s, err := r.Sample(i)
	if err != nil {
		return "", false
	}
	gt, ok := s.Get("GT")
	if !ok {
		return "", false
	}
	return resolveGenotypeBases(gt.String(), r.Ref, r.Alt())
}

func resolveGenotypeBases(gt, ref string, alts []string) (string, bool) {
	var parts []string
	switch {
	case strings.Contains(gt, "/"):
		parts = strings.Split(gt, "/")
	case strings.Contains(gt, "|"):
		parts = strings.Split(gt, "|")
	default:
		return "", false
	}
	if len(parts) != 2 || parts[0] == missing || parts[1] == missing {
		return "", false
	}
	a0, ok := alleleBase(parts[0], ref, alts)
	if !ok {
		return "", false
	}
	a1, ok := alleleBase(parts[1], ref, alts)
	if !ok {
		return "", false
	}
	return a0 + "/" + a1, true
}

func alleleBase(s, ref string, alts []string) (string, bool) {
	n, err := strconv.Atoi(s)
	if err != nil || n < 0 {
		return "", false
	}
	if n == 0 {
		return ref, true
	}
	if n-1 < len(alts) {
		return alts[n-1], true
	}
	return "", false
}

// SetID sets the ID column (use "" to clear it to ".").
func (r *VcfRecord) SetID(id string) {
	r.id = id
	r.idDone = true
	r.dirty = true
}

// ClearID clears the ID column to ".".
func (r *VcfRecord) ClearID() { r.SetID("") }

// AddInfo sets an INFO key to a value.
func (r *VcfRecord) AddInfo(key, value string) {
	r.Info().Set(key, value)
	r.dirty = true
}

// AddInfoFlag sets an INFO flag (bare key, no value).
func (r *VcfRecord) AddInfoFlag(key string) {
	r.Info().SetFlag(key)
	r.dirty = true
}

// AddFormat sets a FORMAT key to a value for the given sample. The key is
// appended to that sample's fields (and becomes part of the FORMAT column on
// write). Annotators that add a per-sample field call this for every sample.
func (r *VcfRecord) AddFormat(sampleIdx int, key, value string) error {
	s, err := r.Sample(sampleIdx)
	if err != nil {
		return err
	}
	s.Set(key, value)
	r.dirty = true
	return nil
}

// AddFilter appends a FILTER code to the record (marking it as failing that
// filter). A record with no codes is PASS.
func (r *VcfRecord) AddFilter(code string) {
	_ = r.Filters() // ensure parsed
	r.filters = append(r.filters, code)
	r.dirty = true
}

// ClearFilters removes all FILTER codes, returning the record to PASS.
func (r *VcfRecord) ClearFilters() {
	r.filters = nil
	r.filtDone = true
	r.dirty = true
}

// SetFilters replaces the FILTER codes (an empty slice clears to PASS).
func (r *VcfRecord) SetFilters(codes []string) {
	r.filters = append([]string(nil), codes...)
	r.filtDone = true
	r.dirty = true
}

// RetainFilters keeps only the FILTER codes for which keep returns true,
// removing the rest in place. Unlike SetFilters it preserves the PASS-vs-"."
// distinction: a PASS record (no codes) is left untouched, while a record that
// already carried codes stays non-PASS (rendering ".") even when emptied. This
// matches ngsutilsj's parse-time filter stripping (vcf-strip). Marks dirty when
// the record carried codes.
func (r *VcfRecord) RetainFilters(keep func(code string) bool) {
	cur := r.Filters()
	if cur == nil {
		return
	}
	kept := cur[:0]
	for _, f := range cur {
		if keep(f) {
			kept = append(kept, f)
		}
	}
	r.filters = kept
	r.filtDone = true
	r.dirty = true
}

// filterField renders the FILTER column: PASS when no filters, "." when the
// column was explicitly missing, else the codes joined by ";". Ports the FILTER
// logic in VCFRecord.write.
func (r *VcfRecord) filterField() string {
	f := r.Filters()
	if f == nil {
		return pass
	}
	if len(f) == 0 {
		return missing
	}
	return strings.Join(f, ";")
}

// serialize reconstructs the full VCF line from the parsed model. It is used for
// records that have been modified; it ports ngsutilsj VCFRecord.write (ID/ALT/
// QUAL/FILTER/INFO rules, FORMAT keys derived from the first sample, per-sample
// trailing-missing trim).
func (r *VcfRecord) serialize() string {
	cols := make([]string, 0, 8+r.NumSamples()+1)
	cols = append(cols, r.Chrom, strconv.Itoa(r.Pos))
	if id := r.ID(); id == "" {
		cols = append(cols, missing)
	} else {
		cols = append(cols, id)
	}
	cols = append(cols, r.Ref)
	if alt := r.Alt(); len(alt) == 0 {
		cols = append(cols, missing)
	} else {
		cols = append(cols, strings.Join(alt, ","))
	}
	if q := r.Qual(); q == -1 {
		cols = append(cols, missing)
	} else {
		cols = append(cols, formatFloat(q))
	}
	cols = append(cols, r.filterField())
	cols = append(cols, r.Info().infoString())

	if n := r.NumSamples(); n > 0 {
		s0, _ := r.Sample(0)
		formatKeys := s0.Keys()
		cols = append(cols, strings.Join(formatKeys, ":"))
		for i := 0; i < n; i++ {
			si, _ := r.Sample(i)
			cols = append(cols, si.formatString(formatKeys))
		}
	}
	return strings.Join(cols, "\t")
}

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
