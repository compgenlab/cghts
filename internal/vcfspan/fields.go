package vcfspan

import "strings"

// VCF column indices, 0-based, for the fixed fields that bear on the span.
const (
	colRef    = 3
	colAlt    = 4
	colInfo   = 7
	colFormat = 8
)

// LineFields adapts a VCF data line already split on tabs. This is what the tabix
// index writer and reader have: they parse lines themselves and must not depend on
// the vcf package, which imports tabix rather than the other way round.
//
// A truncated line is tolerated rather than rejected -- the span is a hint for
// indexing, and a malformed line is the caller's problem to report, not a reason to
// panic here.
type LineFields []string

// FieldsEnd is End over a split line, returning beg+1 when the line is too short
// to carry a REF at all.
func FieldsEnd(fields []string, beg int) int {
	if len(fields) <= colRef {
		return beg + 1
	}
	return End(LineFields(fields), beg)
}

func (f LineFields) Ref() string {
	if len(f) <= colRef {
		return ""
	}
	return f[colRef]
}

func (f LineFields) Alts() []string {
	if len(f) <= colAlt {
		return nil
	}
	return strings.Split(f[colAlt], ",")
}

func (f LineFields) Info(key string) (string, bool) {
	if len(f) <= colInfo {
		return "", false
	}
	info := f[colInfo]
	for len(info) > 0 {
		var kv string
		if i := strings.IndexByte(info, ';'); i >= 0 {
			kv, info = info[:i], info[i+1:]
		} else {
			kv, info = info, ""
		}
		if v, ok := strings.CutPrefix(kv, key+"="); ok {
			return v, true
		}
	}
	return "", false
}

func (f LineFields) SampleValues(key string) []string {
	if len(f) <= colFormat {
		return nil
	}
	pos := -1
	for i, k := range strings.Split(f[colFormat], ":") {
		if k == key {
			pos = i
			break
		}
	}
	if pos < 0 {
		return nil
	}
	out := make([]string, 0, len(f)-colFormat-1)
	for _, sample := range f[colFormat+1:] {
		if v := subfield(sample, pos); v != "" {
			out = append(out, v)
		}
	}
	return out
}

// subfield returns the i-th colon-separated component, or "" when there are fewer
// -- a sample column may legally be truncated where trailing fields are missing.
func subfield(s string, i int) string {
	for ; i > 0; i-- {
		j := strings.IndexByte(s, ':')
		if j < 0 {
			return ""
		}
		s = s[j+1:]
	}
	if j := strings.IndexByte(s, ':'); j >= 0 {
		s = s[:j]
	}
	return s
}
