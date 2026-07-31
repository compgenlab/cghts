package tabix

import (
	"strconv"
	"strings"
)

// Preset values for the tabix index header's `format` field, matching htslib's
// TBX_* constants. The low 16 bits hold the preset; bit 0x10000 is the
// zero-based (BED) coordinate flag.
const (
	PresetGeneric int32 = 0
	PresetSAM     int32 = 1
	PresetVCF     int32 = 2

	presetMask   int32 = 0xffff
	zeroBasedBit int32 = 0x10000
)

// VCF column indices, 0-based, for the fixed fields we inspect.
const (
	vcfColRef    = 3
	vcfColAlt    = 4
	vcfColInfo   = 7
	vcfColFormat = 8
)

// vcfSpanEnd returns the 0-based, half-open end coordinate of a VCF data line
// whose 0-based start is beg.
//
// A VCF record covers more than one reference base far more often than its
// single POS column suggests, and the index has to know that: a query for a
// region inside a long deletion or a gVCF reference block must still find the
// record that spans it. The tabix spec says so directly -- "for the VCF format,
// the end of a region equals POS plus the size of the deletion" -- and the
// `col_end = 0` in the VCF preset is not a claim that records are one base wide.
// That is a separate rule, which applies when col_beg *equals* col_end.
//
// Four sources contribute, and the widest wins. Order and precedence follow
// htslib's tbx_parse1, because interoperating with indexes htslib wrote (and
// with the ones it reads) is the whole point:
//
//   - len(REF), the reference span of a plain record. This is the only source
//     the tabix spec itself documents.
//   - INFO/END, for symbolic alternates and GATK-style gVCF blocks. Ignored
//     when it is not greater than POS, matching htslib, since a malformed END
//     must not shrink the interval.
//   - INFO/SVLEN, consulted only for the symbolic alternates that measure their
//     length against the reference. Taken as an absolute value: a deletion
//     reports it negative, but the span it occupies is positive either way.
//   - FORMAT/LEN, the VCF 4.5 per-sample reference-block length, consulted only
//     when an alternate is <*> or <NON_REF>.
//
// Widening beyond the truth is safe and narrowing is not, which is what settles
// the per-sample case: an index holds one interval per record, so FORMAT/LEN is
// taken as the **maximum** across samples. An interval that is too wide returns
// candidate records a caller then discards; one that is too narrow loses records
// with no way to notice.
func vcfSpanEnd(fields []string, beg int) int {
	end := beg + 1
	if len(fields) <= vcfColRef {
		return end
	}

	refLen := len(fields[vcfColRef])
	if refLen == 0 {
		refLen = 1
	}
	span := refLen

	// The alternates decide which of the remaining sources even apply.
	wantSVLen, wantFmtLen := false, false
	if len(fields) > vcfColAlt {
		for _, alt := range strings.Split(fields[vcfColAlt], ",") {
			switch {
			case alt == "<*>" || alt == "<NON_REF>":
				wantFmtLen = true
			case svLenMeasuredOnRef(alt):
				wantSVLen = true
			}
		}
	}

	infoEnd := -1
	if len(fields) > vcfColInfo {
		info := fields[vcfColInfo]
		if v, ok := infoField(info, "END"); ok {
			if n, err := strconv.Atoi(v); err == nil {
				infoEnd = n
			}
		}
		if wantSVLen {
			if v, ok := infoField(info, "SVLEN"); ok {
				for _, part := range strings.Split(v, ",") {
					n, err := strconv.Atoi(part)
					if err != nil {
						continue
					}
					if n < 0 {
						n = -n
					}
					if n > span {
						span = n
					}
				}
			}
		}
	}

	if wantFmtLen && len(fields) > vcfColFormat {
		if pos := formatKeyIndex(fields[vcfColFormat], "LEN"); pos >= 0 {
			for _, sample := range fields[vcfColFormat+1:] {
				n, err := strconv.Atoi(subfield(sample, pos))
				if err == nil && n > span {
					span = n
				}
			}
		}
	}

	if beg+span > end {
		end = beg + span
	}
	// INFO/END is an absolute coordinate rather than a length, and is 1-based
	// inclusive, which makes it the 0-based exclusive end unchanged. htslib
	// discards it when it does not exceed POS rather than letting it narrow the
	// record.
	if infoEnd > beg && infoEnd > end {
		end = infoEnd
	}
	return end
}

// svLenMeasuredOnRef reports whether a symbolic alternate's SVLEN describes a
// span on the reference. Insertions and breakends do not: an <INS> occupies the
// bases named by REF however long the inserted sequence is, so letting SVLEN
// widen it would index a region the record does not cover.
func svLenMeasuredOnRef(alt string) bool {
	if len(alt) < 2 || alt[0] != '<' || alt[len(alt)-1] != '>' {
		return false
	}
	name := alt[1 : len(alt)-1]
	if i := strings.IndexByte(name, ':'); i >= 0 {
		name = name[:i] // <DEL:ME:ALU> is a deletion
	}
	switch name {
	case "DEL", "DUP", "INV", "CNV":
		return true
	}
	return false
}

// infoField returns the value of a semicolon-separated INFO key.
func infoField(info, key string) (string, bool) {
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

// formatKeyIndex returns the position of a key within a colon-separated FORMAT
// string, or -1 when absent.
func formatKeyIndex(format, key string) int {
	for i, k := range strings.Split(format, ":") {
		if k == key {
			return i
		}
	}
	return -1
}

// subfield returns the i-th colon-separated component, or "" if there are
// fewer -- a sample column may legally be truncated where trailing fields are
// missing.
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
