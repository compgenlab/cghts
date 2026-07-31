// Package vcfspan computes how many reference bases a VCF record covers.
//
// It exists because two layers need the same answer from different shapes of the
// same data. The tabix index writer works on a raw line it has just split into
// fields and must stay independent of the vcf package, which sits above it in the
// import graph. Everything else holds a parsed record. Keeping one copy of the
// precedence rules matters more than either caller's convenience: if the index
// says a record spans 2 kb and a downstream tool says it spans one base, queries
// and output disagree and nothing reports it.
package vcfspan

import "strconv"

// Fields supplies the parts of a VCF record that bear on its reference span.
// Implementations adapt whatever they already have -- split line fields, or a
// parsed record -- so the rules below are written once.
type Fields interface {
	// Ref is the REF allele.
	Ref() string
	// Alts are the ALT alleles, comma-split. A bare "." may be included or
	// omitted; it affects nothing here.
	Alts() []string
	// Info returns an INFO value by key.
	Info(key string) (string, bool)
	// SampleValues returns one FORMAT value per sample for a key, skipping
	// samples that do not carry it. Only called for reference blocks.
	SampleValues(key string) []string
}

// End returns the 0-based, half-open end coordinate of a record whose 0-based
// start is beg.
//
// A VCF record covers more than one reference base far more often than its POS
// column suggests, and an index has to know that: a query for a region inside a
// long deletion or a gVCF reference block must still find the record spanning it.
// The tabix spec says so directly -- "for the VCF format, the end of a region
// equals POS plus the size of the deletion" -- and the `col_end = 0` in the VCF
// preset is not a claim that records are one base wide. That is a separate rule,
// which applies when col_beg *equals* col_end.
//
// Four sources contribute and the widest wins. Order and precedence follow
// htslib's tbx_parse1, because interoperating with indexes htslib wrote, and with
// the ones it reads, is the whole point:
//
//   - len(REF), the reference span of a plain record. The only source the tabix
//     spec itself documents.
//   - INFO/END, for symbolic alternates and GATK-style gVCF blocks. Ignored when
//     not greater than POS, matching htslib, since a malformed END must not shrink
//     the interval.
//   - INFO/SVLEN, only for symbolic alternates measured against the reference.
//     Taken as an absolute value: a deletion reports it negative, but the span it
//     occupies is positive either way.
//   - FORMAT/LEN, the VCF 4.5 per-sample reference-block length, only for <*> or
//     <NON_REF>.
//
// Widening beyond the truth is safe and narrowing is not, which settles the
// per-sample case: an index holds one interval per record, so FORMAT/LEN is taken
// as the maximum across samples. An interval that is too wide returns candidate
// records a caller then discards; one that is too narrow loses records with no way
// to notice.
func End(f Fields, beg int) int {
	span := len(f.Ref())
	if span == 0 {
		span = 1
	}

	// The alternates decide which of the remaining sources even apply.
	wantSVLen, wantFmtLen := false, false
	for _, alt := range f.Alts() {
		switch {
		case IsRefBlockAlt(alt):
			wantFmtLen = true
		case svLenMeasuredOnRef(alt):
			wantSVLen = true
		}
	}

	if wantSVLen {
		if v, ok := f.Info("SVLEN"); ok {
			for _, part := range splitComma(v) {
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

	if wantFmtLen {
		for _, v := range f.SampleValues("LEN") {
			if n, err := strconv.Atoi(v); err == nil && n > span {
				span = n
			}
		}
	}

	end := beg + span

	// INFO/END is an absolute 1-based inclusive coordinate, which makes it the
	// 0-based exclusive end unchanged.
	if v, ok := f.Info("END"); ok {
		if n, err := strconv.Atoi(v); err == nil && n > beg && n > end {
			end = n
		}
	}
	if end <= beg {
		end = beg + 1
	}
	return end
}

// IsRefBlockAlt reports whether an ALT marks a reference block rather than a
// variant: VCF 4.5's <*> or GATK's <NON_REF>. A bare "." is not included -- it
// means "no ALT reported", which a gVCF block written by some tools also uses, so
// callers that must recognise those check AltOrig separately.
func IsRefBlockAlt(alt string) bool {
	return alt == "<*>" || alt == "<NON_REF>"
}

// svLenMeasuredOnRef reports whether a symbolic alternate's SVLEN describes a span
// on the reference. Insertions and breakends do not: an <INS> occupies the bases
// named by REF however long the inserted sequence is, so letting SVLEN widen it
// would index a region the record does not cover.
func svLenMeasuredOnRef(alt string) bool {
	if len(alt) < 2 || alt[0] != '<' || alt[len(alt)-1] != '>' {
		return false
	}
	name := alt[1 : len(alt)-1]
	for i := 0; i < len(name); i++ {
		if name[i] == ':' {
			name = name[:i] // <DEL:ME:ALU> is a deletion
			break
		}
	}
	switch name {
	case "DEL", "DUP", "INV", "CNV":
		return true
	}
	return false
}

// splitComma splits on commas without allocating for the common single value.
func splitComma(s string) []string {
	n := 1
	for i := 0; i < len(s); i++ {
		if s[i] == ',' {
			n++
		}
	}
	if n == 1 {
		return []string{s}
	}
	out := make([]string, 0, n)
	start := 0
	for i := 0; i < len(s); i++ {
		if s[i] == ',' {
			out = append(out, s[start:i])
			start = i + 1
		}
	}
	return append(out, s[start:])
}
