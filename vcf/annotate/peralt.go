package annotate

import (
	"strings"

	"github.com/compgenlab/cghts/vcf"
)

// Attributing a copied value to the allele it came from.
//
// An exact-match annotation asks the source for the record's REF and one of its
// ALTs, so every value it finds belongs to a particular allele. Nothing recorded
// which: values were appended to one slice in the order the source happened to
// return them and joined with commas, which is the separator VCF uses for
// per-allele lists. The result read as a per-allele vector and was not one.
//
// The failure is silent and it is data corruption. For a record with ALTs A,C
// where the source lists C before A, the values came out in the source's order —
// so a reader assigned C's value to A. Where only the second ALT matched, a
// single value was written with nothing in front of it, and read as the first
// allele's. Both produce a well-formed file with the numbers on the wrong rows,
// which no parser can detect and no test that checks "is the value present"
// will catch.
//
// So a value is now placed at its allele's index and the field is padded to the
// ALT count with ".", which is what Number=A means. A position-matched
// annotation is not per-allele — it asked only about the locus — and keeps list
// semantics under Number=".".

// altMatches returns the indices of rec's ALT alleles that src also carries.
//
// Empty when REF differs or no allele is shared, which is the same condition
// altRefMatch reports; this says which rather than whether.
func altMatches(src, rec *vcf.VcfRecord) []int {
	if src.Ref != rec.Ref {
		return nil
	}
	srcAlts := src.Alt()
	var out []int
	for i, want := range rec.Alt() {
		for _, have := range srcAlts {
			if have == want {
				out = append(out, i)
				break
			}
		}
	}
	return out
}

// perAlleleField renders one value per ALT, padding with "." where an allele had
// none.
//
// Several source records can match the same allele — a duplicated entry, or two
// rows describing one variant — and their values are joined with "&" rather than
// a comma, because a comma is what separates the alleles. Dropping the extras
// would be the other option and it is worse: it loses data silently, and which
// value survived would depend on the order the source returned them.
//
// Reports false when no allele matched anything, so the caller writes no field
// at all rather than a row of dots.
func perAlleleField(byAlt [][]string, unique bool) (string, bool) {
	out := make([]string, len(byAlt))
	any := false
	for i, vals := range byAlt {
		if len(vals) == 0 {
			out[i] = "."
			continue
		}
		if unique {
			// Within the allele, not across the record. Deduplicating the whole
			// field would reorder it and break the alignment this exists to
			// keep.
			vals = sortedUnique(vals)
		}
		out[i] = strings.Join(vals, "&")
		any = true
	}
	if !any {
		return "", false
	}
	return strings.Join(out, ","), true
}

// infoNumberFor is the Number an annotation declares.
//
// "A" — one value per ALT — only when the match was exact, which is the case
// where a value belongs to an allele. A position match asks about the locus, so
// its values are a list of whatever the source had there: "." is the count VCF
// uses for that, and it is what the field actually contains.
//
// It used to be "1" in both cases, while several matches were written as a
// comma-separated list. That is invalid on its own terms — a strict reader is
// entitled to reject it — quite apart from the misattribution.
func infoNumberFor(exact bool) string {
	if exact {
		return "A"
	}
	return "."
}
