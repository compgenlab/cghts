package vcf

import (
	"strconv"
	"strings"
)

// ToUCSC converts a chromosome name to UCSC form: a "chr" prefix is added, and
// the mitochondrial contig "MT" becomes "chrM". Names already starting with
// "chr" are returned unchanged. It ports vcf-chrfix's --ucsc mapping.
func ToUCSC(name string) string {
	if strings.HasPrefix(name, "chr") {
		return name
	}
	if name == "MT" {
		return "chrM"
	}
	return "chr" + name
}

// ToEnsembl converts a chromosome name to Ensembl form: a "chr" prefix is
// stripped, and "chrM" becomes "MT". Names without a "chr" prefix are returned
// unchanged. It ports vcf-chrfix's --ensembl mapping.
func ToEnsembl(name string) string {
	if !strings.HasPrefix(name, "chr") {
		return name
	}
	if name == "chrM" {
		return "MT"
	}
	return name[len("chr"):]
}

// PrimaryHumanContigs is the set of primary human contigs (1-22, X, Y, M, MT),
// in Ensembl naming. vcf-chrfix --primary-human keeps a contig when its
// "chr"-stripped name is in this set.
var PrimaryHumanContigs = map[string]bool{
	"1": true, "2": true, "3": true, "4": true, "5": true, "6": true,
	"7": true, "8": true, "9": true, "10": true, "11": true, "12": true,
	"13": true, "14": true, "15": true, "16": true, "17": true, "18": true,
	"19": true, "20": true, "21": true, "22": true,
	"X": true, "Y": true, "M": true, "MT": true,
}

// IsPrimaryHuman reports whether chrom (UCSC or Ensembl named) is a primary
// human contig. Naming is matched case-insensitively, so "chr1", "Chr1" and
// "CHR1" agree; see stripChrPrefix.
func IsPrimaryHuman(chrom string) bool {
	return PrimaryHumanContigs[stripChrPrefix(chrom)]
}

// stripChrPrefix upper-cases a contig name and then removes a leading "CHR".
//
// The order matters, and used to be the other way round. TrimPrefix is
// case-sensitive, so stripping "chr" before upper-casing left "Chr1" and "CHR1"
// untouched: they upper-cased to "CHR1", missed the contig table, fell through
// to the NCBI accession parser and came back as not primary at all. Both
// spellings occur in the wild, and the failure is silent -- CanonicalContig
// returns ok=false, callers fall back to exact-name matching, and a query
// against a differently spelled file quietly matches nothing.
func stripChrPrefix(name string) string {
	return strings.TrimPrefix(strings.ToUpper(name), "CHR")
}

// CanonicalContig returns a canonical chromosome key (Ensembl-style: "1".."22",
// "X", "Y", "MT") for a human primary contig named in UCSC ("chr1"/"chrM"),
// Ensembl ("1"/"MT"/"M"), or NCBI RefSeq ("NC_000001.11") form. ok is false for
// names outside the primary set (unplaced scaffolds, alt loci, NCBI NT_/NW_,
// non-human accessions); callers fall back to exact-name matching. It underpins
// [ContigConverter] and the vcf-annotate --auto-convert option.
func CanonicalContig(name string) (key string, ok bool) {
	// UCSC/Ensembl: strip a "chr" prefix, normalize the mitochondrion to "MT".
	// Shared with IsPrimaryHuman so the two cannot disagree about what "primary"
	// means for a given spelling.
	s := stripChrPrefix(name)
	if s == "M" {
		s = "MT"
	}
	if PrimaryHumanContigs[s] {
		return s, true
	}
	return ncbiCanonical(name)
}

// ncbiCanonical maps a human NCBI RefSeq accession to its canonical key. The
// primary chromosomes are NC_000001..NC_000024 (1-22, then 23=X, 24=Y) and the
// rCRS mitochondrion is NC_012920. The accession version (".11") is ignored, so
// the mapping is assembly-independent (GRCh37 and GRCh38 both resolve).
func ncbiCanonical(name string) (string, bool) {
	const prefix = "NC_"
	if !strings.HasPrefix(name, prefix) {
		return "", false
	}
	digits := name[len(prefix):]
	if i := strings.IndexByte(digits, '.'); i >= 0 {
		digits = digits[:i] // drop the version suffix
	}
	n, err := strconv.Atoi(digits)
	if err != nil {
		return "", false
	}
	switch {
	case n == 12920:
		return "MT", true
	case n >= 1 && n <= 22:
		return strconv.Itoa(n), true
	case n == 23:
		return "X", true
	case n == 24:
		return "Y", true
	}
	return "", false
}
