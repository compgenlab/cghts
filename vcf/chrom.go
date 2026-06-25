package vcf

import "strings"

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
// human contig.
func IsPrimaryHuman(chrom string) bool {
	return PrimaryHumanContigs[strings.TrimPrefix(chrom, "chr")]
}
