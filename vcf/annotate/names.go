package annotate

// Naming the fields an annotator writes.
//
// An annotator that copies a value from a source is told the INFO key to write
// it under — see VcfOptions.Name — because only the caller knows what the field
// should be called. The annotators that *compute* their values were not: they
// wrote fixed names of their own, so an annotation a caller had named `tstv`
// arrived as CG_TSTV and one named `gene` arrived as GTF_GENE.
//
// That is a defect rather than a design. The caller names its output; the
// annotator supplies the value. It held for one half of this package and not the
// other, and callers bridged the gap by translating names afterwards — which
// works only where there is an afterwards, and a streaming VCF has none.
//
// So the computing annotators take a map. A map rather than a single name
// because most of them write more than one field: an indel annotator writes five,
// a GTF annotator seven. The keys are the logical fields, declared as constants
// below, and the values are what to call them in the output. An absent key keeps
// the historical name, so an existing caller sees no change at all.

// FieldNames maps an annotator's logical fields to the INFO or FORMAT keys they
// are written under. A nil map, or an absent key, keeps the default.
type FieldNames map[string]string

// nameOr returns the caller's name for a logical field, or the default.
//
// An empty value counts as absent: a caller building this map from a
// configuration file will produce empty strings for fields nobody named, and
// writing a field under "" would produce a record no parser accepts.
func (f FieldNames) nameOr(logical, dflt string) string {
	if f == nil {
		return dflt
	}
	if v := f[logical]; v != "" {
		return v
	}
	return dflt
}

// FieldNamer is implemented by annotators whose output field names can be
// chosen by the caller.
//
// A method rather than a constructor argument, matching [ContigMatcher]: the
// constructors are many and most callers want the defaults, so the ones that do
// not can say so without every call site growing a parameter. Set it before
// SetupHeader, which is where the names are first used.
type FieldNamer interface {
	SetFieldNames(FieldNames)
}

// The logical fields a GTF annotator writes.
//
// The values are the historical suffixes, so a caller already using them as
// names — which is the only way a GTF annotation could have worked, since the
// name had to match what the annotator wrote — needs no change.
const (
	GtfGeneSymbol = "GENE"      // gene name, e.g. TP53
	GtfGeneID     = "GENEID"    // gene identifier, e.g. ENSG00000141510
	GtfStrand     = "STRAND"    // + or -
	GtfBiotype    = "BIOTYPE"   // protein_coding, lncRNA, …
	GtfRegion     = "REGION"    // genic region: coding_exon, intron, UTR, …
	GtfCoding     = "CODING"    // overlapping coding gene names
	GtfNoncoding  = "NONCODING" // overlapping non-coding gene names
)

// The logical fields the indel annotator writes.
const (
	IndelInsert   = "INSERT"   // flag: the variant is an insertion
	IndelDelete   = "DELETE"   // flag: the variant is a deletion
	IndelInsLen   = "INSLEN"   // inserted bases
	IndelDelLen   = "DELLEN"   // deleted bases
	IndelIndelLen = "INDELLEN" // net length change
)

// The logical fields the flanking-sequence annotator writes.
const (
	FlankingSeq = "FLANKING"     // reference bases either side of the variant
	FlankingSub = "FLANKING_SUB" // the substitution the variant causes
)

// The logical field of each single-field annotator.
//
// Named for uniformity: a caller renaming several annotators should not have to
// remember which of them happen to write one field and so take a bare name.
const (
	TsTvField     = "TSTV"    // transition/transversion call
	VarDistField  = "VARDIST" // distance to the nearest other variant
	CopyLogRatio  = "CNLR"    // copy-number log2 ratio
	DosageField   = "DS"      // per-sample allele dosage
	VAFField      = "VAF"     // per-sample variant allele frequency
	MinorStrand   = "SBPCT"   // per-sample minor-strand percentage
	FisherStrandB = "FSB"     // per-sample Fisher strand bias
)

// namedFields is embedded by the computing annotators to give them settable
// output keys.
//
// Embedded rather than repeated because the rule is the same for all of them:
// the caller's name if there is one, otherwise the historical CG_ name. Nine
// copies of that would be nine chances for one to drift, and the drift would
// show up as a field written under a name its own header never declared.
type namedFields struct{ names FieldNames }

// SetFieldNames chooses the keys this annotator's fields are written under.
func (n *namedFields) SetFieldNames(f FieldNames) { n.names = f }

// key is the output key for a logical field. Every default is the logical name
// behind CG_, which is what these annotators have always written.
func (n *namedFields) key(logical string) string {
	return n.names.nameOr(logical, "CG_"+logical)
}
