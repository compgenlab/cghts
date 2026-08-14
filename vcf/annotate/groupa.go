package annotate

import (
	"strconv"
	"strings"

	"github.com/compgenlab/cghts/vcf"
)

// AutoID sets the ID column to a canonical variant identifier, one per alt,
// ';'-joined:
//
//	1-115256529-T-C
//	X-9720345-CCAGA-C
//
// The chromosome loses its "chr" and the four fields are joined with hyphens.
// That is the form gnomAD, ClinGen and most variant portals use, so an id
// written here can be pasted into one of them and an id from one of them
// recognized here — which is the point of synthesizing one at all.
//
// It was chrom_pos_ref_alt before, which matched nothing in particular.
type AutoID struct{ closeNoop }

// NewAutoID returns an AutoID annotator (--auto-id).
func NewAutoID() *AutoID { return &AutoID{} }

// SetupHeader adds no header definitions: this writes the ID column, which
// every VCF already has.
func (a *AutoID) SetupHeader(*vcf.VcfHeader) error { return nil }

// Annotate sets the record ID.
//
// Multiple identifiers are ';'-joined, which is the ID column's own separator
// for exactly that.
func (a *AutoID) Annotate(rec *vcf.VcfRecord) error {
	alt := rec.Alt()
	ids := make([]string, len(alt))
	for i, al := range alt {
		ids[i] = VariantID(rec.Chrom, rec.Pos, rec.Ref, al)
	}
	rec.SetID(strings.Join(ids, ";"))
	return nil
}

// VariantID renders one variant as "1-115256529-T-C".
//
// Exported because it is a format rather than an implementation detail:
// anything producing or recognizing the same identifier should come through
// here instead of rebuilding the string and drifting from it.
//
// Nothing is normalized beyond the chromosome. REF and ALT are written exactly
// as they arrived and the position is the VCF position — no left-alignment, no
// trimming of a shared first base. An identifier that quietly disagreed with the
// row it labels would be worse than none.
func VariantID(chrom string, pos int, ref, alt string) string {
	return TrimChrom(chrom) + "-" + strconv.Itoa(pos) + "-" + ref + "-" + alt
}

// TrimChrom drops a leading "chr" from a contig name.
//
// Only that exact prefix, in any case, and only with something after it: the
// convention is "chr1" and "chrUn_GL000220v1", and a contig genuinely named
// "chr" is not one to render as the empty string. A name that never had the
// prefix — "1", "GL000220.1" — comes back untouched, which is what makes this
// safe on an assembly that does not use it.
func TrimChrom(chrom string) string {
	if len(chrom) > 3 && strings.EqualFold(chrom[:3], "chr") {
		return chrom[3:]
	}
	return chrom
}

// ConstantTag adds a fixed INFO flag or key=value to every record.
type ConstantTag struct {
	closeNoop
	key   string
	value string
	flag  bool
}

// NewConstantFlag returns a ConstantTag that adds a bare INFO flag (--tag KEY).
func NewConstantFlag(key string) *ConstantTag { return &ConstantTag{key: key, flag: true} }

// NewConstantTag returns a ConstantTag that adds INFO key=value (--tag KEY:VALUE).
func NewConstantTag(key, value string) *ConstantTag {
	return &ConstantTag{key: key, value: value}
}

// SetupHeader declares the INFO field.
func (a *ConstantTag) SetupHeader(h *vcf.VcfHeader) error {
	if a.flag {
		h.AddInfo(infoDef(a.key, "0", "Flag", a.key))
	} else {
		h.AddInfo(infoDef(a.key, ".", "String", a.key))
	}
	return nil
}

// Annotate adds the constant tag.
func (a *ConstantTag) Annotate(rec *vcf.VcfRecord) error {
	if a.flag {
		rec.AddInfoFlag(a.key)
	} else {
		rec.AddInfo(a.key, a.value)
	}
	return nil
}

// Indel flags insertions/deletions and records their lengths.
type Indel struct {
	namedFields
	closeNoop
}

// NewIndel returns an Indel annotator (--indel).
func NewIndel() *Indel { return &Indel{} }

// SetupHeader declares the indel INFO fields. (ngsutilsj registers these as
// FORMAT defs, which is a bug — the values go into INFO; this package registers
// them correctly as INFO.)
func (a *Indel) SetupHeader(h *vcf.VcfHeader) error {
	h.AddInfo(infoDef(a.key(IndelInsert), "0", "Flag", "Variant is an insertion"))
	h.AddInfo(infoDef(a.key(IndelDelete), "0", "Flag", "Variant is an deletion"))
	h.AddInfo(infoDef(a.key(IndelInsLen), "1", "Integer", "Insertion length"))
	h.AddInfo(infoDef(a.key(IndelDelLen), "1", "Integer", "Deletion length"))
	h.AddInfo(infoDef(a.key(IndelIndelLen), "1", "Integer", "In-del length"))
	return nil
}

// Annotate adds the indel flags/lengths.
func (a *Indel) Annotate(rec *vcf.VcfRecord) error {
	insert, deletion := false, false
	insLen, delLen := 0, 0
	if len(rec.Ref) > 1 {
		deletion = true
		delLen = len(rec.Ref) - 1
	}
	for _, alt := range rec.Alt() {
		if len(alt) > 1 {
			insert = true
			if l := len(alt) - 1; l > insLen {
				insLen = l
			}
		}
	}
	if insert {
		rec.AddInfoFlag(a.key(IndelInsert))
		rec.AddInfo(a.key(IndelInsLen), strconv.Itoa(insLen))
		rec.AddInfo(a.key(IndelIndelLen), strconv.Itoa(insLen))
	}
	if deletion {
		rec.AddInfoFlag(a.key(IndelDelete))
		rec.AddInfo(a.key(IndelDelLen), strconv.Itoa(delLen))
		rec.AddInfo(a.key(IndelIndelLen), "-"+strconv.Itoa(delLen))
	}
	return nil
}

// TsTv classifies SNVs as transition (TS) or transversion (TV).
type TsTv struct {
	namedFields
	closeNoop
}

// NewTsTv returns a TsTv annotator (--tstv).
func NewTsTv() *TsTv { return &TsTv{} }

// SetupHeader declares the CG_TSTV INFO field.
func (a *TsTv) SetupHeader(h *vcf.VcfHeader) error {
	h.AddInfo(infoDef(a.key(TsTvField), "1", "String", "Is the variant and transition (TS) or transversion (TV), skips all multi-variants and indels"))
	return nil
}

// Annotate adds CG_TSTV for single-base biallelic SNVs.
func (a *TsTv) Annotate(rec *vcf.VcfRecord) error {
	switch rec.CalcTsTv() {
	case -1:
		rec.AddInfo(a.key(TsTvField), "TS")
	case 1:
		rec.AddInfo(a.key(TsTvField), "TV")
	}
	return nil
}

// Dosage computes per-alt allele dosage from each sample's GT.
type Dosage struct {
	namedFields
	closeNoop
}

// NewDosage returns a Dosage annotator (--dosage).
func NewDosage() *Dosage { return &Dosage{} }

// SetupHeader declares the CG_DS FORMAT field.
func (a *Dosage) SetupHeader(h *vcf.VcfHeader) error {
	h.AddFormat(formatDef(a.key(DosageField), "A", "Integer", "Convert GT to dosage value (0, 1, 2)"))
	return nil
}

// Annotate adds CG_DS to every sample.
func (a *Dosage) Annotate(rec *vcf.VcfRecord) error {
	nalt := len(rec.Alt())
	for i := 0; i < rec.NumSamples(); i++ {
		s, err := rec.Sample(i)
		if err != nil {
			return err
		}
		gt, ok := s.Get("GT")
		if !ok {
			if err := rec.AddFormat(i, a.key(DosageField), "."); err != nil {
				return err
			}
			continue
		}
		gts := strings.FieldsFunc(gt.String(), func(r rune) bool { return r == '/' || r == '|' })
		outs := make([]string, 0, nalt)
		for altNum := 1; altNum <= nalt; altNum++ {
			ds := 0
			target := strconv.Itoa(altNum)
			for _, g := range gts {
				if g == target {
					ds++
				}
			}
			outs = append(outs, strconv.Itoa(ds))
		}
		if err := rec.AddFormat(i, a.key(DosageField), strings.Join(outs, ",")); err != nil {
			return err
		}
	}
	return nil
}
