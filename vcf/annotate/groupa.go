package annotate

import (
	"strconv"
	"strings"

	"github.com/compgenlab/cghts/vcf"
)

// AutoID sets the ID column to chrom_pos_ref_alt (one per alt, ';'-joined).
type AutoID struct{ closeNoop }

// NewAutoID returns an AutoID annotator (--auto-id).
func NewAutoID() *AutoID { return &AutoID{} }

// SetupHeader adds no header definitions.
func (a *AutoID) SetupHeader(*vcf.VcfHeader) error { return nil }

// Annotate sets the record ID.
func (a *AutoID) Annotate(rec *vcf.VcfRecord) error {
	alt := rec.Alt()
	ids := make([]string, len(alt))
	for i, al := range alt {
		ids[i] = rec.Chrom + "_" + strconv.Itoa(rec.Pos) + "_" + rec.Ref + "_" + al
	}
	rec.SetID(strings.Join(ids, ";"))
	return nil
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
type Indel struct{ closeNoop }

// NewIndel returns an Indel annotator (--indel).
func NewIndel() *Indel { return &Indel{} }

// SetupHeader declares the indel INFO fields. (ngsutilsj registers these as
// FORMAT defs, which is a bug — the values go into INFO; this package registers
// them correctly as INFO.)
func (a *Indel) SetupHeader(h *vcf.VcfHeader) error {
	h.AddInfo(infoDef("CG_INSERT", "0", "Flag", "Variant is an insertion"))
	h.AddInfo(infoDef("CG_DELETE", "0", "Flag", "Variant is an deletion"))
	h.AddInfo(infoDef("CG_INSLEN", "1", "Integer", "Insertion length"))
	h.AddInfo(infoDef("CG_DELLEN", "1", "Integer", "Deletion length"))
	h.AddInfo(infoDef("CG_INDELLEN", "1", "Integer", "In-del length"))
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
		rec.AddInfoFlag("CG_INSERT")
		rec.AddInfo("CG_INSLEN", strconv.Itoa(insLen))
		rec.AddInfo("CG_INDELLEN", strconv.Itoa(insLen))
	}
	if deletion {
		rec.AddInfoFlag("CG_DELETE")
		rec.AddInfo("CG_DELLEN", strconv.Itoa(delLen))
		rec.AddInfo("CG_INDELLEN", "-"+strconv.Itoa(delLen))
	}
	return nil
}

// TsTv classifies SNVs as transition (TS) or transversion (TV).
type TsTv struct{ closeNoop }

// NewTsTv returns a TsTv annotator (--tstv).
func NewTsTv() *TsTv { return &TsTv{} }

// SetupHeader declares the CG_TSTV INFO field.
func (a *TsTv) SetupHeader(h *vcf.VcfHeader) error {
	h.AddInfo(infoDef("CG_TSTV", "1", "String", "Is the variant and transition (TS) or transversion (TV), skips all multi-variants and indels"))
	return nil
}

// Annotate adds CG_TSTV for single-base biallelic SNVs.
func (a *TsTv) Annotate(rec *vcf.VcfRecord) error {
	switch rec.CalcTsTv() {
	case -1:
		rec.AddInfo("CG_TSTV", "TS")
	case 1:
		rec.AddInfo("CG_TSTV", "TV")
	}
	return nil
}

// Dosage computes per-alt allele dosage from each sample's GT.
type Dosage struct{ closeNoop }

// NewDosage returns a Dosage annotator (--dosage).
func NewDosage() *Dosage { return &Dosage{} }

// SetupHeader declares the CG_DS FORMAT field.
func (a *Dosage) SetupHeader(h *vcf.VcfHeader) error {
	h.AddFormat(formatDef("CG_DS", "A", "Integer", "Convert GT to dosage value (0, 1, 2)"))
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
			if err := rec.AddFormat(i, "CG_DS", "."); err != nil {
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
		if err := rec.AddFormat(i, "CG_DS", strings.Join(outs, ",")); err != nil {
			return err
		}
	}
	return nil
}
