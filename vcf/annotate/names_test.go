package annotate

import (
	"strings"
	"testing"

	"github.com/compgenlab/cghts/vcf"
)

// An annotator writes its fields under the names the caller chose.
//
// The caller names its output and the annotator supplies the value. That held
// for the annotators that copy a value from a source and not for the ones that
// compute it, so an annotation a caller had named "tstv" arrived as CG_TSTV —
// and callers bridged the gap by translating afterwards, which works only where
// there is an afterwards. A streaming VCF has none.
func TestAComputedFieldTakesTheCallersName(t *testing.T) {
	a := NewTsTv()
	a.SetFieldNames(FieldNames{TsTvField: "tstv"})

	h := vcf.NewVcfHeader()
	if err := a.SetupHeader(h); err != nil {
		t.Fatal(err)
	}
	if _, ok := h.InfoDef("tstv"); !ok {
		t.Error("the header does not declare the caller's name")
	}
	if _, ok := h.InfoDef("CG_TSTV"); ok {
		t.Error("the header still declares the default name as well")
	}

	// And the value is written under it, not merely declared.
	rec := vcf.NewRecord("chr1", 100, "T", []string{"C"})
	if err := a.Annotate(rec); err != nil {
		t.Fatal(err)
	}
	if v, ok := rec.Info().Get("tstv"); !ok || v.String() != "TS" {
		t.Errorf("value not written under the caller's name: %v", rec.Info().Keys())
	}
	if rec.Info().Contains("CG_TSTV") {
		t.Error("the value was also written under the default name")
	}
}

// An unnamed field keeps the name it has always had.
//
// Every existing caller passes nothing, so their output must not move.
func TestAnUnnamedFieldKeepsItsDefault(t *testing.T) {
	a := NewTsTv()
	h := vcf.NewVcfHeader()
	if err := a.SetupHeader(h); err != nil {
		t.Fatal(err)
	}
	if _, ok := h.InfoDef("CG_TSTV"); !ok {
		t.Error("an unnamed field did not keep CG_TSTV")
	}
}

// A partial map renames only what it names.
//
// An annotator writing five fields should not require a caller who wants one of
// them renamed to name the other four.
func TestRenamingIsPerField(t *testing.T) {
	a := NewIndel()
	a.SetFieldNames(FieldNames{IndelInsert: "is_insertion"})

	h := vcf.NewVcfHeader()
	if err := a.SetupHeader(h); err != nil {
		t.Fatal(err)
	}
	if _, ok := h.InfoDef("is_insertion"); !ok {
		t.Error("the named field was not renamed")
	}
	for _, keep := range []string{"CG_DELETE", "CG_INSLEN", "CG_DELLEN", "CG_INDELLEN"} {
		if _, ok := h.InfoDef(keep); !ok {
			t.Errorf("%s moved even though it was not named", keep)
		}
	}
}

// An empty name is not a name.
//
// A caller building this map from a config file produces empty strings for the
// fields nobody named, and a field written under "" is a record no parser reads.
func TestAnEmptyNameFallsBackToTheDefault(t *testing.T) {
	a := NewTsTv()
	a.SetFieldNames(FieldNames{TsTvField: ""})
	h := vcf.NewVcfHeader()
	if err := a.SetupHeader(h); err != nil {
		t.Fatal(err)
	}
	if _, ok := h.InfoDef("CG_TSTV"); !ok {
		t.Error("an empty name should leave the default in place")
	}
}

// Every computing annotator can be renamed, so a caller need not learn which
// ones happen to support it.
//
// The list is spelled out rather than discovered, because the failure this
// guards is an annotator added later that quietly does not implement it — which
// no reflective test over the existing ones would catch either.
func TestEveryComputedAnnotatorIsRenameable(t *testing.T) {
	for name, a := range map[string]any{
		"tstv":          NewTsTv(),
		"indel":         NewIndel(),
		"dosage":        NewDosage(),
		"vaf":           NewVAF(),
		"minor_strand":  NewMinorStrand(),
		"fisher_sb":     NewFisherSB(),
		"vardist":       NewVariantDistance(),
		"copy_logratio": &CopyNumberLogRatio{},
	} {
		if _, ok := a.(FieldNamer); !ok {
			t.Errorf("%s cannot be renamed; its caller has no way to name its output", name)
		}
	}
}

// A GTF annotator names its fields the same way, and a named one ignores the
// prefix — the prefix decorates a default, and a caller who gave the whole key
// has not asked for it to be decorated.
func TestGtfFieldsTakeTheCallersNames(t *testing.T) {
	a := &GtfAnnotator{opts: GtfOptions{Names: FieldNames{GtfGeneSymbol: "gene_symbol"}},
		prefix: "GTF_"}
	if got := a.key(GtfGeneSymbol); got != "gene_symbol" {
		t.Errorf("named field = %q, want the caller's name undecorated", got)
	}
	if got := a.key(GtfGeneID); got != "GTF_GENEID" {
		t.Errorf("unnamed field = %q, want the prefixed default", got)
	}
}

var _ = strings.Contains
