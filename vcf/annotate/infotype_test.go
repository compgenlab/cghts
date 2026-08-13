package annotate

import (
	"strings"
	"testing"

	"github.com/compgenlab/cghts/vcf"
)

// A caller declaring an INFO key can declare its type.
//
// The annotators that copy a value out of a source see text, and text cannot say
// whether it is a score or a label. Without a way to declare it, every copied
// value was String — so a numeric field reached consumers as a string, and a
// downstream tool sorted "10" before "9" on a file that was otherwise correct.
func TestResolveInfoType(t *testing.T) {
	for _, tc := range []struct {
		typ      InfoType
		isNumber bool
		want     string
		why      string
	}{
		{"", false, "String", "unset and not numeric is the old default"},
		{"", true, "Float", "unset falls back to the isNumber flag"},
		{TypeFloat, false, "Float", "a declared type needs no flag"},
		{TypeInteger, false, "Integer", "Integer was unreachable before"},
		{TypeString, true, "String", "an explicit type beats the flag"},
		{TypeFlag, false, "Flag", ""},
		{TypeCharacter, false, "Character", ""},
		{"float", false, "Float", "case-insensitive, for types read from config"},
		{"INTEGER", false, "Integer", ""},
	} {
		got, err := resolveInfoType(tc.typ, tc.isNumber)
		if err != nil {
			t.Errorf("resolveInfoType(%q, %v): %v", tc.typ, tc.isNumber, err)
			continue
		}
		if got != tc.want {
			t.Errorf("resolveInfoType(%q, %v) = %q, want %q — %s",
				tc.typ, tc.isNumber, got, tc.want, tc.why)
		}
	}
}

// An unrecognised type is refused rather than quietly becoming String.
//
// A caller who wrote "Double" meant something. Declaring the opposite of what
// they asked for is how a wrong type reaches a file nobody reads again.
func TestAnUnknownInfoTypeIsRefused(t *testing.T) {
	_, err := resolveInfoType("Double", false)
	if err == nil {
		t.Fatal("an unknown type was accepted")
	}
	if !strings.Contains(err.Error(), "Double") || !strings.Contains(err.Error(), "Float") {
		t.Errorf("the error names neither the mistake nor the alternatives: %v", err)
	}
}

// The declared type reaches the ##INFO line a VCF annotation writes.
//
// Asserted through SetupHeader rather than against the resolver, because the
// resolver being right is no use if the annotator does not consult it — which is
// exactly what was wrong: VcfAnnotation hardcoded "String" and had nowhere to
// take a type from.
func TestAVcfAnnotationDeclaresItsType(t *testing.T) {
	a := &VcfAnnotation{opts: VcfOptions{
		Name: "gnomad_af", Field: "AF", Filename: "g.vcf.gz", Type: TypeFloat,
	}}
	h := vcf.NewVcfHeader()
	if err := a.SetupHeader(h); err != nil {
		t.Fatal(err)
	}
	def, ok := h.InfoDef("gnomad_af")
	if !ok {
		t.Fatal("no ##INFO def was added")
	}
	if def.Type != "Float" {
		t.Errorf("Type = %q, want Float — the declared type did not reach the header", def.Type)
	}
}

// And an undeclared one is still String, so existing callers are unaffected.
func TestAnUndeclaredVcfTypeIsStillString(t *testing.T) {
	a := &VcfAnnotation{opts: VcfOptions{Name: "clinvar_sig", Field: "CLNSIG", Filename: "c.vcf.gz"}}
	h := vcf.NewVcfHeader()
	if err := a.SetupHeader(h); err != nil {
		t.Fatal(err)
	}
	def, _ := h.InfoDef("clinvar_sig")
	if def.Type != "String" {
		t.Errorf("Type = %q, want String", def.Type)
	}
}

// A bad type fails at header setup, which is before any record is written.
func TestABadTypeFailsBeforeAnythingIsWritten(t *testing.T) {
	a := &VcfAnnotation{opts: VcfOptions{
		Name: "x", Field: "X", Filename: "f.vcf.gz", Type: "Numeric",
	}}
	if err := a.SetupHeader(vcf.NewVcfHeader()); err == nil {
		t.Fatal("SetupHeader accepted an unknown type")
	}
}
