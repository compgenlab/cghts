package vcf

import (
	"reflect"
	"testing"
)

// Attributes backs both INFO and FORMAT, and its whole contract is that key order
// is insertion order -- serialize() renders INFO in that order and derives the
// FORMAT column from it, so a map-ordered implementation would produce a
// different file on every run.

func TestAttributesKeepInsertionOrder(t *testing.T) {
	a := newAttributes()
	a.Set("DP", "30")
	a.SetFlag("DB")
	a.SetValue("AF", AttrValue{raw: "0.5"})

	want := []string{"DP", "DB", "AF"}
	if got := a.Keys(); !reflect.DeepEqual(got, want) {
		t.Errorf("Keys() = %v, want %v", got, want)
	}
	if got := a.infoString(); got != "DP=30;DB;AF=0.5" {
		t.Errorf("infoString() = %q, want %q", got, "DP=30;DB;AF=0.5")
	}
}

// Replacing a key must not move it. INFO column order is stable output, and an
// annotator that overwrites DP should not shuffle every field after it.
func TestAttributesSetReplacesInPlace(t *testing.T) {
	a := newAttributes()
	a.Set("DP", "30")
	a.Set("AF", "0.5")
	a.Set("DP", "99")

	if got := a.Keys(); !reflect.DeepEqual(got, []string{"DP", "AF"}) {
		t.Errorf("Keys() = %v, want [DP AF] -- a replaced key should keep its position", got)
	}
	v, ok := a.Get("DP")
	if !ok || v.String() != "99" {
		t.Errorf("Get(DP) = (%q, %v), want (99, true)", v.String(), ok)
	}
}

// The three setters differ only in what they store, and the difference is
// load-bearing: a flag renders as a bare key, "." renders as itself, and the two
// are not the same claim.
func TestAttributesSetterFlavors(t *testing.T) {
	a := newAttributes()
	a.Set("VAL", "7")
	a.SetFlag("FLAG")
	a.SetValue("MISS", AttrValue{raw: missing})

	val, _ := a.Get("VAL")
	if val.IsEmpty() || val.IsMissing() {
		t.Errorf("a plain value read as empty/missing: %+v", val)
	}
	flag, _ := a.Get("FLAG")
	if !flag.IsEmpty() {
		t.Error("SetFlag did not produce an empty (bare-flag) value")
	}
	if flag.IsMissing() {
		t.Error("a bare flag is not the missing marker")
	}
	miss, _ := a.Get("MISS")
	if !miss.IsMissing() {
		t.Error("SetValue did not preserve the missing marker")
	}
	if miss.IsEmpty() {
		t.Error("the missing marker is not an empty value")
	}

	// A present-but-missing key is still present -- which is the distinction
	// Get's boolean exists to report.
	if !a.Contains("MISS") {
		t.Error("a missing-valued key should still be present")
	}
	if got := a.infoString(); got != "VAL=7;FLAG;MISS=." {
		t.Errorf("infoString() = %q, want %q", got, "VAL=7;FLAG;MISS=.")
	}
}

func TestAttributesFindKeys(t *testing.T) {
	a := newAttributes()
	for _, k := range []string{"AD", "AF", "DP", "MIN_DP", "GQ"} {
		a.Set(k, "1")
	}

	for _, tc := range []struct {
		glob string
		want []string
	}{
		{"*", []string{"AD", "AF", "DP", "MIN_DP", "GQ"}},
		{"A*", []string{"AD", "AF"}},
		{"A?", []string{"AD", "AF"}},
		// A trailing glob is not implied: "DP" matches only DP, not MIN_DP.
		{"DP", []string{"DP"}},
		{"*DP", []string{"DP", "MIN_DP"}},
		{"?D", []string{"AD"}},
		{"ZZ*", nil},
	} {
		if got := a.FindKeys(tc.glob); !reflect.DeepEqual(got, tc.want) {
			t.Errorf("FindKeys(%q) = %v, want %v", tc.glob, got, tc.want)
		}
	}
}

// FindKeys reports in insertion order, not glob order or map order, so a caller
// removing matched keys walks them in the order they appear in the column.
func TestAttributesFindKeysFollowsInsertionOrder(t *testing.T) {
	a := newAttributes()
	a.Set("ZED", "1")
	a.Set("ABLE", "2")
	if got := a.FindKeys("*"); !reflect.DeepEqual(got, []string{"ZED", "ABLE"}) {
		t.Errorf("FindKeys(*) = %v, want [ZED ABLE]", got)
	}
}

func TestAttrValueIsMissing(t *testing.T) {
	for _, tc := range []struct {
		raw  string
		want bool
	}{
		{".", true},
		{"", false},   // a bare flag is empty, not missing
		{"0", false},  // zero is a real value
		{"..", false}, // only the exact marker
		{"./.", false},
	} {
		if got := (AttrValue{raw: tc.raw}).IsMissing(); got != tc.want {
			t.Errorf("AttrValue{%q}.IsMissing() = %v, want %v", tc.raw, got, tc.want)
		}
	}
}
