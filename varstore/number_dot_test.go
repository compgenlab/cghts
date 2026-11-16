package varstore

import (
	"path/filepath"
	"strings"
	"testing"
)

// Number=. fields are captured whole, as text, whatever their declared Type: a
// variable-length value has no per-ALT split to trust and no single number to
// keep. These pin the coercion, the guard against bypassing it, and that the
// verbatim text survives the store.

func TestNumberDotIsCapturedAsText(t *testing.T) {
	for _, typ := range []string{"String", "Float", "Integer"} {
		f := InfoFieldFor("XS", typ, ".")
		if f.Type != InfoString || f.Number != "." || f.Column != "info_xs" {
			t.Errorf("InfoFieldFor(XS, %s, .) = %+v, want String/./info_xs", typ, f)
		}
		if err := ValidateInfo([]InfoField{f}); err != nil {
			t.Errorf("InfoFieldFor(%s, .) does not validate: %v", typ, err)
		}
		g := FormatFieldFor("XF", typ, ".")
		if g.Type != InfoString || g.Number != "." || g.Column != "fmt_xf" {
			t.Errorf("FormatFieldFor(XF, %s, .) = %+v, want String/./fmt_xf", typ, g)
		}
		if err := ValidateFormat([]FormatField{g}); err != nil {
			t.Errorf("FormatFieldFor(%s, .) does not validate: %v", typ, err)
		}
	}
	// Other cardinalities keep their declared type.
	if f := InfoFieldFor("R2", "Float", "1"); f.Type != InfoFloat {
		t.Errorf("Number=1 Float became %s", f.Type)
	}
	if f := InfoFieldFor("IMP", "Flag", "0"); f.Type != InfoFlag {
		t.Errorf("a Flag became %s", f.Type)
	}
}

// A caller building the struct by hand must not get a numeric Number=. column,
// which would silently store only the rows that carry one value.
func TestNumberDotMustBeText(t *testing.T) {
	err := ValidateInfo([]InfoField{{Name: "XS", Type: InfoFloat, Number: "."}})
	if err == nil || !strings.Contains(err.Error(), "XS") || !strings.Contains(err.Error(), "verbatim") {
		t.Errorf("a Float Number=. INFO field: %v", err)
	}
	err = ValidateFormat([]FormatField{{Name: "XF", Type: InfoInteger, Number: "."}})
	if err == nil || !strings.Contains(err.Error(), "XF") {
		t.Errorf("an Integer Number=. FORMAT field: %v", err)
	}
	for _, n := range []string{"R", "G", "2"} {
		if ValidateInfo([]InfoField{InfoFieldFor("X", "Float", n)}) == nil {
			t.Errorf("Number=%s INFO was accepted", n)
		}
	}
}

func TestNumberDotRoundTripsVerbatim(t *testing.T) {
	dir := filepath.Join(t.TempDir(), "store")
	w := infoWriter(t, dir, []InfoField{InfoFieldFor("XS", "Float", "."), InfoFieldFor("IMPUTED", "String", ".")})
	for _, alt := range []string{"T", "G"} {
		if err := w.WriteSiteInfo(Site{Chrom: "chr1", Pos: 100, Ref: "A", Alt: alt, AN: 2},
			map[string]any{"XS": "0.91,0.88", "IMPUTED": "beagle,panel5"}); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.WriteSiteInfo(Site{Chrom: "chr1", Pos: 200, Ref: "C", Alt: "G", AN: 2}, nil); err != nil {
		t.Fatal(err)
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
	s, err := OpenParquet(dir)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	for _, f := range s.InfoFields() {
		if f.Type != InfoString || f.Number != "." {
			t.Errorf("manifest records %+v, want String/.", f)
		}
	}
	var n int
	if err := s.SitesInfo(func(site Site, info InfoRow) bool {
		xs, okX := info.String("info_xs")
		imp, okI := info.String("info_imputed")
		switch site.Pos {
		case 100:
			n++
			if xs != "0.91,0.88" || imp != "beagle,panel5" {
				t.Errorf("chr1:100 %s: xs=%q imputed=%q, want the record's full values", site.Alt, xs, imp)
			}
		case 200:
			if okX || okI {
				t.Errorf("chr1:200 carried no values but reads xs=%q imputed=%q", xs, imp)
			}
		}
		return true
	}); err != nil {
		t.Fatal(err)
	}
	if n != 2 {
		t.Errorf("read %d rows at chr1:100, want one per ALT", n)
	}

	fdir := filepath.Join(t.TempDir(), "fstore")
	fw := formatWriter(t, fdir, []FormatField{FormatFieldFor("XF", "Integer", ".")})
	if err := fw.WriteSite(Site{Chrom: "chr1", Pos: 100, Ref: "G", Alt: "A", AN: 4, NCalled: 2}); err != nil {
		t.Fatal(err)
	}
	if err := fw.WriteCallFormat(Call{SampleID: "S1", Chrom: "chr1", Pos: 100, Ref: "G", Alt: "A",
		GT: "0/1", DP: 40, ADRef: 20, ADAlt: 20, GQ: 99}, map[string]any{"XF": "3,1,4"}); err != nil {
		t.Fatal(err)
	}
	if err := fw.Finish(); err != nil {
		t.Fatal(err)
	}
	fs, err := OpenParquet(fdir)
	if err != nil {
		t.Fatal(err)
	}
	defer fs.Close()
	var got string
	if err := fs.CallsFormat(func(c Call, f FormatRow) bool {
		got, _ = f.String("fmt_xf")
		return true
	}); err != nil {
		t.Fatal(err)
	}
	if got != "3,1,4" {
		t.Errorf("fmt_xf = %q, want the call's full value", got)
	}
}
