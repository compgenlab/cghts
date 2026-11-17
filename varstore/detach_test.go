package varstore

import (
	"path/filepath"
	"strings"
	"testing"
	"unsafe"
)

// inside reports whether s's bytes lie within src's backing array -- that is,
// whether holding s keeps all of src alive.
func inside(s, src string) bool {
	if len(s) == 0 || len(src) == 0 {
		return false
	}
	p := uintptr(unsafe.Pointer(unsafe.StringData(s)))
	lo := uintptr(unsafe.Pointer(unsafe.StringData(src)))
	return p >= lo && p < lo+uintptr(len(src))
}

// A buffered row must not keep its source line alive. The VCF reader hands out
// CHROM/REF/ALT and FORMAT values as substrings of the record's line, which in a
// wide cohort is megabytes; 8,192 buffered rows each pinning one held gigabytes,
// and OOM-killed a conversion at 8 GB. Each check fails without detach.go.
func TestBufferedRowsDoNotPinTheirSource(t *testing.T) {
	// A stand-in for a wide record line.
	line := "chr1\t100\t.\tA\tG\t50\tPASS\tXS=0.9,0.8\tGT:PID\t0/1:p1" + strings.Repeat("\t0/0:.", 1<<16)
	at := func(sub string) string {
		i := strings.Index(line, sub)
		if i < 0 {
			t.Fatalf("%q not in the fixture line", sub)
		}
		return line[i : i+len(sub)]
	}
	chrom, ref, alt := at("chr1"), at("\tA\t")[1:2], at("\tG\t")[1:2]
	xs, gt, pid := at("0.9,0.8"), at("0/1"), at("p1")

	dir := filepath.Join(t.TempDir(), "store")
	w, err := NewWriter(dir, WriterOpts{
		Samples: []string{"S1"}, MinDP: 10,
		Info:   []InfoField{InfoFieldFor("XS", "Float", ".")},
		Format: []FormatField{FormatFieldFor("PID", "String", "1")},
	})
	if err != nil {
		t.Fatal(err)
	}
	if err := w.WriteSiteInfo(Site{Chrom: chrom, Pos: 100, Ref: ref, Alt: alt, AN: 2},
		map[string]any{"XS": xs}); err != nil {
		t.Fatal(err)
	}
	if err := w.WriteCallFormat(Call{SampleID: "S1", Chrom: chrom, Pos: 100, Ref: ref, Alt: alt, GT: gt},
		map[string]any{"PID": pid}); err != nil {
		t.Fatal(err)
	}
	if err := w.WriteRegion(CalledSiteRun{SampleID: "S1", Chrom: chrom, Start: 100, End: 100, NSites: 1}); err != nil {
		t.Fatal(err)
	}

	if len(w.siteRows) != 1 || len(w.callRows) != 1 || len(w.catalogRunBuf) != 1 {
		t.Fatalf("expected one buffered row per table, have %d/%d/%d",
			len(w.siteRows), len(w.callRows), len(w.catalogRunBuf))
	}
	check := func(what string, v any) {
		t.Helper()
		var s string
		switch x := v.(type) {
		case string:
			s = x
		case *string:
			s = *x
		default:
			t.Fatalf("%s is %T", what, v)
		}
		if inside(s, line) {
			t.Errorf("buffered %s %q still points into its source line", what, s)
		}
	}
	if line[11:12] != "A" || line[13:14] != "G" {
		t.Fatal("fixture offsets are wrong")
	}
	site := w.siteRows[0].(map[string]any)
	for _, k := range []string{"chrom", "ref", "alt", "info_xs"} {
		check("site "+k, site[k])
	}
	call := w.callRows[0].(map[string]any)
	for _, k := range []string{"chrom", "ref", "alt", "gt", "fmt_pid"} {
		check("call "+k, call[k])
	}
	check("run chrom", w.catalogRunBuf[0].Chrom)
	w.Discard()
}

// The typed (no capture) path buffers structs rather than maps; same rule.
func TestBufferedStructRowsDoNotPinTheirSource(t *testing.T) {
	line := "chr1\t100\t.\tA\tG" + strings.Repeat("\t0/0", 1<<16)
	chrom, ref, alt := line[0:4], line[11:12], line[13:14] // "chr1", "A", "G"
	w, err := NewWriter(filepath.Join(t.TempDir(), "store"), WriterOpts{Samples: []string{"S1"}, MinDP: 10})
	if err != nil {
		t.Fatal(err)
	}
	if err := w.WriteSite(Site{Chrom: chrom, Pos: 100, Ref: ref, Alt: alt, AN: 2}); err != nil {
		t.Fatal(err)
	}
	if err := w.WriteCall(Call{SampleID: "S1", Chrom: chrom, Pos: 100, Ref: ref, Alt: alt, GT: "0/1"}); err != nil {
		t.Fatal(err)
	}
	s, c := w.siteBuf[0], w.callBuf[0]
	for what, v := range map[string]string{"site chrom": s.Chrom, "site ref": s.Ref, "site alt": s.Alt,
		"call chrom": c.Chrom, "call ref": c.Ref, "call alt": c.Alt} {
		if inside(v, line) {
			t.Errorf("buffered %s %q still points into its source line", what, v)
		}
	}
	w.Discard()
}
