package vcf

import (
	"reflect"
	"strings"
	"testing"
)

// The header mutators all maintain two structures at once -- a map for lookup and
// a slice for order -- so every one of them can leave the two disagreeing. These
// assert on both, plus on Lines(), which is what actually reaches the file.

func parseHeader(t *testing.T, lines ...string) *VcfHeader {
	t.Helper()
	var b strings.Builder
	b.WriteString("##fileformat=VCFv4.2\n")
	for _, l := range lines {
		b.WriteString(l + "\n")
	}
	b.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	r, err := NewVcfReader(strings.NewReader(b.String()))
	if err != nil {
		t.Fatal(err)
	}
	h, err := r.Header()
	if err != nil {
		t.Fatal(err)
	}
	return h
}

func headerText(t *testing.T, h *VcfHeader) string {
	t.Helper()
	return strings.Join(h.Lines(), "\n") + "\n"
}

func TestHeaderAddFilterAndFilterDef(t *testing.T) {
	h := NewVcfHeader()
	if _, ok := h.FilterDef("lowqual"); ok {
		t.Fatal("an empty header reported a FILTER definition")
	}

	h.AddFilter(&FilterDef{ID: "lowqual", Description: "Low quality"})
	d, ok := h.FilterDef("lowqual")
	if !ok {
		t.Fatal("FilterDef did not find the definition just added")
	}
	if d.Description != "Low quality" {
		t.Errorf("Description = %q, want %q", d.Description, "Low quality")
	}
	if got := h.FilterIDs(); !reflect.DeepEqual(got, []string{"lowqual"}) {
		t.Errorf("FilterIDs() = %v, want [lowqual]", got)
	}
	// A definition with no OrigLine is rendered, with the description quoted.
	if want := `##FILTER=<ID=lowqual,Description="Low quality">`; !strings.Contains(headerText(t, h), want) {
		t.Errorf("header does not contain %q:\n%s", want, headerText(t, h))
	}
}

// Re-adding an existing ID replaces the definition without appending a second
// entry to the order -- otherwise the header would carry the ID twice and a
// reader would see a duplicate declaration.
func TestHeaderAddFilterReplacesInPlace(t *testing.T) {
	h := NewVcfHeader()
	h.AddFilter(&FilterDef{ID: "a", Description: "first"})
	h.AddFilter(&FilterDef{ID: "b", Description: "other"})
	h.AddFilter(&FilterDef{ID: "a", Description: "second"})

	if got := h.FilterIDs(); !reflect.DeepEqual(got, []string{"a", "b"}) {
		t.Errorf("FilterIDs() = %v, want [a b] -- a replacement should keep its position", got)
	}
	d, _ := h.FilterDef("a")
	if d.Description != "second" {
		t.Errorf("Description = %q, want the replacement", d.Description)
	}
	if n := strings.Count(headerText(t, h), "ID=a,"); n != 1 {
		t.Errorf("filter a appears %d times in the header, want 1", n)
	}
}

func TestHeaderAddAlt(t *testing.T) {
	h := NewVcfHeader()
	h.AddAlt(&AltDef{ID: "NON_REF", Description: "Represents any possible alternative allele"})
	h.AddAlt(&AltDef{ID: "DEL"}) // no description

	if got := h.AltIDs(); !reflect.DeepEqual(got, []string{"NON_REF", "DEL"}) {
		t.Errorf("AltIDs() = %v, want [NON_REF DEL]", got)
	}
	if d, ok := h.AltDef("NON_REF"); !ok || d.ID != "NON_REF" {
		t.Errorf("AltDef(NON_REF) = (%+v, %v)", d, ok)
	}
	if _, ok := h.AltDef("MISSING"); ok {
		t.Error("AltDef reported an absent ID as present")
	}

	text := headerText(t, h)
	// gVCF detection keys on this line, so its exact spelling matters.
	if want := `##ALT=<ID=NON_REF,Description="Represents any possible alternative allele">`; !strings.Contains(text, want) {
		t.Errorf("header does not contain %q:\n%s", want, text)
	}
	// A description-less definition renders without an empty Description="".
	if want := "##ALT=<ID=DEL>"; !strings.Contains(text, want) {
		t.Errorf("header does not contain %q:\n%s", want, text)
	}
}

func TestHeaderAddAltReplacesInPlace(t *testing.T) {
	h := NewVcfHeader()
	h.AddAlt(&AltDef{ID: "NON_REF", Description: "first"})
	h.AddAlt(&AltDef{ID: "NON_REF", Description: "second"})
	if got := h.AltIDs(); !reflect.DeepEqual(got, []string{"NON_REF"}) {
		t.Errorf("AltIDs() = %v, want one entry", got)
	}
	d, _ := h.AltDef("NON_REF")
	if d.Description != "second" {
		t.Errorf("Description = %q, want the replacement", d.Description)
	}
}

// AddLine is the escape hatch for anything the parser does not model -- a
// ##reference, a ##source, a command-line provenance stamp -- and it must survive
// verbatim, since the point is to preserve what cghts does not understand.
func TestHeaderAddLine(t *testing.T) {
	h := parseHeader(t, `##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">`)
	h.AddLine("##reference=file:///ref/GRCh38.fa")
	h.AddLine(`##cgkit_vcf-clearfilterCommand=cgkit vcf-clearfilter in.vcf`)

	want := []string{"##reference=file:///ref/GRCh38.fa", "##cgkit_vcf-clearfilterCommand=cgkit vcf-clearfilter in.vcf"}
	if got := h.OtherLines(); !reflect.DeepEqual(got, want) {
		t.Errorf("OtherLines() = %v, want %v", got, want)
	}

	lines := h.Lines()
	// Added lines come last, after the structured definitions.
	if n := len(lines); n < 2 || lines[n-2] != want[0] || lines[n-1] != want[1] {
		t.Errorf("added lines are not at the end of Lines():\n%s", strings.Join(lines, "\n"))
	}
	// Unrecognized lines are preserved on the way in as well as on the way out.
	round := parseHeader(t, `##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">`, want[0])
	if got := round.OtherLines(); !reflect.DeepEqual(got, []string{want[0]}) {
		t.Errorf("a parsed unrecognized line came back as %v, want %v", got, []string{want[0]})
	}
}

func TestHeaderRemoveContig(t *testing.T) {
	h := parseHeader(t,
		"##contig=<ID=chr1,length=248956422>",
		"##contig=<ID=chr2,length=242193529>",
		"##contig=<ID=chrM,length=16569>",
	)
	if got := h.ContigNames(); !reflect.DeepEqual(got, []string{"chr1", "chr2", "chrM"}) {
		t.Fatalf("ContigNames() = %v", got)
	}

	h.RemoveContig("chr2")
	if got := h.ContigNames(); !reflect.DeepEqual(got, []string{"chr1", "chrM"}) {
		t.Errorf("ContigNames() = %v, want [chr1 chrM]", got)
	}
	if _, ok := h.ContigDef("chr2"); ok {
		t.Error("ContigDef still resolves a removed contig -- the map and the order have drifted")
	}
	if text := headerText(t, h); strings.Contains(text, "ID=chr2") {
		t.Errorf("the removed contig is still in the header:\n%s", text)
	}
	// The survivors keep their definitions, not just their names.
	if d, ok := h.ContigDef("chr1"); !ok || d.Length != 248956422 {
		t.Errorf("ContigDef(chr1) = (%+v, %v) after removing a different contig", d, ok)
	}

	// Removing something absent is a no-op, not a panic and not a truncation.
	h.RemoveContig("chrNope")
	if got := h.ContigNames(); !reflect.DeepEqual(got, []string{"chr1", "chrM"}) {
		t.Errorf("removing an absent contig changed the list to %v", got)
	}

	// Removing the last one leaves an empty list rather than a stale entry.
	h.RemoveContig("chr1")
	h.RemoveContig("chrM")
	if got := h.ContigNames(); len(got) != 0 {
		t.Errorf("ContigNames() = %v, want empty", got)
	}
	if strings.Contains(headerText(t, h), "##contig=") {
		t.Error("a contig line survived removing every contig")
	}
}

// The removers used to compact the order slice in place -- the same slice
// ContigNames, InfoIDs, FormatIDs, FilterIDs and Attributes.Keys hand out. So
// the obvious way to drop a set of entries, ranging over the accessor and
// removing as you go, skipped every second one: range fixes the length up front
// while the elements slide down past the cursor. It failed silently, leaving a
// header that still declared contigs the caller had asked to remove.
func TestRemovalDoesNotDisturbAnAccessorSlice(t *testing.T) {
	t.Run("contigs", func(t *testing.T) {
		h := parseHeader(t,
			"##contig=<ID=chr1,length=100>",
			"##contig=<ID=chr2,length=100>",
			"##contig=<ID=chr3,length=100>",
			"##contig=<ID=chr4,length=100>",
		)
		for _, c := range h.ContigNames() {
			h.RemoveContig(c)
		}
		if got := h.ContigNames(); len(got) != 0 {
			t.Errorf("ContigNames() = %v after removing every contig, want empty", got)
		}
		if text := headerText(t, h); strings.Contains(text, "##contig=") {
			t.Errorf("contig lines survived:\n%s", text)
		}
	})

	t.Run("info", func(t *testing.T) {
		h := parseHeader(t,
			`##INFO=<ID=A,Number=1,Type=Integer,Description="d">`,
			`##INFO=<ID=B,Number=1,Type=Integer,Description="d">`,
			`##INFO=<ID=C,Number=1,Type=Integer,Description="d">`,
			`##INFO=<ID=D,Number=1,Type=Integer,Description="d">`,
		)
		for _, id := range h.InfoIDs() {
			h.RemoveInfo(id)
		}
		if got := h.InfoIDs(); len(got) != 0 {
			t.Errorf("InfoIDs() = %v after removing every field, want empty", got)
		}
	})

	// A slice taken before a removal keeps reading what it read at the time --
	// it is a snapshot, not a view that rewrites itself underneath the caller.
	t.Run("an earlier slice is not rewritten", func(t *testing.T) {
		h := parseHeader(t,
			"##contig=<ID=chr1,length=100>",
			"##contig=<ID=chr2,length=100>",
			"##contig=<ID=chr3,length=100>",
		)
		before := h.ContigNames()
		h.RemoveContig("chr1")
		if !reflect.DeepEqual(before, []string{"chr1", "chr2", "chr3"}) {
			t.Errorf("a slice taken before the removal now reads %v", before)
		}
	})

	t.Run("attributes", func(t *testing.T) {
		a := newAttributes()
		for _, k := range []string{"A", "B", "C", "D"} {
			a.Set(k, "1")
		}
		for _, k := range a.Keys() {
			a.Remove(k)
		}
		if got := a.Keys(); len(got) != 0 {
			t.Errorf("Keys() = %v after removing every key, want empty", got)
		}
		if got := a.infoString(); got != "." {
			t.Errorf("infoString() = %q, want %q", got, ".")
		}
	})
}

func TestHeaderSetFileDate(t *testing.T) {
	h := parseHeader(t, "##fileDate=20200101", "##source=somecaller")
	h.SetFileDate("20260805")

	other := h.OtherLines()
	var dates []string
	for _, l := range other {
		if strings.HasPrefix(l, "##fileDate=") {
			dates = append(dates, l)
		}
	}
	// Exactly one -- two ##fileDate lines is not something a reader has to
	// resolve, so replacing rather than appending is the contract.
	if !reflect.DeepEqual(dates, []string{"##fileDate=20260805"}) {
		t.Errorf("fileDate lines = %v, want exactly the new one", dates)
	}
	if !strings.Contains(strings.Join(other, "\n"), "##source=somecaller") {
		t.Errorf("SetFileDate dropped an unrelated line: %v", other)
	}
	// The replacement is appended, so it moves to the end of the other lines
	// even when the original was first. Order among ## lines is not meaningful,
	// but pinning it keeps output byte-stable.
	if other[len(other)-1] != "##fileDate=20260805" {
		t.Errorf("OtherLines() = %v, want the new date last", other)
	}

	// A slice taken before the call is a snapshot: filtering otherLines in place
	// used to rewrite it under the caller.
	snapshot := parseHeader(t, "##fileDate=20200101", "##source=somecaller")
	was := snapshot.OtherLines()
	snapshot.SetFileDate("20260805")
	if !reflect.DeepEqual(was, []string{"##fileDate=20200101", "##source=somecaller"}) {
		t.Errorf("a slice taken before SetFileDate now reads %v", was)
	}

	// On a header with no date, it is simply added.
	fresh := parseHeader(t, "##source=somecaller")
	fresh.SetFileDate("20260805")
	if got := fresh.OtherLines(); !reflect.DeepEqual(got, []string{"##source=somecaller", "##fileDate=20260805"}) {
		t.Errorf("OtherLines() = %v", got)
	}

	// And on an entirely empty header, where otherLines is nil.
	empty := NewVcfHeader()
	empty.SetFileDate("20260805")
	if got := empty.OtherLines(); !reflect.DeepEqual(got, []string{"##fileDate=20260805"}) {
		t.Errorf("OtherLines() on a fresh header = %v", got)
	}
}

func TestHeaderFormatIDs(t *testing.T) {
	h := parseHeader(t,
		`##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`,
		`##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">`,
		`##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">`,
	)
	// Header order, not the order INFO and FORMAT were interleaved in the file.
	if got := h.FormatIDs(); !reflect.DeepEqual(got, []string{"GT", "AD"}) {
		t.Errorf("FormatIDs() = %v, want [GT AD]", got)
	}
	if got := h.InfoIDs(); !reflect.DeepEqual(got, []string{"DP"}) {
		t.Errorf("InfoIDs() = %v, want [DP]", got)
	}

	h.AddFormat(&AnnotationDef{ID: "GQ", Number: "1", Type: "Integer", Description: "Genotype quality"})
	if got := h.FormatIDs(); !reflect.DeepEqual(got, []string{"GT", "AD", "GQ"}) {
		t.Errorf("FormatIDs() after AddFormat = %v, want [GT AD GQ]", got)
	}
	h.RemoveFormat("AD")
	if got := h.FormatIDs(); !reflect.DeepEqual(got, []string{"GT", "GQ"}) {
		t.Errorf("FormatIDs() after RemoveFormat = %v, want [GT GQ]", got)
	}
	if _, ok := h.FormatDef("AD"); ok {
		t.Error("FormatDef still resolves a removed ID")
	}
}

// The glob matchers are what back vcf-strip's --keep-info/--keep-format, so a
// pattern matching too much silently retains annotations the caller asked to drop.
func TestHeaderMatchIDs(t *testing.T) {
	h := parseHeader(t,
		`##INFO=<ID=DP,Number=1,Type=Integer,Description="d">`,
		`##INFO=<ID=AF,Number=A,Type=Float,Description="d">`,
		`##INFO=<ID=AC,Number=A,Type=Integer,Description="d">`,
		`##INFO=<ID=ANN,Number=.,Type=String,Description="d">`,
		`##FORMAT=<ID=GT,Number=1,Type=String,Description="d">`,
		`##FORMAT=<ID=AD,Number=R,Type=Integer,Description="d">`,
		`##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="d">`,
	)

	for _, tc := range []struct {
		glob string
		want []string
	}{
		{"*", []string{"DP", "AF", "AC", "ANN"}},
		{"A*", []string{"AF", "AC", "ANN"}},
		{"A?", []string{"AF", "AC"}},
		{"DP", []string{"DP"}}, // exact, no implied wildcard
		{"??", []string{"DP", "AF", "AC"}},
		{"Z*", nil},
	} {
		if got := h.MatchInfoIDs(tc.glob); !reflect.DeepEqual(got, tc.want) {
			t.Errorf("MatchInfoIDs(%q) = %v, want %v", tc.glob, got, tc.want)
		}
	}

	for _, tc := range []struct {
		glob string
		want []string
	}{
		{"*", []string{"GT", "AD", "MIN_DP"}},
		{"*DP", []string{"MIN_DP"}}, // DP is an INFO field here, not a FORMAT one
		{"?D", []string{"AD"}},
		{"MIN_*", []string{"MIN_DP"}},
		{"Z*", nil},
	} {
		if got := h.MatchFormatIDs(tc.glob); !reflect.DeepEqual(got, tc.want) {
			t.Errorf("MatchFormatIDs(%q) = %v, want %v", tc.glob, got, tc.want)
		}
	}
}

// MetaLines is Lines() minus the fileformat line: the definition block a caller
// splices into a header it is assembling itself. Emitting the fileformat line
// twice would produce a file bcftools rejects.
func TestHeaderMetaLines(t *testing.T) {
	h := parseHeader(t,
		`##INFO=<ID=DP,Number=1,Type=Integer,Description="d">`,
		`##FILTER=<ID=lowqual,Description="d">`,
		"##contig=<ID=chr1,length=100>",
		"##source=caller",
	)

	meta := h.MetaLines()
	for _, l := range meta {
		if strings.HasPrefix(l, "##fileformat=") {
			t.Errorf("MetaLines() contains the fileformat line: %q", l)
		}
		if strings.HasPrefix(l, "#CHROM") {
			t.Errorf("MetaLines() contains the #CHROM line: %q", l)
		}
	}
	if got, want := len(meta), len(h.Lines())-1; got != want {
		t.Errorf("MetaLines() has %d lines, want %d (Lines() minus fileformat)", got, want)
	}
	// It is Lines() in the same order, so a caller can rely on INFO before
	// FILTER before contig.
	if !reflect.DeepEqual(meta, h.Lines()[1:]) {
		t.Errorf("MetaLines() = %v, want Lines()[1:] = %v", meta, h.Lines()[1:])
	}

	// A header with nothing but a fileformat line yields no meta lines at all,
	// rather than one empty string.
	if got := NewVcfHeader().MetaLines(); len(got) != 0 {
		t.Errorf("MetaLines() on a bare header = %v, want empty", got)
	}
}
