package varstore

import (
	"path/filepath"
	"testing"
)

// Captured INFO fields.
//
// The tests that matter here are the COMPATIBILITY ones, not the round trip. A
// new column that broke every existing reader would be a far worse bug than one
// that failed to store a number, and neither direction announces itself: an old
// reader on a new file and a new reader on an old file both fail silently or not
// at all.

func infoWriter(t *testing.T, dir string, fields []InfoField) *Writer {
	t.Helper()
	w, err := NewWriter(dir, WriterOpts{
		Samples: []string{"S1"}, MinDP: 10, Info: fields,
	})
	if err != nil {
		t.Fatalf("NewWriter: %v", err)
	}
	return w
}

var r2Field = InfoField{Name: "R2", Type: InfoFloat, Number: "1"}
var impField = InfoField{Name: "IMP", Type: InfoFlag, Number: "0"}

func TestInfoRoundTrips(t *testing.T) {
	dir := filepath.Join(t.TempDir(), "store")
	w := infoWriter(t, dir, []InfoField{r2Field, impField})

	if err := w.WriteSiteInfo(
		Site{Chrom: "chr1", Pos: 100, Ref: "A", Alt: "T", AC: 7, AN: 20, NCalled: 10},
		map[string]any{"R2": 0.87, "IMP": true}); err != nil {
		t.Fatal(err)
	}
	if err := w.WriteSiteInfo(
		Site{Chrom: "chr1", Pos: 200, Ref: "C", Alt: "G", AC: 3, AN: 20, NCalled: 10},
		map[string]any{"R2": 0.42, "IMP": false}); err != nil {
		t.Fatal(err)
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}

	s, err := OpenParquet(dir)
	if err != nil {
		t.Fatalf("OpenParquet: %v", err)
	}
	defer s.Close()

	if got := s.InfoFields(); len(got) != 2 {
		t.Fatalf("manifest recorded %d info fields, want 2: %+v", len(got), got)
	}

	type seen struct {
		pos int32
		r2  float64
		imp bool
	}
	var rows []seen
	err = s.SitesInfo(func(site Site, info InfoRow) bool {
		r2, ok := info.Float("info_r2")
		if !ok {
			t.Errorf("chr1:%d has no info_r2", site.Pos)
		}
		rows = append(rows, seen{site.Pos, r2, info.Flag("info_imp")})
		// The Site itself must survive the dynamic read intact.
		if site.Chrom != "chr1" || site.Ref == "" || site.AN != 20 {
			t.Errorf("site rebuilt wrong: %+v", site)
		}
		return true
	})
	if err != nil {
		t.Fatalf("SitesInfo: %v", err)
	}
	want := []seen{{100, 0.87, true}, {200, 0.42, false}}
	if len(rows) != len(want) {
		t.Fatalf("got %d rows, want %d", len(rows), len(want))
	}
	for i := range want {
		if rows[i] != want[i] {
			t.Errorf("row %d = %+v, want %+v", i, rows[i], want[i])
		}
	}
}

// An existing reader must not notice that the file grew columns.
//
// This is the whole premise of putting captured INFO in sites.parquet rather
// than in a sidecar table, so if it ever stops holding, the design is wrong
// rather than the test.
func TestTypedSiteReaderIgnoresCapturedColumns(t *testing.T) {
	dir := filepath.Join(t.TempDir(), "store")
	w := infoWriter(t, dir, []InfoField{r2Field})
	if err := w.WriteSiteInfo(
		Site{Chrom: "chr2", Pos: 500, Ref: "A", Alt: "T", AC: 1, AN: 8, NCalled: 4},
		map[string]any{"R2": 0.5}); err != nil {
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

	var n int
	err = s.Sites(func(site Site) bool {
		n++
		if site.Chrom != "chr2" || site.Pos != 500 || site.AC != 1 || site.AN != 8 {
			t.Errorf("typed read got %+v", site)
		}
		return true
	})
	if err != nil {
		t.Fatalf("Sites: %v", err)
	}
	if n != 1 {
		t.Fatalf("typed read saw %d sites, want 1", n)
	}
}

// And the reverse: reading a store written before capture existed.
func TestSitesInfoOnStoreWithoutCapture(t *testing.T) {
	dir := filepath.Join(t.TempDir(), "store")
	w := infoWriter(t, dir, nil)
	if err := w.WriteSite(Site{Chrom: "chr3", Pos: 7, Ref: "A", Alt: "C", AN: 4}); err != nil {
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

	if got := s.InfoFields(); len(got) != 0 {
		t.Errorf("store captured nothing but reports %+v", got)
	}
	var n int
	err = s.SitesInfo(func(site Site, info InfoRow) bool {
		n++
		if info.Present("info_r2") {
			t.Error("a store that captured nothing reports a value present")
		}
		if _, ok := info.Float("info_r2"); ok {
			t.Error("absent field read as ok")
		}
		return true
	})
	if err != nil {
		t.Fatalf("SitesInfo: %v", err)
	}
	if n != 1 {
		t.Fatalf("saw %d sites, want 1", n)
	}
}

// Absent and zero are different claims, and the optional column exists to keep
// them apart. A site the imputation program said nothing about is not a site
// with R2 = 0 -- one is "unknown", the other is "certainly worthless", and a
// filter for R2 >= 0.3 treats them identically only by accident.
func TestAbsentInfoIsNotZero(t *testing.T) {
	dir := filepath.Join(t.TempDir(), "store")
	w := infoWriter(t, dir, []InfoField{r2Field})

	// Site 1 has a value, site 2 has none, site 3 has a genuine zero.
	if err := w.WriteSiteInfo(Site{Chrom: "chr1", Pos: 1, Ref: "A", Alt: "T"},
		map[string]any{"R2": 0.9}); err != nil {
		t.Fatal(err)
	}
	if err := w.WriteSiteInfo(Site{Chrom: "chr1", Pos: 2, Ref: "A", Alt: "T"}, nil); err != nil {
		t.Fatal(err)
	}
	if err := w.WriteSiteInfo(Site{Chrom: "chr1", Pos: 3, Ref: "A", Alt: "T"},
		map[string]any{"R2": 0.0}); err != nil {
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

	present := map[int32]bool{}
	value := map[int32]float64{}
	if err := s.SitesInfo(func(site Site, info InfoRow) bool {
		present[site.Pos] = info.Present("info_r2")
		v, _ := info.Float("info_r2")
		value[site.Pos] = v
		return true
	}); err != nil {
		t.Fatal(err)
	}

	if !present[1] || value[1] != 0.9 {
		t.Errorf("site 1: present=%v value=%v, want true 0.9", present[1], value[1])
	}
	if present[2] {
		t.Error("site 2 had no R2 in the source but reads as present")
	}
	if !present[3] {
		t.Error("site 3 had R2=0 in the source and must read as present, not absent")
	}
	if value[3] != 0 {
		t.Errorf("site 3 value = %v, want 0", value[3])
	}
}

// The buffer of row maps is reused across batches. A field the previous occupant
// of a slot carried must not survive into a site that lacks it -- inheriting an
// R2 from 8,192 rows ago is wrong in a way that looks entirely plausible, since
// the value is a real number from a real site.
//
// Written to span more than one batch so the reuse actually happens.
func TestReusedRowBufferDoesNotLeakValues(t *testing.T) {
	dir := filepath.Join(t.TempDir(), "store")
	w := infoWriter(t, dir, []InfoField{r2Field, impField})

	const n = batchSize + 100
	for i := 0; i < n; i++ {
		site := Site{Chrom: "chr1", Pos: int32(i + 1), Ref: "A", Alt: "T"}
		// Every other site carries a value; the rest carry none.
		if i%2 == 0 {
			if err := w.WriteSiteInfo(site, map[string]any{"R2": 0.5, "IMP": true}); err != nil {
				t.Fatal(err)
			}
			continue
		}
		if err := w.WriteSiteInfo(site, nil); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}

	s, err := OpenParquet(dir)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	var leaked, seen int
	if err := s.SitesInfo(func(site Site, info InfoRow) bool {
		seen++
		odd := site.Pos%2 == 0 // Pos is i+1, so even Pos means odd i: no value.
		if odd && (info.Present("info_r2") || info.Flag("info_imp")) {
			leaked++
		}
		return true
	}); err != nil {
		t.Fatal(err)
	}
	if seen != n {
		t.Fatalf("read %d sites, wrote %d", seen, n)
	}
	if leaked > 0 {
		t.Errorf("%d sites inherited a value from a reused row buffer", leaked)
	}
}

func TestValidateInfoRefusals(t *testing.T) {
	cases := []struct {
		name  string
		field InfoField
		why   string
	}{
		{"per-allele-and-ref", InfoField{Name: "AD", Type: InfoInteger, Number: "R"},
			"Number=R has a reference value with nowhere to go"},
		{"variable-length", InfoField{Name: "X", Type: InfoFloat, Number: "."},
			"Number=. is variable per record"},
		{"per-genotype", InfoField{Name: "PL", Type: InfoInteger, Number: "G"},
			"Number=G is variable per record"},
		{"unknown-type", InfoField{Name: "Q", Type: InfoType("Character"), Number: "1"},
			"Character is not a storable type"},
		{"flag-with-values", InfoField{Name: "F", Type: InfoFlag, Number: "1"},
			"a Flag carries no values"},
		{"bad-key", InfoField{Name: "not a key", Type: InfoFloat, Number: "1"},
			"not a legal INFO key"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if err := ValidateInfo([]InfoField{c.field}); err == nil {
				t.Fatalf("accepted %+v, but %s", c.field, c.why)
			}
		})
	}
}

// The prefix is what makes `--info AC` safe. sites.parquet's own ac/an/n_called
// are recomputed over the samples in the store rather than copied from the
// source, which is the only reason they survive a subset; a capture that landed
// in the same column would undo that silently.
func TestCapturedColumnsCannotCollideWithTheStoresOwn(t *testing.T) {
	if err := ValidateInfo([]InfoField{{Name: "AC", Type: InfoInteger, Number: "A"}}); err != nil {
		t.Fatalf("AC should be capturable as info_ac: %v", err)
	}
	if got := InfoColumn("AC"); got != "info_ac" {
		t.Fatalf("InfoColumn(AC) = %q, want info_ac", got)
	}
	// And a caller naming a column by hand cannot reach one of the store's.
	err := ValidateInfo([]InfoField{{Name: "AC", Column: "ac", Type: InfoInteger, Number: "A"}})
	if err == nil {
		t.Fatal("accepted a capture into the store's own computed ac column")
	}
}

func TestDuplicateCapturedColumnsRefused(t *testing.T) {
	err := ValidateInfo([]InfoField{
		{Name: "R2", Type: InfoFloat, Number: "1"},
		{Name: "r2", Type: InfoFloat, Number: "1"},
	})
	if err == nil {
		t.Fatal("accepted two fields wanting the same column")
	}
}

// Offering values to a store that captures nothing is an error rather than a
// silent drop: the store would otherwise look like it holds them.
func TestWritingInfoToAStoreThatCapturesNoneFails(t *testing.T) {
	dir := filepath.Join(t.TempDir(), "store")
	w := infoWriter(t, dir, nil)
	err := w.WriteSiteInfo(Site{Chrom: "chr1", Pos: 1, Ref: "A", Alt: "T"},
		map[string]any{"R2": 0.5})
	if err == nil {
		t.Fatal("accepted INFO values into a store with no captured columns")
	}
	_ = w.Close()
}
