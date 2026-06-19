package bed

import (
	"bytes"
	"compress/gzip"
	"errors"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestFormatScoreInt(t *testing.T) {
	cases := []struct {
		in   float64
		want string
	}{
		{1.5, "1"},
		{2.9, "2"},
		{1.9, "1"},
		{-1.9, "-1"},
		{0.0, "0"},
		{10, "10"},
		{1000000, "1000000"},
	}
	for _, c := range cases {
		if got := formatScoreInt(c.in); got != c.want {
			t.Errorf("formatScoreInt(%v) = %q, want %q", c.in, got, c.want)
		}
	}
}

func TestFormatScoreFloat(t *testing.T) {
	cases := []struct {
		in   float64
		want string
	}{
		{1.5, "1.5"},
		{1.0, "1"},
		{0.0, "0"},
		{1234567.89, "1234567.89"},
		{-1.9, "-1.9"},
		{10, "10"},
	}
	for _, c := range cases {
		if got := formatScoreFloat(c.in); got != c.want {
			t.Errorf("formatScoreFloat(%v) = %q, want %q", c.in, got, c.want)
		}
	}
}

func readAll(t *testing.T, r *BedReader) []*BedRecord {
	t.Helper()
	var recs []*BedRecord
	for {
		rec, err := r.NextRecord()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("NextRecord: %v", err)
		}
		recs = append(recs, rec)
	}
	return recs
}

func TestReaderEdgeCases(t *testing.T) {
	input := "#comment\n\nchr1\t100\t200\nchr2\t50\t75\tregionB\t2.9\t-\tfoo\tbar\nbad\tline\n"
	r, err := NewBedReader(strings.NewReader(input))
	if err != nil {
		t.Fatal(err)
	}
	recs := readAll(t, r)

	if len(recs) != 2 {
		t.Fatalf("got %d records, want 2 (comment/blank/<3-col lines skipped)", len(recs))
	}

	// BED3 line: name defaults to "" but HasName is true (file record).
	r0 := recs[0]
	if r0.Ref != "chr1" || r0.Start != 100 || r0.End != 200 {
		t.Errorf("rec0 coords = %s:%d-%d", r0.Ref, r0.Start, r0.End)
	}
	if !r0.HasName || r0.Name != "" || r0.Strand != StrandNone || r0.Score != 0 {
		t.Errorf("rec0 = %+v, want HasName=true Name=\"\" Strand=none Score=0", r0)
	}
	if r0.Extras != nil {
		t.Errorf("rec0 extras = %v, want nil", r0.Extras)
	}

	// Full line with extras.
	r1 := recs[1]
	if r1.Name != "regionB" || r1.Score != 2.9 || r1.Strand != StrandMinus {
		t.Errorf("rec1 = %+v", r1)
	}
	if len(r1.Extras) != 2 || r1.Extras[0] != "foo" || r1.Extras[1] != "bar" {
		t.Errorf("rec1 extras = %v, want [foo bar]", r1.Extras)
	}
}

func TestReaderParseErrorAborts(t *testing.T) {
	r, err := NewBedReader(strings.NewReader("chr1\tnotanumber\t200\n"))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := r.NextRecord(); err == nil || err == io.EOF {
		t.Fatalf("expected a parse error, got %v", err)
	}
}

func TestReaderGzipFile(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "in.bed.gz")

	w, err := OpenBedWriter(path)
	if err != nil {
		t.Fatal(err)
	}
	if err := w.WriteRecord(NewBed6("chr1", 100, 200, "a", 1.0, StrandPlus)); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	// Confirm it really is gzip on disk.
	raw, err := io.ReadAll(mustOpen(t, path))
	if err != nil {
		t.Fatal(err)
	}
	if len(raw) < 2 || raw[0] != 0x1f || raw[1] != 0x8b {
		t.Fatalf("expected gzip magic bytes, got % x", raw[:min(2, len(raw))])
	}

	r, err := NewBedFile(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	recs := readAll(t, r)
	if len(recs) != 1 || recs[0].Name != "a" {
		t.Fatalf("gzip round-trip failed: %+v", recs)
	}
}

func writeRecords(t *testing.T, recs []*BedRecord, opts *BedWriterOpts) string {
	t.Helper()
	var buf bytes.Buffer
	w := NewBedWriter(&buf, opts)
	for _, rec := range recs {
		if err := w.WriteRecord(rec); err != nil {
			t.Fatalf("WriteRecord: %v", err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	return buf.String()
}

func TestWriterColumnModes(t *testing.T) {
	// A record read from a 6-col + extras line, plus a bare BED3 line.
	rec, ok, err := parseBedLine("chr1\t100\t200\tregionA\t1.5\t+\tfoo\tbar", false)
	if err != nil || !ok {
		t.Fatalf("parseBedLine: ok=%v err=%v", ok, err)
	}

	// ColumnsAuto, float score: full passthrough including extras.
	if got := writeRecords(t, []*BedRecord{rec}, NewBedWriterOpts()); got != "chr1\t100\t200\tregionA\t1.5\t+\tfoo\tbar\n" {
		t.Errorf("auto/float = %q", got)
	}
	// ColumnsAuto, forceScoreInt: 1.5 -> 1.
	if got := writeRecords(t, []*BedRecord{rec}, NewBedWriterOpts().ForceScoreInt(true)); got != "chr1\t100\t200\tregionA\t1\t+\tfoo\tbar\n" {
		t.Errorf("auto/int = %q", got)
	}
	// Columns3: drop everything past end.
	if got := writeRecords(t, []*BedRecord{rec}, NewBedWriterOpts().Columns(Columns3)); got != "chr1\t100\t200\n" {
		t.Errorf("cols3 = %q", got)
	}
	// Columns6: drop extras, keep BED6.
	if got := writeRecords(t, []*BedRecord{rec}, NewBedWriterOpts().Columns(Columns6)); got != "chr1\t100\t200\tregionA\t1.5\t+\n" {
		t.Errorf("cols6 = %q", got)
	}

	// Bare BED3 record: ColumnsAuto emits only 3 cols (HasName false).
	bed3 := NewBed3("chr3", 5, 9)
	if got := writeRecords(t, []*BedRecord{bed3}, NewBedWriterOpts()); got != "chr3\t5\t9\n" {
		t.Errorf("bed3 auto = %q", got)
	}
	// Columns6 forces name/score/strand even on a BED3 record (none -> +).
	if got := writeRecords(t, []*BedRecord{bed3}, NewBedWriterOpts().Columns(Columns6)); got != "chr3\t5\t9\t\t0\t+\n" {
		t.Errorf("bed3 cols6 = %q", got)
	}
}

func TestWriterForcesStrand(t *testing.T) {
	// A record with no strand column should write "+" in the 6-col path.
	rec, _, _ := parseBedLine("chr1\t100\t200\tx\t0", false)
	if got := writeRecords(t, []*BedRecord{rec}, NewBedWriterOpts().ForceScoreInt(true)); got != "chr1\t100\t200\tx\t0\t+\n" {
		t.Errorf("got %q, want strand forced to +", got)
	}
}

func TestExtend(t *testing.T) {
	plus := NewBed6("chr1", 100, 200, "p", 0, StrandPlus)
	if e := plus.Extend5(50); e.Start != 50 || e.End != 200 {
		t.Errorf("plus Extend5: %d-%d", e.Start, e.End)
	}
	if e := plus.Extend3(50); e.Start != 100 || e.End != 250 {
		t.Errorf("plus Extend3: %d-%d", e.Start, e.End)
	}

	minus := NewBed6("chr1", 100, 200, "m", 0, StrandMinus)
	if e := minus.Extend5(50); e.Start != 100 || e.End != 250 {
		t.Errorf("minus Extend5: %d-%d", e.Start, e.End)
	}
	if e := minus.Extend3(50); e.Start != 50 || e.End != 200 {
		t.Errorf("minus Extend3: %d-%d", e.Start, e.End)
	}

	// Clamp at 0.
	near := NewBed6("chr1", 10, 20, "n", 0, StrandPlus)
	if e := near.Extend5(100); e.Start != 0 {
		t.Errorf("clamp: start = %d, want 0", e.Start)
	}

	// Extend returns a copy; original is unchanged.
	if plus.Start != 100 || plus.End != 200 {
		t.Errorf("original mutated: %d-%d", plus.Start, plus.End)
	}
}

func TestTabixRoundTrip(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "out.bed.gz")

	recs := []*BedRecord{
		NewBed6("chr1", 100, 200, "a", 1, StrandPlus),
		NewBed6("chr1", 500, 600, "b", 2, StrandMinus),
		NewBed6("chr1", 1000, 1100, "c", 3, StrandPlus),
		NewBed6("chr2", 50, 75, "d", 4, StrandPlus),
	}

	w, err := OpenBedWriter(path, NewBedWriterOpts().ForceScoreInt(true).Tabix(true))
	if err != nil {
		t.Fatal(err)
	}
	for _, rec := range recs {
		if err := w.WriteRecord(rec); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	idx, err := NewIndexedBedReader(path)
	if err != nil {
		t.Fatalf("NewIndexedBedReader: %v", err)
	}
	defer idx.Close()

	// Query a region overlapping only the second chr1 record.
	seq, err := idx.Query("chr1", 450, 650)
	if err != nil {
		t.Fatal(err)
	}
	var got []*BedRecord
	for rec, err := range seq {
		if err != nil {
			t.Fatalf("query iter: %v", err)
		}
		got = append(got, rec)
	}
	if len(got) != 1 || got[0].Name != "b" || got[0].Start != 500 {
		t.Fatalf("query [450,650) = %+v, want single record b@500", got)
	}

	// Query covering all of chr1.
	seq, _ = idx.Query("chr1", 0, 100000)
	n := 0
	for _, err := range seq {
		if err != nil {
			t.Fatal(err)
		}
		n++
	}
	if n != 3 {
		t.Fatalf("chr1 query returned %d records, want 3", n)
	}
}

// mustOpen opens a file for reading or fails the test.
func mustOpen(t *testing.T, path string) io.Reader {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { f.Close() })
	return f
}

// sanity check that a gzip stream we wrote decompresses.
func TestGzipStreamValid(t *testing.T) {
	var buf bytes.Buffer
	gz := gzip.NewWriter(&buf)
	io.WriteString(gz, "chr1\t1\t2\n")
	gz.Close()

	gr, err := gzip.NewReader(&buf)
	if err != nil {
		t.Fatal(err)
	}
	out, err := io.ReadAll(gr)
	if err != nil && !errors.Is(err, io.EOF) {
		t.Fatal(err)
	}
	if string(out) != "chr1\t1\t2\n" {
		t.Fatalf("got %q", out)
	}
}
