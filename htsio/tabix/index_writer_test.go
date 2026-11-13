package tabix

import (
	"errors"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"testing"

	"github.com/compgenlab/cghts/htsio/bgzf"
)

// writeBGZF writes the given lines to a BGZF file with no index.
func writeBGZF(t *testing.T, path string, lines ...string) {
	t.Helper()
	w, err := bgzf.NewBGZipFile(path)
	if err != nil {
		t.Fatal(err)
	}
	for _, l := range lines {
		if _, err := io.WriteString(w, l+"\n"); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
}

func TestIndexWriterBED(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "data.bed.gz")
	writeBGZF(t, path,
		"chr1\t90\t110\tgeneA",
		"chr1\t145\t155\tenhB",
		"chr2\t400\t600\tgeneC") // single chr2 record (and the last line)

	if err := NewIndexWriter(NewWriterOpts().BED()).WriteIndex(path); err != nil {
		t.Fatalf("WriteIndex: %v", err)
	}

	r, err := NewReader(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()

	cases := []struct {
		ref        string
		start, end int
		want       string
	}{
		{"chr1", 99, 100, "chr1\t90\t110\tgeneA"},
		{"chr1", 149, 150, "chr1\t145\t155\tenhB"},
		{"chr2", 499, 500, "chr2\t400\t600\tgeneC"},
	}
	for _, c := range cases {
		got := queryLines(t, r, c.ref, c.start, c.end)
		if len(got) != 1 || got[0] != c.want {
			t.Errorf("query %s:%d-%d = %v, want [%q]", c.ref, c.start, c.end, got, c.want)
		}
	}
}

func TestIndexWriterVCFWithHeader(t *testing.T) {
	// A '#'-prefixed header line must be skipped (meta='#'); 1-based positions.
	dir := t.TempDir()
	path := filepath.Join(dir, "v.vcf.gz")
	writeBGZF(t, path,
		"##fileformat=VCFv4.2",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
		"chr1\t100\t.\tA\tG\t.\tPASS\t.",
		"chr2\t500\t.\tC\tT\t.\tPASS\t.")

	if err := NewIndexWriter(NewWriterOpts().VCF()).WriteIndex(path); err != nil {
		t.Fatalf("WriteIndex: %v", err)
	}
	r, err := NewReader(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	if got := queryLines(t, r, "chr2", 499, 500); len(got) != 1 {
		t.Errorf("chr2 query = %v, want 1 row", got)
	}
	if got := queryLines(t, r, "chr1", 99, 100); len(got) != 1 {
		t.Errorf("chr1 query = %v, want 1 row", got)
	}
}

// TestIndexWriterUnsorted covers the order check: an unsorted file yields an
// UnsortedError and no index, rather than an index whose linear offsets lie.
func TestIndexWriterUnsorted(t *testing.T) {
	header := []string{"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"}
	rec := func(ref string, pos int) string {
		return ref + "\t" + strconv.Itoa(pos) + "\t.\tA\tG\t.\tPASS\t."
	}
	cases := []struct {
		name    string
		lines   []string
		want    UnsortedError
		wantMsg string
	}{
		{
			name:    "position goes backwards",
			lines:   []string{rec("chr1", 300), rec("chr1", 100)},
			want:    UnsortedError{Line: 3, Ref: "chr1", Pos: 100, PrevRef: "chr1", PrevPos: 300},
			wantMsg: "line 3 has chr1:100, which follows chr1:300",
		},
		{
			name:    "reference resumes",
			lines:   []string{rec("chr1", 100), rec("chr2", 500), rec("chr1", 300)},
			want:    UnsortedError{Line: 4, Ref: "chr1", Pos: 300, PrevRef: "chr2", PrevPos: 500, Revisit: true},
			wantMsg: "line 4 has chr1:300, but chr1 already gave way to chr2",
		},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			path := filepath.Join(t.TempDir(), "v.vcf.gz")
			writeBGZF(t, path, append(header, c.lines...)...)

			err := NewIndexWriter(NewWriterOpts().VCF()).WriteIndex(path)
			var ue *UnsortedError
			if !errors.As(err, &ue) {
				t.Fatalf("WriteIndex error = %v, want *UnsortedError", err)
			}
			if *ue != c.want {
				t.Errorf("UnsortedError = %+v, want %+v", *ue, c.want)
			}
			if ue.Error() != "not coordinate sorted: "+c.wantMsg {
				t.Errorf("message = %q, want it to end with %q", ue.Error(), c.wantMsg)
			}
			if _, err := os.Stat(path + ".tbi"); err == nil {
				t.Errorf("an index was written for an unsorted file")
			}
		})
	}
}

// TestIndexWriterEqualPositions covers records sharing a position (multiallelic
// sites split across lines), which is sorted, not a violation.
func TestIndexWriterEqualPositions(t *testing.T) {
	path := filepath.Join(t.TempDir(), "v.vcf.gz")
	writeBGZF(t, path,
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
		"chr1\t100\t.\tA\tG\t.\tPASS\t.",
		"chr1\t100\t.\tA\tT\t.\tPASS\t.")

	if err := NewIndexWriter(NewWriterOpts().VCF()).WriteIndex(path); err != nil {
		t.Fatalf("WriteIndex: %v", err)
	}
	r, err := NewReader(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	if got := queryLines(t, r, "chr1", 99, 100); len(got) != 2 {
		t.Errorf("chr1:100 query = %v, want both records", got)
	}
}

// TestIndexWriterZeroBasedUnsortedPos checks that a BED violation is reported in
// BED's own 0-based coordinates.
func TestIndexWriterZeroBasedUnsortedPos(t *testing.T) {
	path := filepath.Join(t.TempDir(), "data.bed.gz")
	writeBGZF(t, path, "chr1\t300\t400\tb", "chr1\t100\t200\ta")

	err := NewIndexWriter(NewWriterOpts().BED()).WriteIndex(path)
	var ue *UnsortedError
	if !errors.As(err, &ue) {
		t.Fatalf("WriteIndex error = %v, want *UnsortedError", err)
	}
	if ue.Pos != 100 || ue.PrevPos != 300 || ue.Line != 2 {
		t.Errorf("UnsortedError = %+v, want line 2, 100 after 300", *ue)
	}
}

// A GTF compressed as a header stream and a body stream, then concatenated,
// carries an EOF marker in mid-file. htslib reads and indexes straight through
// it; the index writer panicked on it (bgzf.Reader.ReadByte indexed an empty
// buffer), and a query whose chunk crossed the join lost its records.
func TestIndexWriterConcatenatedStreams(t *testing.T) {
	dir := t.TempDir()
	hdr := filepath.Join(dir, "hdr.gz")
	body := filepath.Join(dir, "body.gz")
	writeBGZF(t, hdr, "##description: test", "##format: gtf")
	var lines []string
	for i := 0; i < 3000; i++ {
		lines = append(lines, "chr1\tHAVANA\texon\t"+strconv.Itoa(1000+i*100)+"\t"+strconv.Itoa(1050+i*100)+
			"\t.\t+\t.\tgene_id \"G"+strconv.Itoa(i)+"\"; padding \"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\";")
	}
	writeBGZF(t, body, lines...)
	a, _ := os.ReadFile(hdr)
	b, _ := os.ReadFile(body)
	path := filepath.Join(dir, "concat.gtf.gz")
	if err := os.WriteFile(path, append(a, b...), 0o644); err != nil {
		t.Fatal(err)
	}

	if err := NewIndexWriter(NewWriterOpts().GFF().Meta('#')).WriteIndex(path); err != nil {
		t.Fatalf("WriteIndex: %v", err)
	}
	r, err := NewReader(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()

	// The first record, the last, and the whole contig: none may be lost.
	if got := queryLines(t, r, "chr1", 1000, 1001); len(got) != 1 || got[0] != lines[0] {
		t.Errorf("first record: got %v", got)
	}
	last := 1000 + 2999*100
	if got := queryLines(t, r, "chr1", last, last+1); len(got) != 1 || got[0] != lines[2999] {
		t.Errorf("last record: got %v", got)
	}
	if got := queryLines(t, r, "chr1", 0, 1<<29); len(got) != len(lines) {
		t.Errorf("whole contig: %d records, want %d", len(got), len(lines))
	}
}
