package bam

import (
	"context"
	"fmt"
	"net/http"
	"net/http/httptest"
	"os"
	"os/exec"
	"path/filepath"
	"testing"

	"github.com/compgenlab/cghts/htsio"
	"github.com/compgenlab/cghts/iosource"
)

// writeSortedBAM writes a coordinate-sorted BAM with its .bai, so region
// queries have something to index.
func writeSortedBAM(t *testing.T, dir string) string {
	t.Helper()
	path := filepath.Join(dir, "test.bam")
	h := htsio.NewSamHeader()
	h.AddLine("@HD\tVN:1.6\tSO:coordinate")
	h.AddLine("@SQ\tSN:chr1\tLN:10000")
	h.AddLine("@SQ\tSN:chr2\tLN:10000")
	w, err := NewSortedWriter(path, h, true)
	if err != nil {
		t.Fatal(err)
	}
	for _, ref := range []string{"chr1", "chr2"} {
		for pos := 1; pos <= 500; pos++ {
			rec := &htsio.SamRecord{
				ReadName: fmt.Sprintf("%s_r%d", ref, pos),
				Flag:     0, RefName: ref, Pos: pos, MapQ: 60,
				Cigar: "10M", RefNext: "*", Seq: "ACGTACGTAC", Qual: "**********",
			}
			if err := w.Write(rec); err != nil {
				t.Fatal(err)
			}
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	// cghts has no BAI writer, so the index comes from samtools. Without it
	// there is nothing to query and the test skips rather than pretending.
	if _, err := exec.LookPath("samtools"); err != nil {
		t.Skip("samtools not available to build the .bai fixture")
	}
	if out, err := exec.Command("samtools", "index", path).CombinedOutput(); err != nil {
		t.Skipf("samtools index failed: %v: %s", err, out)
	}
	return path
}

func queryNames(t *testing.T, r *Reader, ref string, start, end int) []string {
	t.Helper()
	seq, err := r.Query(ref, start, end)
	if err != nil {
		t.Fatal(err)
	}
	var out []string
	for rec, err := range seq {
		if err != nil {
			t.Fatal(err)
		}
		out = append(out, rec.ReadName)
	}
	return out
}

// An indexed region query over a BAM served by HTTP must match the same file on
// disk — the .bai resolved through the same transport as the data.
func TestRemoteBamIndexedQuery(t *testing.T) {
	dir := t.TempDir()
	path := writeSortedBAM(t, dir)
	if _, err := os.Stat(path + ".bai"); err != nil {
		t.Skipf("no .bai produced by the sorted writer: %v", err)
	}

	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	local, err := NewReader(f, path, nil)
	if err != nil {
		t.Fatal(err)
	}
	defer local.Close()

	ts := httptest.NewServer(http.FileServer(http.Dir(dir)))
	defer ts.Close()

	src, err := iosource.Open(context.Background(), ts.URL+"/test.bam")
	if err != nil {
		t.Fatal(err)
	}
	remote, err := NewReaderFromSource(src, ts.URL+"/test.bam", nil)
	if err != nil {
		t.Fatalf("NewReaderFromSource: %v", err)
	}
	defer remote.Close()

	for _, q := range []struct {
		ref        string
		start, end int
	}{
		{"chr1", 100, 200},
		{"chr2", 400, 450},
		{"chr1", 1, 10},
	} {
		want := queryNames(t, local, q.ref, q.start, q.end)
		got := queryNames(t, remote, q.ref, q.start, q.end)
		if len(want) == 0 {
			t.Fatalf("%s:%d-%d: local returned nothing; the comparison proves nothing", q.ref, q.start, q.end)
		}
		if len(want) != len(got) {
			t.Fatalf("%s:%d-%d: remote %d records, local %d", q.ref, q.start, q.end, len(got), len(want))
		}
		for i := range want {
			if want[i] != got[i] {
				t.Fatalf("%s:%d-%d record %d: remote %q, local %q", q.ref, q.start, q.end, i, got[i], want[i])
			}
		}
	}
}

// The dispatcher must route a remote locator to the source-backed constructor,
// so callers get indexed access without naming the format themselves.
func TestOpenSamReaderRemote(t *testing.T) {
	dir := t.TempDir()
	writeSortedBAM(t, dir)
	ts := httptest.NewServer(http.FileServer(http.Dir(dir)))
	defer ts.Close()

	r, err := htsio.OpenSamReader(context.Background(), ts.URL+"/test.bam")
	if err != nil {
		t.Fatalf("OpenSamReader: %v", err)
	}
	defer r.Close()

	br, ok := r.(*Reader)
	if !ok {
		t.Fatalf("got %T, want *bam.Reader", r)
	}
	if names := queryNames(t, br, "chr1", 100, 200); len(names) == 0 {
		t.Error("indexed query through OpenSamReader returned nothing")
	}
}
