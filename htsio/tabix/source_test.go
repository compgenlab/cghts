package tabix

import (
	"fmt"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"sync/atomic"
	"testing"

	"github.com/compgenlab/cghts/iosource"
)

// countingFS serves a directory over HTTP with Range support and records how
// many body bytes it handed out, so a test can prove an indexed query fetched
// only what it needed rather than quietly pulling the whole file.
type countingFS struct {
	dir   string
	bytes atomic.Int64
	reqs  atomic.Int64
}

type countingWriter struct {
	http.ResponseWriter
	n *atomic.Int64
}

func (w countingWriter) Write(p []byte) (int, error) {
	n, err := w.ResponseWriter.Write(p)
	w.n.Add(int64(n))
	return n, err
}

func (c *countingFS) ServeHTTP(w http.ResponseWriter, r *http.Request) {
	c.reqs.Add(1)
	http.FileServer(http.Dir(c.dir)).ServeHTTP(countingWriter{w, &c.bytes}, r)
}

// writeIndexedFixture writes a bgzipped, tabix-indexed VCF with enough records
// to make "fetched only part of it" a meaningful claim.
func writeIndexedFixture(t *testing.T, path string) { writeIndexedFixtureN(t, path, 20000) }

func writeIndexedFixtureN(t *testing.T, path string, perChrom int) {
	t.Helper()
	w := NewWriter(path, NewWriterOpts().VCF().AutoIndex())
	w.WriteHeader("##fileformat=VCFv4.2")
	w.WriteHeader("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
	for _, chrom := range []string{"chr1", "chr2", "chr17"} {
		for pos := 1; pos <= perChrom; pos++ {
			line := fmt.Sprintf("%s\t%d\t.\tA\tG\t.\t.\tDP=%d;PAD=%s",
				chrom, pos*10, pos, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
			if err := w.Write(line); err != nil {
				t.Fatal(err)
			}
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
}

// collect runs a query and returns its records, or the error it produced. Both
// are returned so a remote reader can be compared against a local one on
// failures as well as successes — an absent reference is an error, and the two
// must agree about that too.
func collect(t *testing.T, r *Reader, ref string, start, end int) ([]string, error) {
	t.Helper()
	seq, err := r.Query(ref, start, end)
	if err != nil {
		return nil, err
	}
	var out []string
	for rec, err := range seq {
		if err != nil {
			return nil, err
		}
		out = append(out, rec.Line)
	}
	return out, nil
}

// The acceptance test for reading a tabix file without downloading it: identical
// records to the same file on disk, and only a small part of the file fetched.
func TestQueryOverHTTPRangeMatchesLocal(t *testing.T) {
	dir := t.TempDir()
	local := filepath.Join(dir, "test.vcf.gz")
	writeIndexedFixture(t, local)

	fi, err := os.Stat(local)
	if err != nil {
		t.Fatal(err)
	}
	idxSize := int64(0)
	if s, err := os.Stat(local + ".tbi"); err == nil {
		idxSize = s.Size()
	}

	cfs := &countingFS{dir: dir}
	ts := httptest.NewServer(cfs)
	defer ts.Close()

	lr, err := NewReader(local)
	if err != nil {
		t.Fatal(err)
	}
	defer lr.Close()

	src, err := iosource.NewHTTPRange(ts.URL + "/test.vcf.gz")
	if err != nil {
		t.Fatal(err)
	}
	rs, err := iosource.ReadSeeker(src)
	if err != nil {
		t.Fatal(err)
	}
	idxRC, suffix, err := iosource.ResolveSibling(ts.URL+"/test.vcf.gz",
		[]string{".tbi", ".csi"}, iosource.HTTPSibling)
	if err != nil {
		t.Fatal(err)
	}
	defer idxRC.Close()
	if suffix != ".tbi" {
		t.Fatalf("resolved sibling %q, want .tbi", suffix)
	}

	rr, err := NewReaderFromSource(rs, idxRC, WithCloser(src))
	if err != nil {
		t.Fatal(err)
	}
	defer rr.Close()

	afterIndex := cfs.bytes.Load()

	for _, q := range []struct {
		ref         string
		start, end  int
		wantRecords bool
	}{
		{"chr1", 1000, 2000, true},
		{"chr17", 50000, 51000, true},
		{"chr2", 199990, 200010, true},
		{"chr1", 100000000, 100001000, false}, // past the end
		{"chrX", 1, 1000, false},              // absent reference: an error, and both must say so
	} {
		name := fmt.Sprintf("%s:%d-%d", q.ref, q.start, q.end)
		want, wantErr := collect(t, lr, q.ref, q.start, q.end)
		got, gotErr := collect(t, rr, q.ref, q.start, q.end)

		switch {
		case wantErr == nil && gotErr != nil:
			t.Fatalf("%s: remote failed where local succeeded: %v", name, gotErr)
		case wantErr != nil && gotErr == nil:
			t.Fatalf("%s: remote succeeded where local failed with %v", name, wantErr)
		case wantErr != nil && gotErr.Error() != wantErr.Error():
			t.Fatalf("%s: errors differ\n remote: %v\n  local: %v", name, gotErr, wantErr)
		case wantErr != nil:
			continue // both failed the same way
		}

		if len(want) != len(got) {
			t.Fatalf("%s: got %d records, local gave %d", name, len(got), len(want))
		}
		for i := range want {
			if want[i] != got[i] {
				t.Fatalf("%s: record %d differs\n remote: %s\n  local: %s", name, i, got[i], want[i])
			}
		}
		if q.wantRecords && len(want) == 0 {
			t.Fatalf("%s: fixture produced no records; the test proves nothing", name)
		}
	}

	// Sanity only: these queries span three chromosomes of a small fixture, so
	// touching much of it is legitimate. The real economy claim is measured on a
	// large file with one narrow query, in TestQueryOverHTTPRangeFetchesOnlyWhatItNeeds.
	dataFetched := cfs.bytes.Load() - afterIndex
	if dataFetched >= fi.Size()*2 {
		t.Errorf("queries fetched %d bytes of a %d-byte file — that is not range reading at all",
			dataFetched, fi.Size())
	}
	t.Logf("index %d B, data fetched %d B of %d B total", idxSize, dataFetched, fi.Size())
}

// The economy claim, stated where it can actually be violated: one narrow query
// against a file far larger than the region asked for must move a small fraction
// of it. A silent fallback to a full-object GET would still return correct
// records, so correctness alone cannot catch it.
func TestQueryOverHTTPRangeFetchesOnlyWhatItNeeds(t *testing.T) {
	dir := t.TempDir()
	local := filepath.Join(dir, "big.vcf.gz")
	writeIndexedFixtureN(t, local, 300000)

	fi, err := os.Stat(local)
	if err != nil {
		t.Fatal(err)
	}

	cfs := &countingFS{dir: dir}
	ts := httptest.NewServer(cfs)
	defer ts.Close()

	src, err := iosource.NewHTTPRange(ts.URL + "/big.vcf.gz")
	if err != nil {
		t.Fatal(err)
	}
	rs, err := iosource.ReadSeeker(src)
	if err != nil {
		t.Fatal(err)
	}
	idxRC, _, err := iosource.ResolveSibling(ts.URL+"/big.vcf.gz",
		[]string{".tbi", ".csi"}, iosource.HTTPSibling)
	if err != nil {
		t.Fatal(err)
	}
	defer idxRC.Close()

	rr, err := NewReaderFromSource(rs, idxRC, WithCloser(src))
	if err != nil {
		t.Fatal(err)
	}
	defer rr.Close()

	baseline := cfs.bytes.Load()
	got, err := collect(t, rr, "chr2", 1000000, 1001000)
	if err != nil {
		t.Fatal(err)
	}
	if len(got) == 0 {
		t.Fatal("query returned no records; the measurement would be meaningless")
	}
	fetched := cfs.bytes.Load() - baseline

	// A generous ceiling: the point is to catch a full-file read, not to pin the
	// exact block count, which depends on BGZF block packing.
	if limit := fi.Size() / 10; fetched > limit {
		t.Errorf("one narrow query fetched %d bytes of a %d-byte file (>%d); "+
			"ranges are far wider than the region, or it fell back to a full GET",
			fetched, fi.Size(), limit)
	}
	t.Logf("%d records; fetched %d B of %d B (%.2f%%)",
		len(got), fetched, fi.Size(), 100*float64(fetched)/float64(fi.Size()))
}

// A CSI-indexed file must work through the same entry point, since the caller
// resolves the sidecar by trying both suffixes and does not tell us which won.
func TestNewReaderFromSourceDetectsIndexFormat(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "t.vcf.gz")
	writeIndexedFixture(t, path)

	data, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer data.Close()
	idx, err := os.Open(path + ".tbi")
	if err != nil {
		t.Fatal(err)
	}
	defer idx.Close()

	r, err := NewReaderFromSource(data, idx)
	if err != nil {
		t.Fatalf("NewReaderFromSource: %v", err)
	}
	defer r.Close()
	recs, err := collect(t, r, "chr1", 1000, 2000)
	if err != nil {
		t.Fatal(err)
	}
	if len(recs) == 0 {
		t.Error("no records from a TBI-indexed source")
	}
}

// Something that is not an index at all must be rejected clearly, rather than
// being misparsed into an index that silently returns nothing.
func TestNewReaderFromSourceRejectsNonIndex(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "t.vcf.gz")
	writeIndexedFixture(t, path)

	data, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer data.Close()
	// The data file is valid BGZF but is not an index.
	notIdx, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer notIdx.Close()

	if _, err := NewReaderFromSource(data, notIdx); err == nil {
		t.Fatal("a non-index stream was accepted as an index")
	}
}

// Close must release the source only when the caller handed it over.
func TestCloserOwnership(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "t.vcf.gz")
	writeIndexedFixture(t, path)

	open := func(t *testing.T) (*os.File, *os.File) {
		t.Helper()
		d, err := os.Open(path)
		if err != nil {
			t.Fatal(err)
		}
		i, err := os.Open(path + ".tbi")
		if err != nil {
			t.Fatal(err)
		}
		return d, i
	}

	// Without WithCloser the caller still owns the file: it stays usable.
	d, i := open(t)
	r, err := NewReaderFromSource(d, i)
	if err != nil {
		t.Fatal(err)
	}
	if err := r.Close(); err != nil {
		t.Fatalf("Close without an owned resource: %v", err)
	}
	if _, err := d.Seek(0, 0); err != nil {
		t.Errorf("reader closed a source it did not own: %v", err)
	}
	d.Close()
	i.Close()

	// With WithCloser it is released.
	d, i = open(t)
	r, err = NewReaderFromSource(d, i, WithCloser(d))
	if err != nil {
		t.Fatal(err)
	}
	if err := r.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	if _, err := d.Seek(0, 0); err == nil {
		t.Error("reader did not close a source it was given ownership of")
	}
	i.Close()
}
