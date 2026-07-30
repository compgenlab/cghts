package annotate

import (
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"testing"

	"github.com/compgenlab/cghts/htsio/tabix"
	"github.com/compgenlab/cghts/iosource"
	"github.com/compgenlab/cghts/vcf"
)

// writeAnnSource writes a small bgzipped, tabix-indexed annotation source.
func writeAnnSource(t *testing.T, path string) {
	t.Helper()
	w := tabix.NewWriter(path, tabix.NewWriterOpts().VCF().AutoIndex())
	w.WriteHeader("##fileformat=VCFv4.2")
	w.WriteHeader("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
	for _, line := range []string{
		"chr1\t100\t.\tA\tG\t.\t.\tSIG=Benign",
		"chr1\t200\t.\tC\tT\t.\t.\tSIG=Pathogenic",
		"chr17\t7676154\t.\tC\tT\t.\t.\tSIG=Pathogenic",
	} {
		if err := w.Write(line); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
}

// annotateOne runs one record through an annotator and returns the INFO value.
func annotateOne(t *testing.T, a interface {
	Annotate(*vcf.VcfRecord) error
}, chrom string, pos int, ref, alt string) string {
	t.Helper()
	rec := vcf.NewRecord(chrom, pos, ref, []string{alt})
	if err := a.Annotate(rec); err != nil {
		t.Fatal(err)
	}
	v, ok := rec.InfoValue("sig")
	if !ok {
		return ""
	}
	return v.String()
}

// The end of the chain this whole line of work is for: an annotator reading its
// source over HTTP range requests must produce exactly what it produces from
// the same file on disk.
func TestTabixAnnotatorOverRemoteSource(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "src.vcf.gz")
	writeAnnSource(t, path)

	ts := httptest.NewServer(http.FileServer(http.Dir(dir)))
	defer ts.Close()

	opts := TabixOptions{Name: "sig", Filename: path, Col: 8}

	local, err := NewTabixAnnotator(opts)
	if err != nil {
		t.Fatal(err)
	}
	defer local.Close()

	src, err := iosource.NewHTTPRange(ts.URL + "/src.vcf.gz")
	if err != nil {
		t.Fatal(err)
	}
	rs, err := iosource.ReadSeeker(src)
	if err != nil {
		t.Fatal(err)
	}
	idx, _, err := iosource.ResolveSibling(ts.URL+"/src.vcf.gz",
		[]string{".tbi", ".csi"}, iosource.HTTPSibling)
	if err != nil {
		t.Fatal(err)
	}
	defer idx.Close()
	rdr, err := tabix.NewReaderFromSource(rs, idx, tabix.WithCloser(src))
	if err != nil {
		t.Fatal(err)
	}

	remoteOpts := opts
	remoteOpts.Reader = rdr
	remote, err := NewTabixAnnotator(remoteOpts)
	if err != nil {
		t.Fatal(err)
	}
	defer remote.Close()

	for _, q := range []struct {
		chrom, ref, alt string
		pos             int
	}{
		{"chr1", "A", "G", 100},
		{"chr1", "C", "T", 200},
		{"chr17", "C", "T", 7676154},
		{"chr1", "A", "G", 999}, // no hit
	} {
		want := annotateOne(t, local, q.chrom, q.pos, q.ref, q.alt)
		got := annotateOne(t, remote, q.chrom, q.pos, q.ref, q.alt)
		if want != got {
			t.Errorf("%s:%d %s>%s: remote %q, local %q", q.chrom, q.pos, q.ref, q.alt, got, want)
		}
	}
}

// A supplied reader is owned by the annotator, matching what Close does for one
// this package opened. Leaving that to the caller would leak by omission.
func TestSuppliedReaderIsClosed(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "src.vcf.gz")
	writeAnnSource(t, path)

	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	idx, err := os.Open(path + ".tbi")
	if err != nil {
		t.Fatal(err)
	}
	defer idx.Close()

	rdr, err := tabix.NewReaderFromSource(f, idx, tabix.WithCloser(f))
	if err != nil {
		t.Fatal(err)
	}
	a, err := NewTabixAnnotator(TabixOptions{Name: "sig", Filename: path, Col: 8, Reader: rdr})
	if err != nil {
		t.Fatal(err)
	}
	if err := a.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	if _, err := f.Seek(0, 0); err == nil {
		t.Error("the annotator did not close the reader it was given")
	}
}
