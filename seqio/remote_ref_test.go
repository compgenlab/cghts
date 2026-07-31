package seqio

import (
	"context"
	"fmt"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/iosource"
	s3client "github.com/compgenlab/cghts/iosource/s3"
)

// writeFastaWithFai writes a small FASTA and its .fai index.
func writeFastaWithFai(t *testing.T, dir string) string {
	t.Helper()
	const lineLen = 60
	seqs := map[string]string{
		"chr1":  strings.Repeat("ACGT", 45), // 180 bases
		"chr17": strings.Repeat("GGCC", 30), // 120 bases
	}
	order := []string{"chr1", "chr17"}

	path := filepath.Join(dir, "ref.fa")
	var body, fai strings.Builder
	offset := 0
	for _, name := range order {
		hdr := ">" + name + "\n"
		body.WriteString(hdr)
		offset += len(hdr)
		seq := seqs[name]
		fai.WriteString(fmt.Sprintf("%s\t%d\t%d\t%d\t%d\n", name, len(seq), offset, lineLen, lineLen+1))
		for i := 0; i < len(seq); i += lineLen {
			end := i + lineLen
			if end > len(seq) {
				end = len(seq)
			}
			line := seq[i:end] + "\n"
			body.WriteString(line)
			offset += len(line)
		}
	}
	if err := os.WriteFile(path, []byte(body.String()), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(path+".fai", []byte(fai.String()), 0o644); err != nil {
		t.Fatal(err)
	}
	return path
}

func mustRange(t *testing.T, r ReferenceReader, name string, start, end int) string {
	t.Helper()
	b, err := r.GetSequenceRange(name, start, end)
	if err != nil {
		t.Fatal(err)
	}
	return string(b)
}

// OpenReference must dispatch on any registered scheme, not just http(s).
func TestOpenReferenceOverHTTPAndS3(t *testing.T) {
	dir := t.TempDir()
	path := writeFastaWithFai(t, dir)

	local, err := OpenReference(path)
	if err != nil {
		t.Fatal(err)
	}
	defer local.Close()
	want := mustRange(t, local, "chr1", 10, 40)
	if want == "" {
		t.Fatal("local reference returned nothing; the comparison proves nothing")
	}

	ts := httptest.NewServer(http.FileServer(http.Dir(dir)))
	defer ts.Close()
	remote, err := OpenReference(ts.URL + "/ref.fa")
	if err != nil {
		t.Fatalf("open over http: %v", err)
	}
	defer remote.Close()
	if got := mustRange(t, remote, "chr1", 10, 40); got != want {
		t.Errorf("http range = %q, want %q", got, want)
	}
	if got, w := mustRange(t, remote, "chr17", 0, 20), mustRange(t, local, "chr17", 0, 20); got != w {
		t.Errorf("http chr17 = %q, want %q", got, w)
	}

	bucket := os.Getenv("VARHUB_TEST_S3_BUCKET")
	if bucket == "" {
		t.Skip("set VARHUB_TEST_S3_BUCKET (and AWS_ENDPOINT_URL) for the S3 half")
	}
	ctx := context.Background()
	for _, suffix := range []string{"", ".fai"} {
		f, err := os.Open(path + suffix)
		if err != nil {
			t.Fatal(err)
		}
		err = s3client.PutForTest(ctx, bucket, "seqio-test/ref.fa"+suffix, f)
		f.Close()
		if err != nil {
			t.Fatal(err)
		}
	}
	s3ref, err := OpenReference("s3://" + bucket + "/seqio-test/ref.fa")
	if err != nil {
		t.Fatalf("open over s3: %v", err)
	}
	defer s3ref.Close()
	if got := mustRange(t, s3ref, "chr1", 10, 40); got != want {
		t.Errorf("s3 range = %q, want %q", got, want)
	}
	if names := s3ref.Names(); len(names) != 2 {
		t.Errorf("s3 names = %v, want 2 sequences", names)
	}
	_ = iosource.Schemes()
}
