package cram

import (
	"context"
	"fmt"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"testing"

	"github.com/compgenlab/cghts/htsio"
	"github.com/compgenlab/cghts/iosource"
)

// serveTestdata copies the named fixtures into a temp dir and serves it.
func serveTestdata(t *testing.T, names ...string) (*httptest.Server, string) {
	t.Helper()
	dir := t.TempDir()
	for _, n := range names {
		b, err := os.ReadFile(filepath.Join("testdata", n))
		if err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(filepath.Join(dir, n), b, 0o644); err != nil {
			t.Fatal(err)
		}
	}
	ts := httptest.NewServer(http.FileServer(http.Dir(dir)))
	t.Cleanup(ts.Close)
	return ts, dir
}

func readAll(t *testing.T, r *Reader) []string {
	t.Helper()
	var out []string
	for rec, err := range r.Records() {
		if err != nil {
			t.Fatalf("Records: %v", err)
		}
		out = append(out, fmt.Sprintf("%s\t%d\t%s\t%d\t%s", rec.ReadName, rec.Flag, rec.RefName, rec.Pos, rec.Seq))
	}
	return out
}

// A remote CRAM with a LOCAL reference. The reference is independent of where
// the data lives, so this combination has to work as readily as the others.
func TestRemoteCramLocalReference(t *testing.T) {
	ts, _ := serveTestdata(t, "test_raw.cram", "test_raw.cram.crai")
	localRef := filepath.Join("testdata", "ref.fa")

	local, err := NewReader(filepath.Join("testdata", "test_raw.cram"), localRef)
	if err != nil {
		t.Fatal(err)
	}
	defer local.Close()
	want := readAll(t, local)
	if len(want) == 0 {
		t.Fatal("local CRAM yielded no records; the comparison proves nothing")
	}

	src, err := iosource.Open(context.Background(), ts.URL+"/test_raw.cram")
	if err != nil {
		t.Fatal(err)
	}
	remote, err := NewReaderFromSource(src, ts.URL+"/test_raw.cram", localRef, nil)
	if err != nil {
		t.Fatalf("remote CRAM with local reference: %v", err)
	}
	defer remote.Close()

	got := readAll(t, remote)
	if len(got) != len(want) {
		t.Fatalf("remote read %d records, local read %d", len(got), len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("record %d: remote %q, local %q", i, got[i], want[i])
		}
	}
}

// A remote CRAM with a REMOTE reference — the fourth combination.
func TestRemoteCramRemoteReference(t *testing.T) {
	ts, _ := serveTestdata(t, "test_raw.cram", "test_raw.cram.crai", "ref.fa", "ref.fa.fai")

	local, err := NewReader(filepath.Join("testdata", "test_raw.cram"), filepath.Join("testdata", "ref.fa"))
	if err != nil {
		t.Fatal(err)
	}
	defer local.Close()
	want := readAll(t, local)

	src, err := iosource.Open(context.Background(), ts.URL+"/test_raw.cram")
	if err != nil {
		t.Fatal(err)
	}
	opts := htsio.NewSamReaderOpts()
	remote, err := NewReaderFromSource(src, ts.URL+"/test_raw.cram", ts.URL+"/ref.fa", opts)
	if err != nil {
		t.Fatalf("remote CRAM with remote reference: %v", err)
	}
	defer remote.Close()

	got := readAll(t, remote)
	if len(got) != len(want) {
		t.Fatalf("remote/remote read %d records, local read %d", len(got), len(want))
	}
}

// The .crai must be resolved through the same transport as the data.
func TestRemoteCramIndexedQuery(t *testing.T) {
	ts, _ := serveTestdata(t, "test_raw.cram", "test_raw.cram.crai")
	localRef := filepath.Join("testdata", "ref.fa")

	local, err := NewReader(filepath.Join("testdata", "test_raw.cram"), localRef)
	if err != nil {
		t.Fatal(err)
	}
	defer local.Close()
	hdr, err := local.Header()
	if err != nil {
		t.Fatal(err)
	}
	var refName string
	for rec, err := range local.Records() {
		if err == nil && rec.RefName != "" && rec.RefName != "*" {
			refName = rec.RefName
			break
		}
	}
	if refName == "" {
		t.Skip("fixture has no mapped records to query")
	}
	_ = hdr

	src, err := iosource.Open(context.Background(), ts.URL+"/test_raw.cram")
	if err != nil {
		t.Fatal(err)
	}
	remote, err := NewReaderFromSource(src, ts.URL+"/test_raw.cram", localRef, nil)
	if err != nil {
		t.Fatal(err)
	}
	defer remote.Close()

	lseq, lerr := local.Query(refName, 0, 1<<30)
	rseq, rerr := remote.Query(refName, 0, 1<<30)
	if (lerr == nil) != (rerr == nil) {
		t.Fatalf("Query error mismatch: local %v, remote %v", lerr, rerr)
	}
	if lerr != nil {
		t.Skipf("fixture does not support indexed query: %v", lerr)
	}
	var ln, rn int
	for rec, err := range lseq {
		if err != nil || rec == nil {
			break
		}
		ln++
	}
	for rec, err := range rseq {
		if err != nil || rec == nil {
			break
		}
		rn++
	}
	if ln != rn {
		t.Errorf("indexed query: remote %d records, local %d", rn, ln)
	}
	if ln == 0 {
		t.Error("indexed query returned nothing; the comparison proves nothing")
	}
}
