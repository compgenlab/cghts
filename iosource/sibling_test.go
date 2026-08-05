package iosource

import (
	"errors"
	"io"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestFileSibling(t *testing.T) {
	dir := t.TempDir()
	data := filepath.Join(dir, "x.vcf.gz")
	if err := os.WriteFile(data, []byte("data"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(data+".tbi", []byte("index"), 0o644); err != nil {
		t.Fatal(err)
	}
	rc, err := FileSibling(data, ".tbi")
	if err != nil {
		t.Fatal(err)
	}
	defer rc.Close()
	got, _ := io.ReadAll(rc)
	if string(got) != "index" {
		t.Errorf("read %q, want %q", got, "index")
	}
}

// The error must name every suffix tried, not just the last failure. Reporting
// only the last one meant a missing index complained about ".csi" and never
// mentioned ".tbi", which reads as though the wrong index kind was expected
// when in fact neither was there.
func TestResolveSiblingNamesEverySuffixTried(t *testing.T) {
	dir := t.TempDir()
	data := filepath.Join(dir, "x.vcf.gz")
	if err := os.WriteFile(data, []byte("data"), 0o644); err != nil {
		t.Fatal(err)
	}
	_, _, err := ResolveSibling(data, []string{".tbi", ".csi"}, FileSibling)
	if err == nil {
		t.Fatal("ResolveSibling succeeded with no sibling present")
	}
	for _, want := range []string{".tbi", ".csi", "x.vcf.gz"} {
		if !strings.Contains(err.Error(), want) {
			t.Errorf("error %q does not mention %q", err, want)
		}
	}
}

// First match wins, in the order given -- so a caller can express a preference
// between index kinds rather than getting whichever the filesystem returns.
func TestResolveSiblingPrefersTheEarlierSuffix(t *testing.T) {
	dir := t.TempDir()
	data := filepath.Join(dir, "x.vcf.gz")
	for _, f := range []struct{ suffix, body string }{{".tbi", "tbi"}, {".csi", "csi"}} {
		if err := os.WriteFile(data+f.suffix, []byte(f.body), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	if err := os.WriteFile(data, []byte("data"), 0o644); err != nil {
		t.Fatal(err)
	}

	rc, suffix, err := ResolveSibling(data, []string{".csi", ".tbi"}, FileSibling)
	if err != nil {
		t.Fatal(err)
	}
	defer rc.Close()
	if suffix != ".csi" {
		t.Errorf("suffix = %q, want the first one listed (.csi)", suffix)
	}
	got, _ := io.ReadAll(rc)
	if string(got) != "csi" {
		t.Errorf("read %q, want %q", got, "csi")
	}
}

func TestHTTPSibling(t *testing.T) {
	ts := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		if strings.HasSuffix(r.URL.Path, ".tbi") {
			w.Write([]byte("index"))
			return
		}
		http.Error(w, "nope", http.StatusNotFound)
	}))
	defer ts.Close()

	rc, suffix, err := ResolveSibling(ts.URL+"/x.vcf.gz", []string{".csi", ".tbi"}, HTTPSibling)
	if err != nil {
		t.Fatal(err)
	}
	defer rc.Close()
	if suffix != ".tbi" {
		t.Errorf("suffix = %q, want .tbi -- .csi is absent", suffix)
	}
	got, _ := io.ReadAll(rc)
	if string(got) != "index" {
		t.Errorf("read %q, want %q", got, "index")
	}
}

func TestHTTPSiblingReportsAbsence(t *testing.T) {
	ts := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		http.Error(w, "nope", http.StatusNotFound)
	}))
	defer ts.Close()

	_, _, err := ResolveSibling(ts.URL+"/x.vcf.gz", []string{".tbi"}, HTTPSibling)
	if err == nil {
		t.Fatal("ResolveSibling succeeded against a server with no sibling")
	}
	if errors.Is(err, io.EOF) {
		t.Errorf("absence should not surface as EOF: %v", err)
	}
}
