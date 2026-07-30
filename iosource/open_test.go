package iosource

import (
	"context"
	"io"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestSchemeOf(t *testing.T) {
	for _, tc := range []struct{ in, want string }{
		{"/var/data/x.vcf.gz", ""},
		{"relative/x.gz", ""},
		{"s3://bucket/key", "s3"},
		{"S3://bucket/key", "s3"},
		{"https://host/x", "https"},
		{"http://host/x", "http"},
		{`C:\data\x.gz`, ""}, // a drive letter is not a scheme
	} {
		if got := schemeOf(tc.in); got != tc.want {
			t.Errorf("schemeOf(%q) = %q, want %q", tc.in, got, tc.want)
		}
	}
}

func TestOpenDispatch(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "d.bin")
	body := []byte("varhub varhub varhub")
	if err := os.WriteFile(path, body, 0o644); err != nil {
		t.Fatal(err)
	}
	ctx := context.Background()

	// A plain path opens a file.
	src, err := Open(ctx, path)
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()
	if n, _ := src.Size(); n != int64(len(body)) {
		t.Errorf("file size = %d, want %d", n, len(body))
	}

	// http(s) is built in.
	ts := httptest.NewServer(http.FileServer(http.Dir(dir)))
	defer ts.Close()
	hsrc, err := Open(ctx, ts.URL+"/d.bin")
	if err != nil {
		t.Fatal(err)
	}
	defer hsrc.Close()
	p := make([]byte, 6)
	if _, err := hsrc.ReadAt(p, 7); err != nil {
		t.Fatal(err)
	}
	if string(p) != string(body[7:13]) {
		t.Errorf("http read = %q", p)
	}
}

// An unregistered scheme must say so, and say what to do about it — the whole
// point of the registry is that the transport lives in another package.
func TestOpenUnregisteredScheme(t *testing.T) {
	_, err := Open(context.Background(), "gs://bucket/object")
	if err == nil {
		t.Fatal("an unregistered scheme was accepted")
	}
	if !strings.Contains(err.Error(), "gs://") {
		t.Errorf("error does not name the scheme: %v", err)
	}
}

func TestRegisterAndOpen(t *testing.T) {
	Register("teststub", func(_ context.Context, locator string) (ByteSource, error) {
		return &stubSource{data: []byte(strings.TrimPrefix(locator, "teststub://"))}, nil
	})
	src, err := Open(context.Background(), "teststub://hello")
	if err != nil {
		t.Fatal(err)
	}
	n, _ := src.Size()
	if n != 5 {
		t.Errorf("stub size = %d, want 5", n)
	}
	found := false
	for _, s := range Schemes() {
		if s == "teststub" {
			found = true
		}
	}
	if !found {
		t.Error("registered scheme not listed by Schemes()")
	}
}

type stubSource struct{ data []byte }

func (s *stubSource) ReadAt(p []byte, off int64) (int, error) {
	if off >= int64(len(s.data)) {
		return 0, io.EOF
	}
	n := copy(p, s.data[off:])
	if n < len(p) {
		return n, io.EOF
	}
	return n, nil
}
func (s *stubSource) Size() (int64, error) { return int64(len(s.data)), nil }
func (s *stubSource) Close() error         { return nil }
