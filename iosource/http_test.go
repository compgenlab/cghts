package iosource

import (
	"errors"
	"io"
	"io/fs"
	"net/http"
	"net/http/httptest"
	"strings"
	"testing"
	"time"
)

// rangeServer serves body over HTTP with Range support, and reports how the
// client behaved so a test can assert on the request pattern rather than only
// on the bytes.
type rangeServer struct {
	body       string
	headStatus int // 0 means 200
	heads      int
	gets       int
}

func (s *rangeServer) handler() http.Handler {
	return http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		if r.Method == http.MethodHead {
			s.heads++
			if s.headStatus != 0 && s.headStatus != http.StatusOK {
				w.WriteHeader(s.headStatus)
				return
			}
		} else {
			s.gets++
		}
		http.ServeContent(w, r, "f.bin", timeZero, strings.NewReader(s.body))
	})
}

func TestHTTPRangeReadsAndSizes(t *testing.T) {
	srv := &rangeServer{body: "0123456789"}
	ts := httptest.NewServer(srv.handler())
	defer ts.Close()

	h, err := NewHTTPRange(ts.URL + "/f.bin")
	if err != nil {
		t.Fatal(err)
	}
	if n, err := h.Size(); err != nil || n != 10 {
		t.Errorf("Size() = (%d, %v), want (10, nil)", n, err)
	}
	// The size came from the HEAD, so no ranged GET was needed for it.
	if srv.heads != 1 {
		t.Errorf("heads = %d, want 1", srv.heads)
	}

	p := make([]byte, 4)
	if _, err := h.ReadAt(p, 3); err != nil && err != io.EOF {
		t.Fatalf("ReadAt: %v", err)
	}
	if string(p) != "3456" {
		t.Errorf("ReadAt(4, 3) = %q, want %q", p, "3456")
	}
}

// The bug: every HEAD failure was swallowed, so Open on a missing URL succeeded
// and the error surfaced later from whatever first tried to read -- as a
// confusing complaint about the file's contents rather than a plain statement
// that the URL is not there.
func TestHTTPRangeFailsOnMissingURL(t *testing.T) {
	ts := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		http.Error(w, "nope", http.StatusNotFound)
	}))
	defer ts.Close()

	h, err := NewHTTPRange(ts.URL + "/missing.bin")
	if err == nil {
		t.Fatalf("NewHTTPRange succeeded on a 404 (got %v)", h)
	}
	// Absence reads the same way as a missing local file, so a caller that
	// treats "not there" as an answer need not know its transport.
	if !errors.Is(err, fs.ErrNotExist) {
		t.Errorf("error %v does not satisfy errors.Is(err, fs.ErrNotExist)", err)
	}
	if !strings.Contains(err.Error(), "missing.bin") {
		t.Errorf("error %q should name the URL", err)
	}
}

func TestHTTPRangeFailsOnServerError(t *testing.T) {
	ts := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		http.Error(w, "boom", http.StatusInternalServerError)
	}))
	defer ts.Close()

	_, err := NewHTTPRange(ts.URL + "/f.bin")
	if err == nil {
		t.Fatal("NewHTTPRange succeeded on a 500")
	}
	// A server error is not an absence, and must not be mistaken for one --
	// varstore treats a missing optional member as "written empty".
	if errors.Is(err, fs.ErrNotExist) {
		t.Errorf("a 500 should not read as fs.ErrNotExist: %v", err)
	}
}

// HEAD is not universally implemented, and refusing the method says nothing
// about the resource. Those servers must still work, with the length recovered
// from the first ranged GET.
func TestHTTPRangeToleratesServersWithoutHEAD(t *testing.T) {
	for _, status := range []int{http.StatusMethodNotAllowed, http.StatusNotImplemented} {
		srv := &rangeServer{body: "0123456789", headStatus: status}
		ts := httptest.NewServer(srv.handler())

		h, err := NewHTTPRange(ts.URL + "/f.bin")
		if err != nil {
			ts.Close()
			t.Fatalf("status %d: NewHTTPRange failed: %v", status, err)
		}
		p := make([]byte, 4)
		if _, err := h.ReadAt(p, 0); err != nil && err != io.EOF {
			ts.Close()
			t.Fatalf("status %d: ReadAt: %v", status, err)
		}
		if string(p) != "0123" {
			t.Errorf("status %d: read %q, want %q", status, p, "0123")
		}
		if n, err := h.Size(); err != nil || n != 10 {
			t.Errorf("status %d: Size() = (%d, %v), want (10, nil); it should be "+
				"recovered from the ranged GET", status, n, err)
		}
		ts.Close()
	}
}

func TestHTTPRangeFailsOnUnreachableHost(t *testing.T) {
	// Reserved TEST-NET-1, and a port nothing listens on.
	_, err := NewHTTPRange("http://192.0.2.1:9/f.bin", WithClient(&http.Client{Timeout: shortTimeout}))
	if err == nil {
		t.Fatal("NewHTTPRange succeeded against an unreachable host")
	}
}

var timeZero = time.Time{}

const shortTimeout = 2 * time.Second
