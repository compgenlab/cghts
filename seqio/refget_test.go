package seqio

import (
	"fmt"
	"net/http"
	"net/http/httptest"
	"strings"
	"testing"
)

func TestRefgetReader(t *testing.T) {
	// Create a fake refget server with known sequences.
	sequences := map[string]string{
		"abc123": "ACGTACGTNNACGTACGT",
		"def456": "TTTTGGGGCCCCAAAA",
	}

	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		path := r.URL.Path
		// Handle /md5/{hash} and /md5/{hash}/metadata
		parts := strings.Split(strings.TrimPrefix(path, "/"), "/")
		if len(parts) < 1 {
			http.NotFound(w, r)
			return
		}

		md5 := parts[0]
		if len(parts) == 2 && parts[1] == "metadata" {
			seq, ok := sequences[md5]
			if !ok {
				http.NotFound(w, r)
				return
			}
			w.Header().Set("Content-Type", "application/json")
			fmt.Fprintf(w, `{"metadata":{"length":%d,"md5":"%s"}}`, len(seq), md5)
			return
		}

		seq, ok := sequences[md5]
		if !ok {
			http.NotFound(w, r)
			return
		}

		start := 0
		end := len(seq)
		q := r.URL.Query()
		if s := q.Get("start"); s != "" {
			fmt.Sscanf(s, "%d", &start)
		}
		if e := q.Get("end"); e != "" {
			fmt.Sscanf(e, "%d", &end)
		}

		if start > len(seq) {
			start = len(seq)
		}
		if end > len(seq) {
			end = len(seq)
		}

		w.Header().Set("Content-Type", "text/plain")
		w.Write([]byte(seq[start:end]))
	}))
	defer srv.Close()

	md5s := map[string]string{
		"chr1": "abc123",
		"chr2": "def456",
	}
	lengths := map[string]int{
		"chr1": 18,
		"chr2": 16,
	}
	names := []string{"chr1", "chr2"}

	rr, err := NewRefgetReader(md5s, lengths, names, RefgetServer(srv.URL))
	if err != nil {
		t.Fatalf("NewRefgetReader: %v", err)
	}
	defer rr.Close()

	// Names.
	if got := rr.Names(); len(got) != 2 || got[0] != "chr1" || got[1] != "chr2" {
		t.Errorf("Names() = %v", got)
	}

	// SequenceLength.
	if l, ok := rr.SequenceLength("chr1"); !ok || l != 18 {
		t.Errorf("SequenceLength(chr1) = %d, %v", l, ok)
	}

	// Full sequence.
	seq, err := rr.GetSequence("chr1")
	if err != nil {
		t.Fatalf("GetSequence(chr1): %v", err)
	}
	if string(seq) != "ACGTACGTNNACGTACGT" {
		t.Errorf("GetSequence(chr1) = %q", string(seq))
	}

	// Range.
	sub, err := rr.GetSequenceRange("chr1", 4, 8)
	if err != nil {
		t.Fatalf("GetSequenceRange: %v", err)
	}
	if string(sub) != "ACGT" {
		t.Errorf("GetSequenceRange(4,8) = %q, want ACGT", string(sub))
	}

	// Missing sequence.
	_, err = rr.GetSequenceRange("chrX", 0, 10)
	if err == nil {
		t.Error("expected error for missing sequence")
	}
}

func TestRefgetReaderNoMD5s(t *testing.T) {
	_, err := NewRefgetReader(nil, nil, nil)
	if err == nil {
		t.Error("expected error with empty MD5 map")
	}
}

func TestParseRefgetLength(t *testing.T) {
	tests := []struct {
		body string
		want int
	}{
		{`{"metadata":{"length":12345,"md5":"abc"}}`, 12345},
		{`{"length": 999}`, 999},
	}
	for _, tc := range tests {
		got, err := parseRefgetLength(tc.body)
		if err != nil {
			t.Errorf("parseRefgetLength(%q): %v", tc.body, err)
		}
		if got != tc.want {
			t.Errorf("parseRefgetLength(%q) = %d, want %d", tc.body, got, tc.want)
		}
	}
}
