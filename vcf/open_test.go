package vcf

import (
	"context"
	"net/http"
	"net/http/httptest"
	"testing"
)

// recordsOf drains a reader into comparable strings.
func recordsOf(t *testing.T, r *VcfReader) []string {
	t.Helper()
	var out []string
	for {
		rec, err := r.NextRecord()
		if err != nil || rec == nil {
			return out
		}
		out = append(out, rec.Chrom+":"+itoa(rec.Pos)+":"+rec.Ref)
	}
}

func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	var b []byte
	for n > 0 {
		b = append([]byte{byte('0' + n%10)}, b...)
		n /= 10
	}
	return string(b)
}

// A locator with no scheme must reach the existing constructor, not a second
// implementation of it. Otherwise the 1 MB buffer, the gzip magic test and the
// close-on-error path get re-derived here and are free to drift from the
// originals -- silently, since both would still parse a well-formed VCF.
func TestOpenVcfFileMatchesNewVcfFile(t *testing.T) {
	for _, path := range []string{"testdata/sample.vcf", "testdata/sample.vcf.gz"} {
		want, err := NewVcfFile(path)
		if err != nil {
			t.Fatal(err)
		}
		wantRecs := recordsOf(t, want)
		want.Close()
		if len(wantRecs) == 0 {
			t.Fatalf("%s yielded no records; the fixture proves nothing", path)
		}

		got, err := OpenVcfFile(context.Background(), path)
		if err != nil {
			t.Fatal(err)
		}
		gotRecs := recordsOf(t, got)
		got.Close()

		if len(gotRecs) != len(wantRecs) {
			t.Fatalf("%s: OpenVcfFile read %d records, NewVcfFile read %d",
				path, len(gotRecs), len(wantRecs))
		}
		for i := range wantRecs {
			if gotRecs[i] != wantRecs[i] {
				t.Errorf("%s record %d: got %s, want %s", path, i, gotRecs[i], wantRecs[i])
			}
		}
	}
}

func TestOpenVcfFileOverHTTP(t *testing.T) {
	srv := httptest.NewServer(http.FileServer(http.Dir("testdata")))
	defer srv.Close()

	local, err := NewVcfFile("testdata/sample.vcf.gz")
	if err != nil {
		t.Fatal(err)
	}
	want := recordsOf(t, local)
	local.Close()

	r, err := OpenVcfFile(context.Background(), srv.URL+"/sample.vcf.gz")
	if err != nil {
		t.Fatal(err)
	}
	got := recordsOf(t, r)
	r.Close()

	if len(got) != len(want) {
		t.Fatalf("over http read %d records, locally %d", len(got), len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("record %d: got %s, want %s", i, got[i], want[i])
		}
	}
}

func TestOpenIndexedVcfReaderOverHTTP(t *testing.T) {
	srv := httptest.NewServer(http.FileServer(http.Dir("testdata")))
	defer srv.Close()

	local, err := NewIndexedVcfReader("testdata/sample.vcf.gz")
	if err != nil {
		t.Fatal(err)
	}
	defer local.Close()
	lh, err := local.Header()
	if err != nil {
		t.Fatal(err)
	}

	remote, err := OpenIndexedVcfReader(context.Background(), srv.URL+"/sample.vcf.gz")
	if err != nil {
		t.Fatal(err)
	}
	defer remote.Close()
	rh, err := remote.Header()
	if err != nil {
		t.Fatal(err)
	}

	if got, want := len(rh.Samples()), len(lh.Samples()); got != want {
		t.Fatalf("remote header has %d samples, local %d", got, want)
	}
	for i, s := range lh.Samples() {
		if rh.Samples()[i] != s {
			t.Errorf("sample %d: got %s, want %s", i, rh.Samples()[i], s)
		}
	}
}

// A remote file with no index must say so in a way that names both suffixes it
// looked for; a bare 404 for ".csi" reads as the wrong index kind expected.
func TestIndexedOpenNamesEverySuffixTried(t *testing.T) {
	mux := http.NewServeMux()
	mux.HandleFunc("/only.vcf.gz", func(w http.ResponseWriter, r *http.Request) {
		http.ServeFile(w, r, "testdata/sample.vcf.gz")
	})
	srv := httptest.NewServer(mux)
	defer srv.Close()

	_, err := OpenIndexedVcfReader(context.Background(), srv.URL+"/only.vcf.gz")
	if err == nil {
		t.Fatal("an unindexed remote VCF opened for random access")
	}
	msg := err.Error()
	for _, want := range []string{".tbi", ".csi"} {
		if !contains(msg, want) {
			t.Errorf("error does not mention %s: %v", want, err)
		}
	}
}

func contains(h, n string) bool {
	for i := 0; i+len(n) <= len(h); i++ {
		if h[i:i+len(n)] == n {
			return true
		}
	}
	return false
}
