package vcf

import (
	"os"
	"reflect"
	"testing"

	"github.com/compgenlab/cghts/htsio/tabix"
)

func openIndexed(t *testing.T) *IndexedVcfReader {
	t.Helper()
	r, err := NewIndexedVcfReader("testdata/sample.vcf.gz")
	if err != nil {
		t.Fatalf("NewIndexedVcfReader: %v", err)
	}
	t.Cleanup(func() { r.Close() })
	return r
}

// RefNames comes from the index, not the header, and is what a contig converter
// is built from -- tabix matches reference names verbatim, so querying "22"
// against a "chr22" index fails with "unknown reference" unless the caller
// resolves the spelling against this list first.
func TestIndexedRefNames(t *testing.T) {
	r := openIndexed(t)
	got := r.RefNames()
	if want := []string{"chr1", "chr2"}; !reflect.DeepEqual(got, want) {
		t.Errorf("RefNames() = %v, want %v", got, want)
	}
}

// HasRef is the cheap pre-check that turns an absent contig into an empty result
// instead of an error -- a contig a file lacks is an absence, not a failure.
func TestIndexedHasRef(t *testing.T) {
	r := openIndexed(t)
	for _, tc := range []struct {
		ref  string
		want bool
	}{
		{"chr1", true},
		{"chr2", true},
		{"chr3", false},
		// Verbatim, not canonical: HasRef reports what the index literally
		// holds, which is precisely why a converter is needed above it.
		{"1", false},
		{"CHR1", false},
		{"", false},
	} {
		if got := r.HasRef(tc.ref); got != tc.want {
			t.Errorf("HasRef(%q) = %v, want %v", tc.ref, got, tc.want)
		}
	}
}

// NewIndexedVcfReaderFrom is the route for a VCF that is not a local file -- an
// HTTP-Range source or an S3 object. It has to behave identically to the
// by-filename constructor, since that is the whole point of remote inputs.
func TestNewIndexedVcfReaderFrom(t *testing.T) {
	data, err := os.Open("testdata/sample.vcf.gz")
	if err != nil {
		t.Fatal(err)
	}
	defer data.Close()
	index, err := os.Open("testdata/sample.vcf.gz.tbi")
	if err != nil {
		t.Fatal(err)
	}
	defer index.Close()

	tr, err := tabix.NewReaderFromSource(data, index)
	if err != nil {
		t.Fatalf("tabix.NewReaderFromSource: %v", err)
	}
	r := NewIndexedVcfReaderFrom(tr, "s3://bucket/sample.vcf.gz")
	defer r.Close()

	h, err := r.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if got := h.Samples(); !reflect.DeepEqual(got, []string{"NORMAL", "TUMOR"}) {
		t.Errorf("Samples() = %v, want [NORMAL TUMOR]", got)
	}
	if got := r.RefNames(); !reflect.DeepEqual(got, []string{"chr1", "chr2"}) {
		t.Errorf("RefNames() = %v, want [chr1 chr2]", got)
	}
	if !r.HasRef("chr1") {
		t.Error("HasRef(chr1) = false on a source-backed reader")
	}

	// A query returns the same rows the local constructor does.
	seq, err := r.Query("chr1", 0, 1000)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}
	var pos []int
	for rec, err := range seq {
		if err != nil {
			t.Fatalf("Query iteration: %v", err)
		}
		pos = append(pos, rec.Pos)
	}
	if want := []int{100, 200, 300}; !reflect.DeepEqual(pos, want) {
		t.Errorf("positions = %v, want %v", pos, want)
	}
}

// The label is provenance only -- it never has to be openable, which is what
// lets a remote locator stand in for a filename in error messages.
func TestIndexedReaderFromUsesLabelInErrors(t *testing.T) {
	data, err := os.Open("testdata/sample.vcf.gz")
	if err != nil {
		t.Fatal(err)
	}
	defer data.Close()
	index, err := os.Open("testdata/sample.vcf.gz.tbi")
	if err != nil {
		t.Fatal(err)
	}
	defer index.Close()

	tr, err := tabix.NewReaderFromSource(data, index)
	if err != nil {
		t.Fatal(err)
	}
	r := NewIndexedVcfReaderFrom(tr, "https://example.invalid/x.vcf.gz")
	defer r.Close()

	// It reads fine despite the label naming nothing that exists locally.
	if _, err := r.Header(); err != nil {
		t.Errorf("Header on a labelled source reader: %v", err)
	}
}
