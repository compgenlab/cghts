package varstore

import (
	"context"
	"fmt"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"testing"

	s3client "github.com/compgenlab/cghts/iosource/s3"
)

// buildStore writes a small store and returns its base path.
func buildStore(t *testing.T, dir string) string {
	t.Helper()
	base := filepath.Join(dir, "cohort")
	w, err := NewWriter(base, WriterOpts{Samples: []string{"S1", "S2"}, MinDP: 10, RowGroupSize: 500})
	if err != nil {
		t.Fatal(err)
	}
	for _, chrom := range []string{"chr1", "chr17"} {
		for pos := int32(100); pos <= 20000; pos += 100 {
			for _, s := range []string{"S1", "S2"} {
				if err := w.WriteCall(Call{
					SampleID: s, Chrom: chrom, Pos: pos, Ref: "A", Alt: "G",
					GT: "0/1", DP: 30, ADRef: 15, ADAlt: 15, GQ: 99,
				}); err != nil {
					t.Fatal(err)
				}
			}
			if err := w.WriteSite(Site{
				Chrom: chrom, Pos: pos, Ref: "A", Alt: "G",
				AC: 2, AN: 4, NCarriers: 2, NCalled: 2,
			}); err != nil {
				t.Fatal(err)
			}
		}
		if err := w.WriteRegion(CalledSiteRun{
			SampleID: "S1", Chrom: chrom, Start: 100, End: 20000, NSites: 200,
		}); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
	return base
}

func callsAt(t *testing.T, s Store, l Locus) []Call {
	t.Helper()
	got, err := CollectCalls(s, Query{Loci: []Locus{l}})
	if err != nil {
		t.Fatal(err)
	}
	return got
}

// A store served over HTTP must answer exactly as it does from disk. This is
// the whole claim: Parquet's footer statistics let a locus query skip row
// groups, and skipping them remotely means never transferring them.
func TestParquetStoreOverHTTP(t *testing.T) {
	dir := t.TempDir()
	base := buildStore(t, dir)

	ts := httptest.NewServer(http.FileServer(http.Dir(dir)))
	defer ts.Close()

	local, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer local.Close()

	remote, err := OpenParquetContext(context.Background(), ts.URL+"/cohort")
	if err != nil {
		t.Fatalf("open over http: %v", err)
	}
	defer remote.Close()

	ls, err := local.Samples()
	if err != nil {
		t.Fatal(err)
	}
	rs, err := remote.Samples()
	if err != nil {
		t.Fatal(err)
	}
	if fmt.Sprint(ls) != fmt.Sprint(rs) {
		t.Errorf("samples: remote %v, local %v", rs, ls)
	}

	for _, l := range []Locus{
		{Chrom: "chr1", Pos: 5000, Ref: "A", Alt: "G"},
		{Chrom: "chr17", Pos: 19900, Ref: "A", Alt: "G"},
		{Chrom: "chr1", Pos: 12345, Ref: "A", Alt: "G"}, // no such site
	} {
		want := callsAt(t, local, l)
		got := callsAt(t, remote, l)
		if len(want) != len(got) {
			t.Fatalf("%s: %d calls, local gave %d", l, len(got), len(want))
		}
		for i := range want {
			if want[i] != got[i] {
				t.Fatalf("%s call %d: remote %+v, local %+v", l, i, got[i], want[i])
			}
		}
	}

	// Classify exercises all three members together.
	lst, lerr := local.Classify(Locus{Chrom: "chr1", Pos: 5000, Ref: "A", Alt: "G"}, Gate{MinDP: 10})
	rst, rerr := remote.Classify(Locus{Chrom: "chr1", Pos: 5000, Ref: "A", Alt: "G"}, Gate{MinDP: 10})
	if (lerr == nil) != (rerr == nil) {
		t.Fatalf("Classify error mismatch: local %v, remote %v", lerr, rerr)
	}
	if len(lst) != len(rst) {
		t.Fatalf("Classify returned %d states, local gave %d", len(rst), len(lst))
	}
	for i := range lst {
		// Compare by value: SampleState.Call is a pointer, so the structs
		// themselves are never equal across two readers.
		if lst[i].SampleID != rst[i].SampleID || lst[i].State != rst[i].State {
			t.Errorf("state %d: remote %+v, local %+v", i, rst[i], lst[i])
			continue
		}
		switch {
		case lst[i].Call == nil && rst[i].Call == nil:
		case lst[i].Call == nil || rst[i].Call == nil:
			t.Errorf("state %d: one Call is nil (remote %v, local %v)", i, rst[i].Call, lst[i].Call)
		case *lst[i].Call != *rst[i].Call:
			t.Errorf("state %d call: remote %+v, local %+v", i, *rst[i].Call, *lst[i].Call)
		}
	}
	if len(lst) == 0 {
		t.Error("Classify returned nothing; the comparison proves nothing")
	}
}

// The same over s3://, which additionally proves the scheme registry works.
func TestParquetStoreOverS3(t *testing.T) {
	bucket := os.Getenv("VARHUB_TEST_S3_BUCKET")
	if bucket == "" {
		t.Skip("set VARHUB_TEST_S3_BUCKET (and AWS_ENDPOINT_URL) to run S3 integration tests")
	}
	ctx := context.Background()
	dir := t.TempDir()
	base := buildStore(t, dir)

	local, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer local.Close()

	// Upload the three members with the plain AWS CLI-equivalent: a PUT each.
	prefix := "varstore-test"
	for _, m := range []string{CallsMember, SitesMember, RegionsMember} {
		if err := putObject(ctx, bucket, prefix+"/cohort."+m+".parquet",
			MemberPath(base, m)); err != nil {
			t.Fatal(err)
		}
	}

	remote, err := OpenParquetContext(ctx, "s3://"+bucket+"/"+prefix+"/cohort")
	if err != nil {
		t.Fatalf("open over s3: %v", err)
	}
	defer remote.Close()

	for _, l := range []Locus{
		{Chrom: "chr1", Pos: 5000, Ref: "A", Alt: "G"},
		{Chrom: "chr17", Pos: 19900, Ref: "A", Alt: "G"},
	} {
		want := callsAt(t, local, l)
		got := callsAt(t, remote, l)
		if len(want) != len(got) {
			t.Fatalf("%s: %d calls over s3, local gave %d", l, len(got), len(want))
		}
		for i := range want {
			if want[i] != got[i] {
				t.Fatalf("%s call %d: s3 %+v, local %+v", l, i, got[i], want[i])
			}
		}
		if len(want) == 0 {
			t.Fatalf("%s: local returned nothing; the comparison proves nothing", l)
		}
	}
}

// putObject uploads a file, so the s3 test can stage its own fixture without
// varstore gaining a write-side S3 dependency.
func putObject(ctx context.Context, bucket, key, path string) error {
	c, err := s3client.New(ctx)
	if err != nil {
		return err
	}
	_ = c
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()
	return s3client.PutForTest(ctx, bucket, key, f)
}
