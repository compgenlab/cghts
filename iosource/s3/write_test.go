package s3

import (
	"bytes"
	"context"
	"fmt"
	"io"
	"net/http"
	"net/http/httptest"
	"sort"
	"strconv"
	"strings"
	"sync"
	"testing"
)

// The write path, against an in-process server that speaks enough of the S3 REST
// API to run a real multipart upload. Same reasoning as the read tests: the
// requests are really signed, really sent and really parsed, so what is under
// test is this package's use of the SDK rather than a restatement of it.
//
// Multipart is the part worth testing at all. A single PutObject is one call
// that either works or does not; a multipart upload is a protocol -- initiate,
// N parts each returning an ETag, complete with those ETags in ascending order
// -- and every one of its failure modes shows up only after the bytes have
// moved.

type fakeWriteS3 struct {
	mu      sync.Mutex
	objects map[string][]byte          // "bucket/key" -> contents
	uploads map[string]map[int32][]byte // uploadId -> part number -> bytes

	initiated int
	completed int
	aborted   int
	simplePut int
}

func newFakeWriteS3() *fakeWriteS3 {
	return &fakeWriteS3{
		objects: map[string][]byte{},
		uploads: map[string]map[int32][]byte{},
	}
}

func (f *fakeWriteS3) ServeHTTP(w http.ResponseWriter, r *http.Request) {
	name := strings.TrimPrefix(r.URL.Path, "/")
	q := r.URL.Query()

	f.mu.Lock()
	defer f.mu.Unlock()

	switch {
	// Initiate.
	case r.Method == http.MethodPost && q.Has("uploads"):
		f.initiated++
		id := fmt.Sprintf("upload-%d", f.initiated)
		f.uploads[id] = map[int32][]byte{}
		bucket, key, _ := strings.Cut(name, "/")
		fmt.Fprintf(w, `<?xml version="1.0" encoding="UTF-8"?>`+
			`<InitiateMultipartUploadResult><Bucket>%s</Bucket><Key>%s</Key><UploadId>%s</UploadId></InitiateMultipartUploadResult>`,
			bucket, key, id)

	// One part.
	case r.Method == http.MethodPut && q.Get("uploadId") != "":
		id := q.Get("uploadId")
		n, _ := strconv.Atoi(q.Get("partNumber"))
		body, _ := io.ReadAll(r.Body)
		if f.uploads[id] == nil {
			s3Error(w, http.StatusNotFound, "NoSuchUpload", "no such upload")
			return
		}
		f.uploads[id][int32(n)] = body
		w.Header().Set("ETag", fmt.Sprintf(`"part-%d"`, n))
		w.WriteHeader(http.StatusOK)

	// Complete: assemble in part order.
	case r.Method == http.MethodPost && q.Get("uploadId") != "":
		id := q.Get("uploadId")
		parts, ok := f.uploads[id]
		if !ok {
			s3Error(w, http.StatusNotFound, "NoSuchUpload", "no such upload")
			return
		}
		nums := make([]int, 0, len(parts))
		for n := range parts {
			nums = append(nums, int(n))
		}
		sort.Ints(nums)
		var whole []byte
		for _, n := range nums {
			whole = append(whole, parts[int32(n)]...)
		}
		f.objects[name] = whole
		delete(f.uploads, id)
		f.completed++
		bucket, key, _ := strings.Cut(name, "/")
		fmt.Fprintf(w, `<?xml version="1.0" encoding="UTF-8"?>`+
			`<CompleteMultipartUploadResult><Bucket>%s</Bucket><Key>%s</Key><ETag>"whole"</ETag></CompleteMultipartUploadResult>`,
			bucket, key)

	// Abandon.
	case r.Method == http.MethodDelete && q.Get("uploadId") != "":
		delete(f.uploads, q.Get("uploadId"))
		f.aborted++
		w.WriteHeader(http.StatusNoContent)

	case r.Method == http.MethodPut:
		f.simplePut++
		body, _ := io.ReadAll(r.Body)
		f.objects[name] = body
		w.WriteHeader(http.StatusOK)

	case r.Method == http.MethodDelete:
		delete(f.objects, name)
		w.WriteHeader(http.StatusNoContent)

	case r.Method == http.MethodHead:
		body, ok := f.objects[name]
		if !ok {
			s3Error(w, http.StatusNotFound, "NotFound", "no such key")
			return
		}
		w.Header().Set("Content-Length", strconv.Itoa(len(body)))
		w.WriteHeader(http.StatusOK)

	default:
		s3Error(w, http.StatusNotFound, "NoSuchKey", "no such key")
	}
}

func writeClient(t *testing.T, f *fakeWriteS3) *Client {
	t.Helper()
	srv := httptest.NewServer(f)
	t.Cleanup(srv.Close)

	// Static credentials, isolated from whatever the developer has in ~/.aws.
	// Without these the SDK reaches for the ambient chain and stalls on
	// instance metadata for five seconds before failing.
	t.Setenv("AWS_ACCESS_KEY_ID", "testkey")
	t.Setenv("AWS_SECRET_ACCESS_KEY", "testsecret")
	t.Setenv("AWS_SESSION_TOKEN", "")
	t.Setenv("AWS_REGION", "us-east-1")
	t.Setenv("AWS_PROFILE", "")
	t.Setenv("AWS_CONFIG_FILE", t.TempDir()+"/config")
	t.Setenv("AWS_SHARED_CREDENTIALS_FILE", t.TempDir()+"/credentials")
	t.Setenv("AWS_EC2_METADATA_DISABLED", "true")

	c, err := New(context.Background(),
		WithEndpoint(srv.URL), WithRegion("us-east-1"), WithHTTPClient(srv.Client()))
	if err != nil {
		t.Fatal(err)
	}
	return c
}

// An object small enough to fit in one buffer must not create a multipart
// upload at all.
//
// Not an optimisation: an abandoned multipart upload is invisible to a bucket
// listing and is billed until something removes it, so the cheapest way to
// handle one is not to start it. A manifest is a few kilobytes and there is one
// per store.
func TestASmallObjectIsOnePut(t *testing.T) {
	f := newFakeWriteS3()
	c := writeClient(t, f)

	w, err := c.Create(context.Background(), "s3://bucket/small.txt")
	if err != nil {
		t.Fatal(err)
	}
	if _, err := io.WriteString(w, "hello"); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	if f.initiated != 0 {
		t.Errorf("started %d multipart upload(s) for five bytes", f.initiated)
	}
	if f.simplePut != 1 {
		t.Errorf("made %d PutObject calls, want 1", f.simplePut)
	}
	if got := string(f.objects["bucket/small.txt"]); got != "hello" {
		t.Errorf("stored %q", got)
	}
}

// Past one part, the object goes up in parts and is reassembled in order.
func TestALargeObjectIsUploadedInParts(t *testing.T) {
	f := newFakeWriteS3()
	c := writeClient(t, f)

	w, err := c.Create(context.Background(), "s3://bucket/big.bin", WithPartSize(5<<20))
	if err != nil {
		t.Fatal(err)
	}
	// Two and a half parts, written in chunks that do not align with the part
	// size -- a writer that only worked on aligned writes would pass a tidier
	// test and fail on a real Parquet stream.
	want := make([]byte, 0, 13<<20)
	chunk := bytes.Repeat([]byte("abcdefg"), 1024)
	for len(want) < 13<<20 {
		if _, err := w.Write(chunk); err != nil {
			t.Fatal(err)
		}
		want = append(want, chunk...)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	if f.initiated != 1 || f.completed != 1 {
		t.Errorf("initiated %d and completed %d uploads, want 1 and 1", f.initiated, f.completed)
	}
	got := f.objects["bucket/big.bin"]
	if len(got) != len(want) {
		t.Fatalf("stored %d bytes, wrote %d", len(got), len(want))
	}
	if !bytes.Equal(got, want) {
		t.Error("the reassembled object does not match what was written; parts are out of order or overlapping")
	}
	if w.Written() != int64(len(want)) {
		t.Errorf("Written() reports %d, wrote %d", w.Written(), len(want))
	}
}

// Abandoning must leave nothing: no object, and no upload still holding parts.
func TestAbortLeavesNothingBehind(t *testing.T) {
	f := newFakeWriteS3()
	c := writeClient(t, f)

	w, err := c.Create(context.Background(), "s3://bucket/doomed.bin", WithPartSize(5<<20))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := w.Write(bytes.Repeat([]byte("x"), 6<<20)); err != nil {
		t.Fatal(err)
	}
	if f.initiated != 1 {
		t.Fatal("the fixture proves nothing: no multipart upload was started")
	}
	if err := w.Abort(); err != nil {
		t.Fatal(err)
	}

	if _, exists := f.objects["bucket/doomed.bin"]; exists {
		t.Error("an abandoned upload produced an object")
	}
	if f.aborted != 1 {
		t.Errorf("aborted %d uploads, want 1", f.aborted)
	}
	if len(f.uploads) != 0 {
		t.Errorf("%d upload(s) still holding parts; they are invisible in a listing and still billed", len(f.uploads))
	}
}

// Stat has to distinguish "not there" from "could not tell", because the
// overwrite guard treats the two differently: absent means go ahead.
func TestStatSaysAbsentRatherThanFailing(t *testing.T) {
	f := newFakeWriteS3()
	c := writeClient(t, f)

	if _, ok, err := c.Stat(context.Background(), "s3://bucket/missing"); err != nil || ok {
		t.Fatalf("Stat on a missing key returned ok=%v err=%v; want false, nil", ok, err)
	}

	w, _ := c.Create(context.Background(), "s3://bucket/present")
	io.WriteString(w, "1234567890")
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	size, ok, err := c.Stat(context.Background(), "s3://bucket/present")
	if err != nil || !ok {
		t.Fatalf("Stat on a present key returned ok=%v err=%v", ok, err)
	}
	if size != 10 {
		t.Errorf("size %d, want 10", size)
	}
}
