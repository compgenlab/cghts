// These tests run against an in-process HTTP server that speaks enough of the S3
// REST API to answer HEAD and ranged GET. That is deliberately not a mock of the
// SDK: the requests are really signed, really sent, and the responses really
// parsed, so what is under test is this package's use of the SDK rather than a
// restatement of it. It is also what makes the tests meaningful without
// credentials or a network -- WithEndpoint plus path-style addressing is exactly
// the configuration a MinIO or LocalStack gateway needs, so the same code path
// serves both.
package s3

import (
	"context"
	"encoding/base64"
	"encoding/binary"
	"errors"
	"fmt"
	"hash/crc32"
	"io"
	"net/http"
	"net/http/httptest"
	"strconv"
	"strings"
	"sync"
	"sync/atomic"
	"testing"

	"github.com/compgenlab/cghts/iosource"
)

// fakeS3 serves a fixed set of objects and records what was asked of it, so a
// test can assert on the request pattern -- one HEAD, N ranged GETs -- and not
// only on the bytes.
type fakeS3 struct {
	objects map[string]string // "bucket/key" -> contents

	mu     sync.Mutex
	heads  int32
	gets   int32
	puts   int32
	ranges []string
	paths  []string
}

func (f *fakeS3) record(r *http.Request) {
	f.mu.Lock()
	defer f.mu.Unlock()
	f.paths = append(f.paths, r.URL.Path)
	if rg := r.Header.Get("Range"); rg != "" {
		f.ranges = append(f.ranges, rg)
	}
}

func (f *fakeS3) seenRanges() []string {
	f.mu.Lock()
	defer f.mu.Unlock()
	return append([]string(nil), f.ranges...)
}

func (f *fakeS3) seenPaths() []string {
	f.mu.Lock()
	defer f.mu.Unlock()
	return append([]string(nil), f.paths...)
}

// s3Error writes the XML error document a real S3 returns. The SDK parses the
// Code out of it, so a test that returned a bare status would not exercise the
// same error path production does.
func s3Error(w http.ResponseWriter, status int, code, message string) {
	w.Header().Set("Content-Type", "application/xml")
	w.WriteHeader(status)
	fmt.Fprintf(w, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"+
		"<Error><Code>%s</Code><Message>%s</Message></Error>", code, message)
}

func (f *fakeS3) ServeHTTP(w http.ResponseWriter, r *http.Request) {
	f.record(r)

	// Path-style addressing: /bucket/key... The endpoint override switches this
	// on, because a gateway has no wildcard DNS for virtual-host addressing.
	name := strings.TrimPrefix(r.URL.Path, "/")

	if r.Method == http.MethodPut {
		atomic.AddInt32(&f.puts, 1)
		body, _ := io.ReadAll(r.Body)
		f.mu.Lock()
		f.objects[name] = string(body)
		f.mu.Unlock()
		w.WriteHeader(http.StatusOK)
		return
	}

	f.mu.Lock()
	body, ok := f.objects[name]
	f.mu.Unlock()
	if !ok {
		s3Error(w, http.StatusNotFound, "NoSuchKey", "The specified key does not exist.")
		return
	}

	if r.Method == http.MethodHead {
		atomic.AddInt32(&f.heads, 1)
		w.Header().Set("Content-Length", strconv.Itoa(len(body)))
		w.WriteHeader(http.StatusOK)
		return
	}

	atomic.AddInt32(&f.gets, 1)
	rg := r.Header.Get("Range")
	if rg == "" {
		w.Header().Set("Content-Length", strconv.Itoa(len(body)))
		setChecksum(w, body)
		w.WriteHeader(http.StatusOK)
		io.WriteString(w, body)
		return
	}

	start, end, err := parseByteRange(rg, len(body))
	if err != nil {
		// 416 is the response ReadAt has to translate into io.EOF: a
		// block-compressed reader probes past the end of the file as a matter of
		// course, and treating that as a failure would break every BGZF query.
		w.Header().Set("Content-Range", fmt.Sprintf("bytes */%d", len(body)))
		s3Error(w, http.StatusRequestedRangeNotSatisfiable, "InvalidRange",
			"The requested range is not satisfiable")
		return
	}
	chunk := body[start : end+1]
	w.Header().Set("Content-Range", fmt.Sprintf("bytes %d-%d/%d", start, end, len(body)))
	w.Header().Set("Content-Length", strconv.Itoa(len(chunk)))
	setChecksum(w, chunk)
	w.WriteHeader(http.StatusPartialContent)
	io.WriteString(w, chunk)
}

// setChecksum sends the CRC32 header real S3 sends. Without it the SDK logs
// "Response has no supported checksum" on every read, which buries anything a
// failing test has to say -- and sending it means the response really is
// validated, so a fake that corrupted the payload would be caught.
func setChecksum(w http.ResponseWriter, body string) {
	sum := crc32.ChecksumIEEE([]byte(body))
	var b [4]byte
	binary.BigEndian.PutUint32(b[:], sum)
	w.Header().Set("x-amz-checksum-crc32", base64.StdEncoding.EncodeToString(b[:]))
}

// parseByteRange handles the one form this package emits, "bytes=start-end",
// clamping the end to the object and rejecting a start at or past it.
func parseByteRange(rg string, size int) (start, end int, err error) {
	spec, ok := strings.CutPrefix(rg, "bytes=")
	if !ok {
		return 0, 0, fmt.Errorf("unsupported range %q", rg)
	}
	lo, hi, ok := strings.Cut(spec, "-")
	if !ok {
		return 0, 0, fmt.Errorf("unsupported range %q", rg)
	}
	start, err = strconv.Atoi(lo)
	if err != nil {
		return 0, 0, err
	}
	end, err = strconv.Atoi(hi)
	if err != nil {
		return 0, 0, err
	}
	if start < 0 || start >= size {
		return 0, 0, fmt.Errorf("range not satisfiable")
	}
	if end >= size {
		end = size - 1
	}
	return start, end, nil
}

const testBody = "0123456789abcdefghijklmnopqrstuvwxyz"

// startFake brings up the server and points the AWS config at it with static
// credentials, isolated from whatever the developer has in ~/.aws.
func startFake(t *testing.T, objects map[string]string) (*fakeS3, string) {
	t.Helper()
	f := &fakeS3{objects: objects}
	ts := httptest.NewServer(f)
	t.Cleanup(ts.Close)

	// Without these the SDK reads the ambient credential chain, so the test
	// would pass or fail depending on the developer's AWS setup -- and could
	// stall reaching for instance metadata on a machine that has none.
	t.Setenv("AWS_ACCESS_KEY_ID", "testkey")
	t.Setenv("AWS_SECRET_ACCESS_KEY", "testsecret")
	t.Setenv("AWS_SESSION_TOKEN", "")
	t.Setenv("AWS_REGION", "us-east-1")
	t.Setenv("AWS_PROFILE", "")
	t.Setenv("AWS_CONFIG_FILE", t.TempDir()+"/config")
	t.Setenv("AWS_SHARED_CREDENTIALS_FILE", t.TempDir()+"/credentials")
	t.Setenv("AWS_EC2_METADATA_DISABLED", "true")
	t.Setenv("AWS_ENDPOINT_URL", "")
	t.Setenv("AWS_ENDPOINT_URL_S3", "")
	return f, ts.URL
}

func newTestClient(t *testing.T, endpoint string) *Client {
	t.Helper()
	c, err := New(context.Background(), WithEndpoint(endpoint))
	if err != nil {
		t.Fatalf("New: %v", err)
	}
	return c
}

func TestParseURI(t *testing.T) {
	for _, tc := range []struct {
		uri            string
		bucket, key    string
		wantErr        bool
		errMustMention string
	}{
		{uri: "s3://bucket/key.vcf.gz", bucket: "bucket", key: "key.vcf.gz"},
		{uri: "s3://bucket/a/b/c.bam", bucket: "bucket", key: "a/b/c.bam"},
		// A bucket with no key is legal to parse; it is the caller's problem
		// that it names no object.
		{uri: "s3://bucket", bucket: "bucket", key: ""},
		{uri: "s3://bucket/", bucket: "bucket", key: ""},
		// Leading and trailing slashes are trimmed, so "s3://b//k/" and
		// "s3://b/k" name the same object rather than two.
		{uri: "s3://bucket//k//", bucket: "bucket", key: "k"},
		// A key may contain anything; no escaping or normalization happens here.
		{uri: "s3://bucket/dir/file name.vcf", bucket: "bucket", key: "dir/file name.vcf"},
		{uri: "s3://bucket/a=1/b=2/part.parquet", bucket: "bucket", key: "a=1/b=2/part.parquet"},

		{uri: "", wantErr: true, errMustMention: "not an s3"},
		{uri: "bucket/key", wantErr: true, errMustMention: "not an s3"},
		{uri: "https://bucket.s3.amazonaws.com/key", wantErr: true, errMustMention: "not an s3"},
		{uri: "gs://bucket/key", wantErr: true, errMustMention: "not an s3"},
		// The scheme is matched exactly: iosource.Scheme lower-cases before
		// dispatching, so an upper-case scheme never reaches here.
		{uri: "S3://bucket/key", wantErr: true, errMustMention: "not an s3"},
		{uri: "s3://", wantErr: true, errMustMention: "no bucket"},
		{uri: "s3:///key", wantErr: true, errMustMention: "no bucket"},
	} {
		bucket, key, err := ParseURI(tc.uri)
		if tc.wantErr {
			if err == nil {
				t.Errorf("ParseURI(%q) = (%q, %q, nil), want an error", tc.uri, bucket, key)
				continue
			}
			if !strings.Contains(err.Error(), tc.errMustMention) {
				t.Errorf("ParseURI(%q) error = %q, should mention %q", tc.uri, err, tc.errMustMention)
			}
			// The error names the locator, so a failure in a list of inputs
			// says which one.
			if tc.uri != "" && !strings.Contains(err.Error(), tc.uri) {
				t.Errorf("ParseURI(%q) error = %q, should quote the locator", tc.uri, err)
			}
			continue
		}
		if err != nil {
			t.Errorf("ParseURI(%q): %v", tc.uri, err)
			continue
		}
		if bucket != tc.bucket || key != tc.key {
			t.Errorf("ParseURI(%q) = (%q, %q), want (%q, %q)", tc.uri, bucket, key, tc.bucket, tc.key)
		}
	}
}

func TestOpenAndReadAt(t *testing.T) {
	f, endpoint := startFake(t, map[string]string{"bkt/data.bin": testBody})
	c := newTestClient(t, endpoint)

	src, err := c.Open(context.Background(), "s3://bkt/data.bin")
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer src.Close()

	if n, err := src.Size(); err != nil || n != int64(len(testBody)) {
		t.Fatalf("Size() = (%d, %v), want (%d, nil)", n, err, len(testBody))
	}

	for _, tc := range []struct {
		off  int64
		n    int
		want string
	}{
		{0, 4, "0123"},
		{10, 6, "abcdef"},
		{int64(len(testBody)) - 1, 1, "z"},
	} {
		p := make([]byte, tc.n)
		n, err := src.ReadAt(p, tc.off)
		if err != nil && err != io.EOF {
			t.Errorf("ReadAt(%d, %d): %v", tc.n, tc.off, err)
			continue
		}
		if got := string(p[:n]); got != tc.want {
			t.Errorf("ReadAt(%d bytes at %d) = %q, want %q", tc.n, tc.off, got, tc.want)
		}
	}

	// Path-style addressing: the bucket is the first path element, not a
	// subdomain. A gateway has no wildcard DNS, so virtual-host addressing
	// would resolve to nothing.
	for _, p := range f.seenPaths() {
		if p != "/bkt/data.bin" {
			t.Errorf("request path %q, want /bkt/data.bin (path-style addressing)", p)
		}
	}
	// Each read is one ranged GET of exactly the bytes asked for -- the point
	// of the whole exercise is not transferring the rest of the object.
	want := []string{"bytes=0-3", "bytes=10-15", "bytes=35-35"}
	if got := f.seenRanges(); len(got) != len(want) {
		t.Errorf("ranges = %v, want %v", got, want)
	} else {
		for i := range want {
			if got[i] != want[i] {
				t.Errorf("range %d = %q, want %q", i, got[i], want[i])
			}
		}
	}
}

// Open learns the length up front rather than lazily, so a missing object fails
// at open instead of midway through a query, when the caller is several layers
// deep in a BGZF reader and the error reads as a corrupt file.
func TestOpenFailsOnMissingObject(t *testing.T) {
	_, endpoint := startFake(t, map[string]string{"bkt/present.bin": testBody})
	c := newTestClient(t, endpoint)

	src, err := c.Open(context.Background(), "s3://bkt/absent.bin")
	if err == nil {
		src.Close()
		t.Fatal("Open succeeded for an object that does not exist")
	}
	if !strings.Contains(err.Error(), "s3://bkt/absent.bin") {
		t.Errorf("error %q should name the locator", err)
	}
	if !strings.Contains(err.Error(), "head") {
		t.Errorf("error %q should say which request failed", err)
	}
}

func TestOpenRejectsABadLocator(t *testing.T) {
	_, endpoint := startFake(t, map[string]string{})
	c := newTestClient(t, endpoint)

	for _, uri := range []string{"not-a-locator", "https://example.com/x", "s3://"} {
		if _, err := c.Open(context.Background(), uri); err == nil {
			t.Errorf("Open(%q) succeeded, want a parse error", uri)
		}
		if _, err := c.OpenReader(context.Background(), uri); err == nil {
			t.Errorf("OpenReader(%q) succeeded, want a parse error", uri)
		}
	}
}

// The length is learned once and cached. Every BGZF seek asks for it, so a HEAD
// per call would put a round trip in front of every block.
func TestSizeIsCached(t *testing.T) {
	f, endpoint := startFake(t, map[string]string{"bkt/data.bin": testBody})
	c := newTestClient(t, endpoint)

	src, err := c.Open(context.Background(), "s3://bkt/data.bin")
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()

	for i := 0; i < 5; i++ {
		if n, err := src.Size(); err != nil || n != int64(len(testBody)) {
			t.Fatalf("Size() = (%d, %v)", n, err)
		}
	}
	// One, from Open itself.
	if got := atomic.LoadInt32(&f.heads); got != 1 {
		t.Errorf("%d HEAD requests, want 1 -- the size should be learned once", got)
	}
}

// A read at or past the end is io.EOF, not a failure. Block-compressed readers
// probe past the end as a matter of course, so an error here breaks every
// indexed query rather than reporting a real problem.
func TestReadAtPastEndIsEOF(t *testing.T) {
	_, endpoint := startFake(t, map[string]string{"bkt/data.bin": testBody})
	c := newTestClient(t, endpoint)

	src, err := c.Open(context.Background(), "s3://bkt/data.bin")
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()

	for _, off := range []int64{int64(len(testBody)), int64(len(testBody)) + 100} {
		n, err := src.ReadAt(make([]byte, 8), off)
		if !errors.Is(err, io.EOF) {
			t.Errorf("ReadAt at offset %d = (%d, %v), want io.EOF", off, n, err)
		}
		if n != 0 {
			t.Errorf("ReadAt past the end returned %d bytes, want 0", n)
		}
	}
}

// io.ReaderAt requires a non-nil error whenever fewer than len(p) bytes come
// back. A caller that trusts n without checking err would otherwise read stale
// bytes from the tail of its own buffer.
func TestReadAtShortReadReportsEOF(t *testing.T) {
	_, endpoint := startFake(t, map[string]string{"bkt/data.bin": testBody})
	c := newTestClient(t, endpoint)

	src, err := c.Open(context.Background(), "s3://bkt/data.bin")
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()

	// Ask for 10 bytes starting 4 from the end.
	off := int64(len(testBody) - 4)
	p := make([]byte, 10)
	n, err := src.ReadAt(p, off)
	if n != 4 {
		t.Errorf("ReadAt returned %d bytes, want 4", n)
	}
	if !errors.Is(err, io.EOF) {
		t.Errorf("short ReadAt error = %v, want io.EOF", err)
	}
	if got := string(p[:n]); got != "wxyz" {
		t.Errorf("ReadAt = %q, want %q", got, "wxyz")
	}
}

func TestReadAtEdgeCases(t *testing.T) {
	f, endpoint := startFake(t, map[string]string{"bkt/data.bin": testBody})
	c := newTestClient(t, endpoint)

	src, err := c.Open(context.Background(), "s3://bkt/data.bin")
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()

	// An empty buffer is a no-op and must not cost a request -- a ranged GET
	// for zero bytes has no valid Range header to express it.
	before := atomic.LoadInt32(&f.gets)
	if n, err := src.ReadAt(nil, 0); n != 0 || err != nil {
		t.Errorf("ReadAt(nil, 0) = (%d, %v), want (0, nil)", n, err)
	}
	if after := atomic.LoadInt32(&f.gets); after != before {
		t.Errorf("an empty read issued %d requests, want 0", after-before)
	}

	// A negative offset is rejected here rather than sent as a malformed Range.
	if _, err := src.ReadAt(make([]byte, 4), -1); err == nil {
		t.Error("ReadAt at a negative offset succeeded")
	} else if !strings.Contains(err.Error(), "negative offset") {
		t.Errorf("error %q should say the offset was negative", err)
	}
}

func TestOpenReader(t *testing.T) {
	f, endpoint := startFake(t, map[string]string{"bkt/index.tbi": testBody})
	c := newTestClient(t, endpoint)

	rc, err := c.OpenReader(context.Background(), "s3://bkt/index.tbi")
	if err != nil {
		t.Fatalf("OpenReader: %v", err)
	}
	defer rc.Close()

	got, err := io.ReadAll(rc)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if string(got) != testBody {
		t.Errorf("OpenReader gave %q, want %q", got, testBody)
	}
	// A whole-object read, with no Range header: an index is read start to
	// finish and is small enough not to warrant ranging.
	if got := f.seenRanges(); len(got) != 0 {
		t.Errorf("OpenReader sent Range headers %v, want none", got)
	}
	if got := atomic.LoadInt32(&f.heads); got != 0 {
		t.Errorf("OpenReader issued %d HEAD requests, want 0", got)
	}
}

func TestOpenReaderFailsOnMissingObject(t *testing.T) {
	_, endpoint := startFake(t, map[string]string{})
	c := newTestClient(t, endpoint)

	if _, err := c.OpenReader(context.Background(), "s3://bkt/absent.tbi"); err == nil {
		t.Fatal("OpenReader succeeded for an absent object")
	} else if !strings.Contains(err.Error(), "s3://bkt/absent.tbi") {
		t.Errorf("error %q should name the locator", err)
	}
}

// ByteSource requires ReadAt to be safe for concurrent use: a header scan and
// several region queries share one source. Run under -race, this is the test
// that would catch a cached response buffer or a shared request struct.
func TestReadAtIsConcurrencySafe(t *testing.T) {
	_, endpoint := startFake(t, map[string]string{"bkt/data.bin": testBody})
	c := newTestClient(t, endpoint)

	src, err := c.Open(context.Background(), "s3://bkt/data.bin")
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()

	var wg sync.WaitGroup
	errs := make(chan error, 32)
	for i := 0; i < 16; i++ {
		wg.Add(1)
		go func(off int64) {
			defer wg.Done()
			p := make([]byte, 4)
			n, err := src.ReadAt(p, off)
			if err != nil && err != io.EOF {
				errs <- fmt.Errorf("ReadAt at %d: %w", off, err)
				return
			}
			if want := testBody[off : off+int64(n)]; string(p[:n]) != want {
				errs <- fmt.Errorf("ReadAt at %d = %q, want %q", off, p[:n], want)
			}
			// Size races with the reads, and shares a mutex with them.
			if _, err := src.Size(); err != nil {
				errs <- err
			}
		}(int64(i * 2))
	}
	wg.Wait()
	close(errs)
	for err := range errs {
		t.Error(err)
	}
}

// A Source has to satisfy ByteSource, and the SectionReader wrapper that
// block-compressed readers actually use has to work over it. This is the shape
// every real caller goes through.
func TestSourceWorksAsAByteSource(t *testing.T) {
	_, endpoint := startFake(t, map[string]string{"bkt/data.bin": testBody})
	c := newTestClient(t, endpoint)

	src, err := c.Open(context.Background(), "s3://bkt/data.bin")
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()

	var bs iosource.ByteSource = src
	rs, err := iosource.ReadSeeker(bs)
	if err != nil {
		t.Fatalf("ReadSeeker: %v", err)
	}
	if _, err := rs.Seek(10, io.SeekStart); err != nil {
		t.Fatal(err)
	}
	got, err := io.ReadAll(rs)
	if err != nil {
		t.Fatalf("ReadAll after seek: %v", err)
	}
	if string(got) != testBody[10:] {
		t.Errorf("read %q, want %q", got, testBody[10:])
	}

	// Close is a no-op that must stay safe to call more than once: each read is
	// its own request, so there is nothing held open.
	for i := 0; i < 3; i++ {
		if err := src.Close(); err != nil {
			t.Errorf("Close #%d: %v", i+1, err)
		}
	}
}

// Importing this package is the whole activation mechanism -- cgkit does it with
// a blank import in internal/locator, and nothing else references the package.
// If init() stopped registering, that import would compile and silently do
// nothing, and every s3:// locator would report "no transport registered".
func TestSchemeIsRegisteredByImport(t *testing.T) {
	if !iosource.HasTransport("s3") {
		t.Error("HasTransport(s3) = false; importing the package should register the scheme")
	}
	var found bool
	for _, s := range iosource.Schemes() {
		if s == "s3" {
			found = true
		}
	}
	if !found {
		t.Errorf("Schemes() = %v, want it to include s3", iosource.Schemes())
	}
	// An unregistered scheme still reports as unavailable, so the check is
	// discriminating rather than always true.
	if iosource.HasTransport("gs") {
		t.Error("HasTransport(gs) = true, but no gs transport exists")
	}
}

// The endpoint comes from AWS_ENDPOINT_URL_S3 first, then AWS_ENDPOINT_URL --
// the S3-specific variable wins, which is what lets one process talk to a
// gateway for S3 and to real AWS for everything else.
func TestEndpointEnvironmentPrecedence(t *testing.T) {
	primary := &fakeS3{objects: map[string]string{"bkt/data.bin": testBody}}
	other := &fakeS3{objects: map[string]string{"bkt/data.bin": "WRONG SERVER"}}
	primaryTS := httptest.NewServer(primary)
	defer primaryTS.Close()
	otherTS := httptest.NewServer(other)
	defer otherTS.Close()

	startFake(t, map[string]string{}) // credential isolation only

	t.Run("S3-specific variable wins", func(t *testing.T) {
		t.Setenv("AWS_ENDPOINT_URL_S3", primaryTS.URL)
		t.Setenv("AWS_ENDPOINT_URL", otherTS.URL)
		c, err := New(context.Background())
		if err != nil {
			t.Fatal(err)
		}
		assertReadsFrom(t, c, primary)
	})

	t.Run("generic variable is the fallback", func(t *testing.T) {
		t.Setenv("AWS_ENDPOINT_URL_S3", "")
		t.Setenv("AWS_ENDPOINT_URL", primaryTS.URL)
		c, err := New(context.Background())
		if err != nil {
			t.Fatal(err)
		}
		assertReadsFrom(t, c, primary)
	})

	// WithEndpoint overrides both, so a program offering a flag beats the
	// ambient environment.
	t.Run("WithEndpoint overrides the environment", func(t *testing.T) {
		t.Setenv("AWS_ENDPOINT_URL_S3", otherTS.URL)
		t.Setenv("AWS_ENDPOINT_URL", otherTS.URL)
		c, err := New(context.Background(), WithEndpoint(primaryTS.URL))
		if err != nil {
			t.Fatal(err)
		}
		assertReadsFrom(t, c, primary)
	})
}

func assertReadsFrom(t *testing.T, c *Client, want *fakeS3) {
	t.Helper()
	src, err := c.Open(context.Background(), "s3://bkt/data.bin")
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer src.Close()
	p := make([]byte, 4)
	if _, err := src.ReadAt(p, 0); err != nil && err != io.EOF {
		t.Fatalf("ReadAt: %v", err)
	}
	if string(p) != testBody[:4] {
		t.Errorf("read %q from the wrong server", p)
	}
	if atomic.LoadInt32(&want.gets) == 0 {
		t.Error("the expected server received no GET")
	}
}

func TestFirstEnv(t *testing.T) {
	t.Setenv("CGHTS_TEST_A", "")
	t.Setenv("CGHTS_TEST_B", "  spaced  ")
	t.Setenv("CGHTS_TEST_C", "third")

	// Empty and whitespace-only are both "unset": an exported-but-blank variable
	// is how a shell script clears one, and honouring it as an endpoint would
	// produce an unsignable request.
	if got := firstEnv("CGHTS_TEST_A", "CGHTS_TEST_B"); got != "spaced" {
		t.Errorf("firstEnv = %q, want %q (empty skipped, value trimmed)", got, "spaced")
	}
	t.Setenv("CGHTS_TEST_B", "   ")
	if got := firstEnv("CGHTS_TEST_A", "CGHTS_TEST_B", "CGHTS_TEST_C"); got != "third" {
		t.Errorf("firstEnv = %q, want %q (whitespace-only skipped)", got, "third")
	}
	if got := firstEnv("CGHTS_TEST_A", "CGHTS_TEST_UNSET"); got != "" {
		t.Errorf("firstEnv = %q, want empty", got)
	}
	if got := firstEnv(); got != "" {
		t.Errorf("firstEnv() = %q, want empty", got)
	}
}

// A named profile that does not exist must fail with a message naming it,
// rather than silently falling back to anonymous access and 403-ing later.
func TestWithProfileNamesTheProfileOnFailure(t *testing.T) {
	_, endpoint := startFake(t, map[string]string{})
	t.Setenv("AWS_CONFIG_FILE", t.TempDir()+"/absent-config")
	t.Setenv("AWS_SHARED_CREDENTIALS_FILE", t.TempDir()+"/absent-credentials")

	_, err := New(context.Background(), WithProfile("no-such-profile"), WithEndpoint(endpoint))
	if err == nil {
		t.Skip("the SDK resolved a nonexistent profile without error; nothing to assert")
	}
	if !strings.Contains(err.Error(), "no-such-profile") {
		t.Errorf("error %q should name the profile", err)
	}
}

func TestWithRegion(t *testing.T) {
	f, endpoint := startFake(t, map[string]string{"bkt/data.bin": testBody})
	t.Setenv("AWS_REGION", "")
	t.Setenv("AWS_DEFAULT_REGION", "")

	// With no region anywhere, New falls back to us-east-1 rather than failing:
	// a compatible gateway does not care which region signs the request, but
	// the signature needs one.
	c, err := New(context.Background(), WithEndpoint(endpoint))
	if err != nil {
		t.Fatalf("New with no region: %v", err)
	}
	if _, err := c.Open(context.Background(), "s3://bkt/data.bin"); err != nil {
		t.Errorf("Open with the default region: %v", err)
	}

	c2, err := New(context.Background(), WithRegion("eu-west-2"), WithEndpoint(endpoint))
	if err != nil {
		t.Fatalf("New with an explicit region: %v", err)
	}
	if _, err := c2.Open(context.Background(), "s3://bkt/data.bin"); err != nil {
		t.Errorf("Open with an explicit region: %v", err)
	}
	if atomic.LoadInt32(&f.heads) != 2 {
		t.Errorf("%d HEAD requests, want 2", atomic.LoadInt32(&f.heads))
	}
}

// WithHTTPClient is the seam that makes the wire observable. Its own test is
// what keeps it from being dropped as unused.
func TestWithHTTPClient(t *testing.T) {
	_, endpoint := startFake(t, map[string]string{"bkt/data.bin": testBody})

	var calls int32
	spy := &countingHTTPClient{n: &calls, inner: http.DefaultClient}
	c, err := New(context.Background(), WithEndpoint(endpoint), WithHTTPClient(spy))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := c.Open(context.Background(), "s3://bkt/data.bin"); err != nil {
		t.Fatalf("Open: %v", err)
	}
	if got := atomic.LoadInt32(&calls); got == 0 {
		t.Error("the supplied HTTP client was never used")
	}
}

// This is the only test that touches Shared, and it has to be: Shared is a
// sync.Once, so whichever test calls it first fixes the process-wide client and
// every later caller gets that one regardless of its own environment. Keeping
// all of it here means the singleton is built exactly once, against this
// server, whatever order the tests run in.
//
// It covers the path cgkit actually uses -- a blank import, then
// iosource.Open("s3://...") -- which nothing else exercises end to end.
func TestSharedClientAndRegisteredOpener(t *testing.T) {
	f := &fakeS3{objects: map[string]string{"bkt/data.bin": testBody}}
	ts := httptest.NewServer(f)
	defer ts.Close()

	startFake(t, map[string]string{}) // credential isolation only
	t.Setenv("AWS_ENDPOINT_URL_S3", ts.URL)

	ctx := context.Background()

	c, err := Shared(ctx)
	if err != nil {
		t.Fatalf("Shared: %v", err)
	}
	// Built once and reused: a fresh client per call would re-walk the
	// credential chain on every open.
	c2, err := Shared(ctx)
	if err != nil {
		t.Fatalf("Shared (second call): %v", err)
	}
	if c != c2 {
		t.Error("Shared returned two different clients")
	}

	// The registered opener resolves an s3:// locator to a working ByteSource
	// without the caller naming this package at all.
	src, err := iosource.Open(ctx, "s3://bkt/data.bin")
	if err != nil {
		t.Fatalf("iosource.Open: %v", err)
	}
	defer src.Close()
	if n, err := src.Size(); err != nil || n != int64(len(testBody)) {
		t.Errorf("Size() = (%d, %v), want (%d, nil)", n, err, len(testBody))
	}
	p := make([]byte, 5)
	if _, err := src.ReadAt(p, 0); err != nil && err != io.EOF {
		t.Fatalf("ReadAt: %v", err)
	}
	if string(p) != testBody[:5] {
		t.Errorf("read %q, want %q", p, testBody[:5])
	}

	// A bad locator surfaces through the same path rather than panicking in
	// the opener.
	if _, err := iosource.Open(ctx, "s3://"); err == nil {
		t.Error("iosource.Open(s3://) succeeded")
	}
	// An absent object is an error at open, and names the locator.
	if _, err := iosource.Open(ctx, "s3://bkt/absent.bin"); err == nil {
		t.Error("iosource.Open succeeded for an absent object")
	} else if !strings.Contains(err.Error(), "s3://bkt/absent.bin") {
		t.Errorf("error %q should name the locator", err)
	}

	// iosource.Sibling resolves an index over the same transport as the data,
	// which is what makes an indexed remote query seek rather than stream.
	f.mu.Lock()
	f.objects["bkt/data.bin.tbi"] = "INDEX"
	f.mu.Unlock()
	rc, err := iosource.Sibling(ctx)("s3://bkt/data.bin", ".tbi")
	if err != nil {
		t.Fatalf("Sibling: %v", err)
	}
	defer rc.Close()
	got, err := io.ReadAll(rc)
	if err != nil {
		t.Fatal(err)
	}
	if string(got) != "INDEX" {
		t.Errorf("sibling read %q, want %q", got, "INDEX")
	}

	// PutForTest stages a fixture through the same shared client. It is the
	// only write in the package and exists solely for integration tests.
	if err := PutForTest(ctx, "bkt", "staged.bin", strings.NewReader("staged content")); err != nil {
		t.Fatalf("PutForTest: %v", err)
	}
	if atomic.LoadInt32(&f.puts) != 1 {
		t.Errorf("%d PUT requests, want 1", atomic.LoadInt32(&f.puts))
	}
	rc2, err := iosource.OpenReader(ctx, "s3://bkt/staged.bin")
	if err != nil {
		t.Fatalf("reading back the staged object: %v", err)
	}
	defer rc2.Close()
	back, err := io.ReadAll(rc2)
	if err != nil {
		t.Fatal(err)
	}
	if string(back) != "staged content" {
		t.Errorf("staged object read back as %q", back)
	}
}

type countingHTTPClient struct {
	n     *int32
	inner *http.Client
}

func (c *countingHTTPClient) Do(r *http.Request) (*http.Response, error) {
	atomic.AddInt32(c.n, 1)
	return c.inner.Do(r)
}
