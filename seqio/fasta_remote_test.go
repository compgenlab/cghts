package seqio

import (
	"fmt"
	"net"
	"net/http"
	"os"
	"path/filepath"
	"testing"
)

// startTestServer serves a local directory over HTTP for testing.
func startTestServer(t *testing.T, dir string) string {
	t.Helper()
	listener, err := net.Listen("tcp", "127.0.0.1:0")
	if err != nil {
		t.Fatal(err)
	}
	addr := listener.Addr().String()

	srv := &http.Server{Handler: http.FileServer(http.Dir(dir))}
	go srv.Serve(listener)
	t.Cleanup(func() { srv.Close() })

	return "http://" + addr
}

func TestRemoteFastaReader(t *testing.T) {
	// Create a test FASTA + .fai in a temp directory and serve it.
	dir := t.TempDir()
	fastaPath := createTestFasta(t, dir)

	baseURL := startTestServer(t, dir)
	fastaURL := baseURL + "/" + filepath.Base(fastaPath)

	r, err := NewRemoteFastaReader(fastaURL)
	if err != nil {
		t.Fatalf("NewRemoteFastaReader: %v", err)
	}
	defer r.Close()

	// Verify names.
	names := r.Names()
	if len(names) != 2 || names[0] != "seq1" || names[1] != "seq2" {
		t.Errorf("Names() = %v, want [seq1, seq2]", names)
	}

	// Verify lengths.
	if l, ok := r.SequenceLength("seq1"); !ok || l != 2000 {
		t.Errorf("SequenceLength(seq1) = %d, %v", l, ok)
	}

	// Get a range.
	seq, err := r.GetSequenceRange("seq1", 0, 10)
	if err != nil {
		t.Fatalf("GetSequenceRange: %v", err)
	}
	if string(seq) != "ACGTACGTNN" {
		t.Errorf("got %q, want ACGTACGTNN", string(seq))
	}

	// Cross-validate: full sequence from remote matches local.
	localR, err := NewIndexedFastaReader(fastaPath)
	if err != nil {
		t.Fatalf("local reader: %v", err)
	}
	defer localR.Close()

	localSeq, _ := localR.GetSequence("seq1")
	remoteSeq, err := r.GetSequence("seq1")
	if err != nil {
		t.Fatalf("remote GetSequence: %v", err)
	}
	if string(remoteSeq) != string(localSeq) {
		t.Errorf("remote/local mismatch: remote=%d bytes, local=%d bytes", len(remoteSeq), len(localSeq))
	}

	// Range in the middle.
	partial, err := r.GetSequenceRange("seq1", 500, 600)
	if err != nil {
		t.Fatalf("GetSequenceRange(500,600): %v", err)
	}
	localPartial, _ := localR.GetSequenceRange("seq1", 500, 600)
	if string(partial) != string(localPartial) {
		t.Errorf("partial mismatch")
	}
}

func TestOpenReferenceURL(t *testing.T) {
	dir := t.TempDir()
	fastaPath := createTestFasta(t, dir)

	baseURL := startTestServer(t, dir)
	fastaURL := baseURL + "/" + filepath.Base(fastaPath)

	// OpenReference should detect HTTP URL and use RemoteFastaReader.
	ref, err := OpenReference(fastaURL)
	if err != nil {
		t.Fatalf("OpenReference(%s): %v", fastaURL, err)
	}
	defer ref.Close()

	if _, ok := ref.(*RemoteFastaReader); !ok {
		t.Errorf("expected RemoteFastaReader, got %T", ref)
	}

	seq, err := ref.GetSequenceRange("seq1", 0, 10)
	if err != nil {
		t.Fatalf("GetSequenceRange: %v", err)
	}
	if string(seq) != "ACGTACGTNN" {
		t.Errorf("got %q", string(seq))
	}
}

func TestRemoteFastaReaderWithCRAMTestdata(t *testing.T) {
	// Serve the CRAM testdata directory.
	refDir := "../htsio/cram/testdata"
	fastaPath := filepath.Join(refDir, "ref.fa")
	if _, err := os.Stat(fastaPath); err != nil {
		t.Skip("CRAM testdata not available")
	}

	baseURL := startTestServer(t, refDir)
	fastaURL := baseURL + "/ref.fa"

	r, err := NewRemoteFastaReader(fastaURL)
	if err != nil {
		t.Fatalf("NewRemoteFastaReader: %v", err)
	}
	defer r.Close()

	l, ok := r.SequenceLength("chr1")
	if !ok || l != 100000 {
		t.Errorf("chr1 length = %d, %v", l, ok)
	}

	// Fetch a range and validate against local.
	local, err := NewIndexedFastaReader(fastaPath)
	if err != nil {
		t.Fatalf("local: %v", err)
	}
	defer local.Close()

	for _, tc := range []struct{ start, end int }{
		{0, 100},
		{5000, 6000},
		{99900, 100000},
	} {
		t.Run(fmt.Sprintf("%d-%d", tc.start, tc.end), func(t *testing.T) {
			remote, err := r.GetSequenceRange("chr1", tc.start, tc.end)
			if err != nil {
				t.Fatalf("remote: %v", err)
			}
			lseq, _ := local.GetSequenceRange("chr1", tc.start, tc.end)
			if string(remote) != string(lseq) {
				t.Errorf("mismatch at [%d,%d): remote=%d local=%d", tc.start, tc.end, len(remote), len(lseq))
			}
		})
	}
}
