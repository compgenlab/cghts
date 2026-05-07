package htsio_test

import (
	"testing"

	"github.com/compgen-io/cgltk/htsio"
	_ "github.com/compgen-io/cgltk/htsio/bam"
	_ "github.com/compgen-io/cgltk/htsio/cram"
	_ "github.com/compgen-io/cgltk/htsio/sam"
)

const testBAMPath = "testdata/test.bam"

func TestBamQuerySingleResult(t *testing.T) {
	reader, err := htsio.NewSamReader(testBAMPath)
	if err != nil {
		t.Fatalf("NewSamReader: %v", err)
	}
	defer reader.Close()

	records, err := reader.Query("chr1", 400, 600)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var names []string
	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter: %v", err)
		}
		names = append(names, rec.ReadName)
	}

	if len(names) != 1 || names[0] != "read2" {
		t.Errorf("expected [read2], got %v", names)
	}
}

func TestBamQueryMultipleResults(t *testing.T) {
	reader, err := htsio.NewSamReader(testBAMPath)
	if err != nil {
		t.Fatalf("NewSamReader: %v", err)
	}
	defer reader.Close()

	records, err := reader.Query("chr1", 0, 300)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var names []string
	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter: %v", err)
		}
		names = append(names, rec.ReadName)
	}

	if len(names) != 1 || names[0] != "read1" {
		t.Errorf("expected [read1], got %v", names)
	}
}

func TestBamQueryNoResults(t *testing.T) {
	reader, err := htsio.NewSamReader(testBAMPath)
	if err != nil {
		t.Fatalf("NewSamReader: %v", err)
	}
	defer reader.Close()

	records, err := reader.Query("chr1", 2000, 3000)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	count := 0
	for _, err := range records {
		if err != nil {
			t.Fatalf("iter: %v", err)
		}
		count++
	}
	if count != 0 {
		t.Errorf("expected 0 records, got %d", count)
	}
}

func TestBamQueryUnknownRef(t *testing.T) {
	reader, err := htsio.NewSamReader(testBAMPath)
	if err != nil {
		t.Fatalf("NewSamReader: %v", err)
	}
	defer reader.Close()

	_, err = reader.Query("chrX", 0, 100)
	if err == nil {
		t.Error("expected error for unknown reference")
	}
}

func TestBamQueryMultipleQueries(t *testing.T) {
	reader, err := htsio.NewSamReader(testBAMPath)
	if err != nil {
		t.Fatalf("NewSamReader: %v", err)
	}
	defer reader.Close()

	// First query: get read1.
	records, err := reader.Query("chr1", 0, 200)
	if err != nil {
		t.Fatalf("Query 1: %v", err)
	}
	var name string
	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter 1: %v", err)
		}
		name = rec.ReadName
		break
	}
	if name != "read1" {
		t.Errorf("query 1: got %q, want read1", name)
	}

	// Second query: get read4 (tests re-seeking with cache).
	records, err = reader.Query("chr1", 4900, 5100)
	if err != nil {
		t.Fatalf("Query 2: %v", err)
	}
	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter 2: %v", err)
		}
		name = rec.ReadName
		break
	}
	if name != "read4" {
		t.Errorf("query 2: got %q, want read4", name)
	}
}

func TestBamQueryChr2(t *testing.T) {
	reader, err := htsio.NewSamReader(testBAMPath)
	if err != nil {
		t.Fatalf("NewSamReader: %v", err)
	}
	defer reader.Close()

	records, err := reader.Query("chr2", 0, 1000)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var names []string
	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter: %v", err)
		}
		names = append(names, rec.ReadName)
	}

	if len(names) != 1 || names[0] != "read5" {
		t.Errorf("expected [read5], got %v", names)
	}
}

func TestBamQueryHeader(t *testing.T) {
	reader, err := htsio.NewSamReader(testBAMPath)
	if err != nil {
		t.Fatalf("NewSamReader: %v", err)
	}
	defer reader.Close()

	hdr, err := reader.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if hdr == nil {
		t.Fatal("expected non-nil header")
	}

	refs := hdr.References()
	if len(refs) != 2 {
		t.Fatalf("expected 2 refs, got %d", len(refs))
	}
	if refs[0].Name != "chr1" || refs[1].Name != "chr2" {
		t.Errorf("refs: %v", refs)
	}
}
