package htsio

import (
	"io"
	"testing"
)

const testBAMPath = "testdata/test.bam"

func TestIndexedSamReaderQuery(t *testing.T) {
	isr, err := NewIndexedSamReader(testBAMPath)
	if err != nil {
		t.Fatalf("NewIndexedSamReader: %v", err)
	}
	defer isr.Close()

	// Query [400, 600) — should match read2 (pos 500, 50M → [499, 549) 0-based).
	iter, err := isr.Query("chr1", 400, 600)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var names []string
	for {
		rec, err := iter.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Next: %v", err)
		}
		names = append(names, rec.ReadName)
	}

	if len(names) != 1 || names[0] != "read2" {
		t.Errorf("expected [read2], got %v", names)
	}
}

func TestIndexedSamReaderQueryMultiple(t *testing.T) {
	isr, err := NewIndexedSamReader(testBAMPath)
	if err != nil {
		t.Fatalf("NewIndexedSamReader: %v", err)
	}
	defer isr.Close()

	// Query [0, 300) — should match read1 (pos 100).
	iter, err := isr.Query("chr1", 0, 300)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var names []string
	for {
		rec, err := iter.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Next: %v", err)
		}
		names = append(names, rec.ReadName)
	}

	if len(names) != 1 || names[0] != "read1" {
		t.Errorf("expected [read1], got %v", names)
	}
}

func TestIndexedSamReaderQueryNoResults(t *testing.T) {
	isr, err := NewIndexedSamReader(testBAMPath)
	if err != nil {
		t.Fatalf("NewIndexedSamReader: %v", err)
	}
	defer isr.Close()

	// Query a region between reads.
	iter, err := isr.Query("chr1", 2000, 3000)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	rec, err := iter.Next()
	if err != io.EOF {
		t.Errorf("expected io.EOF, got rec=%v err=%v", rec, err)
	}
}

func TestIndexedSamReaderUnknownRef(t *testing.T) {
	isr, err := NewIndexedSamReader(testBAMPath)
	if err != nil {
		t.Fatalf("NewIndexedSamReader: %v", err)
	}
	defer isr.Close()

	_, err = isr.Query("chrX", 0, 100)
	if err == nil {
		t.Error("expected error for unknown reference")
	}
}

func TestIndexedSamReaderMultipleQueries(t *testing.T) {
	isr, err := NewIndexedSamReader(testBAMPath)
	if err != nil {
		t.Fatalf("NewIndexedSamReader: %v", err)
	}
	defer isr.Close()

	// First query: get read1.
	iter, err := isr.Query("chr1", 0, 200)
	if err != nil {
		t.Fatalf("Query 1: %v", err)
	}
	rec, err := iter.Next()
	if err != nil {
		t.Fatalf("Next 1: %v", err)
	}
	if rec.ReadName != "read1" {
		t.Errorf("query 1: got %q, want read1", rec.ReadName)
	}

	// Second query: get read4 (tests re-seeking with cache).
	iter, err = isr.Query("chr1", 4900, 5100)
	if err != nil {
		t.Fatalf("Query 2: %v", err)
	}
	rec, err = iter.Next()
	if err != nil {
		t.Fatalf("Next 2: %v", err)
	}
	if rec.ReadName != "read4" {
		t.Errorf("query 2: got %q, want read4", rec.ReadName)
	}
}

func TestIndexedSamReaderQueryChr2(t *testing.T) {
	isr, err := NewIndexedSamReader(testBAMPath)
	if err != nil {
		t.Fatalf("NewIndexedSamReader: %v", err)
	}
	defer isr.Close()

	// Query chr2 — should get read5.
	iter, err := isr.Query("chr2", 0, 1000)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var names []string
	for {
		rec, err := iter.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Next: %v", err)
		}
		names = append(names, rec.ReadName)
	}

	if len(names) != 1 || names[0] != "read5" {
		t.Errorf("expected [read5], got %v", names)
	}
}

func TestIndexedSamReaderHeader(t *testing.T) {
	isr, err := NewIndexedSamReader(testBAMPath)
	if err != nil {
		t.Fatalf("NewIndexedSamReader: %v", err)
	}
	defer isr.Close()

	hdr, err := isr.Header()
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
