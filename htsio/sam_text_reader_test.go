package htsio

import (
	"io"
	"testing"
)

func TestSamTextReaderBasic(t *testing.T) {
	reader, err := NewSamTextReader("testdata/test.sam")
	if err != nil {
		t.Fatalf("NewSamTextReader: %v", err)
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
	if len(refs) != 2 || refs[0].Name != "chr1" || refs[1].Name != "chr2" {
		t.Errorf("refs: %v", refs)
	}

	var names []string
	for {
		rec, err := reader.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Next: %v", err)
		}
		names = append(names, rec.ReadName)
	}

	expected := []string{"read1", "read2", "read3", "read4", "read5"}
	if len(names) != len(expected) {
		t.Fatalf("expected %d records, got %d: %v", len(expected), len(names), names)
	}
	for i, name := range names {
		if name != expected[i] {
			t.Errorf("record[%d]: got %q, want %q", i, name, expected[i])
		}
	}
}

func TestSamTextReaderViaNewSamReader(t *testing.T) {
	reader, err := NewSamReader("testdata/test.sam")
	if err != nil {
		t.Fatalf("NewSamReader: %v", err)
	}
	defer reader.Close()

	count := 0
	for {
		_, err := reader.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Next: %v", err)
		}
		count++
	}

	if count != 5 {
		t.Errorf("expected 5 records, got %d", count)
	}
}

func TestSamTextReaderFlagFilter(t *testing.T) {
	opts := NewSamReaderOpts().FlagFilter(0x100) // exclude secondary
	reader, err := NewSamTextReader("testdata/test.sam", opts)
	if err != nil {
		t.Fatalf("NewSamTextReader: %v", err)
	}
	defer reader.Close()

	count := 0
	for {
		_, err := reader.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Next: %v", err)
		}
		count++
	}

	// All 5 records are primary, so all should pass.
	if count != 5 {
		t.Errorf("expected 5 records, got %d", count)
	}
}

func TestSamTextReaderQueryNotSupported(t *testing.T) {
	reader, err := NewSamTextReader("testdata/test.sam")
	if err != nil {
		t.Fatalf("NewSamTextReader: %v", err)
	}
	defer reader.Close()

	_, err = reader.Query("chr1", 0, 100)
	if err == nil {
		t.Error("expected error for Query on SAM text reader")
	}
}
