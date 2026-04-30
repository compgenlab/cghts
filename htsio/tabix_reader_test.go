package htsio

import (
	"io"
	"testing"
)

const testBEDPath = "testdata/test.bed.gz"
const testVCFPath = "testdata/test.vcf.gz"
const testCSIBEDPath = "testdata/test_csi.bed.gz"
const testCSIVCFPath = "testdata/test_csi.vcf.gz"

func TestTabixReaderBEDQuery(t *testing.T) {
	tr, err := NewTabixReader(testBEDPath)
	if err != nil {
		t.Fatalf("NewTabixReader: %v", err)
	}
	defer tr.Close()

	// Query [400, 600) — should match feature2 (400-600).
	iter, err := tr.Query("chr1", 400, 600)
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
		// BED: 4th field is name
		fields := splitTabs(rec.Line)
		if len(fields) >= 4 {
			names = append(names, fields[3])
		}
	}

	if len(names) != 1 || names[0] != "feature2" {
		t.Errorf("expected [feature2], got %v", names)
	}
}

func TestTabixReaderBEDQueryMultiple(t *testing.T) {
	tr, err := NewTabixReader(testBEDPath)
	if err != nil {
		t.Fatalf("NewTabixReader: %v", err)
	}
	defer tr.Close()

	// Query [0, 700) — should match feature1 (50-150) and feature2 (400-600).
	iter, err := tr.Query("chr1", 0, 700)
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
		fields := splitTabs(rec.Line)
		if len(fields) >= 4 {
			names = append(names, fields[3])
		}
	}

	if len(names) != 2 || names[0] != "feature1" || names[1] != "feature2" {
		t.Errorf("expected [feature1, feature2], got %v", names)
	}
}

func TestTabixReaderBEDQueryNoResults(t *testing.T) {
	tr, err := NewTabixReader(testBEDPath)
	if err != nil {
		t.Fatalf("NewTabixReader: %v", err)
	}
	defer tr.Close()

	// Query a gap between features.
	iter, err := tr.Query("chr1", 2000, 3000)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	_, err = iter.Next()
	if err != io.EOF {
		t.Errorf("expected io.EOF, got %v", err)
	}
}

func TestTabixReaderBEDQueryChr2(t *testing.T) {
	tr, err := NewTabixReader(testBEDPath)
	if err != nil {
		t.Fatalf("NewTabixReader: %v", err)
	}
	defer tr.Close()

	iter, err := tr.Query("chr2", 0, 1000)
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
		fields := splitTabs(rec.Line)
		if len(fields) >= 4 {
			names = append(names, fields[3])
		}
	}

	if len(names) != 1 || names[0] != "feature5" {
		t.Errorf("expected [feature5], got %v", names)
	}
}

func TestTabixReaderVCFQuery(t *testing.T) {
	tr, err := NewTabixReader(testVCFPath)
	if err != nil {
		t.Fatalf("NewTabixReader: %v", err)
	}
	defer tr.Close()

	// VCF is 1-based. Query [400, 600) 0-based should match
	// the variant at VCF pos 500 (0-based 499).
	iter, err := tr.Query("chr1", 400, 600)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var positions []int
	for {
		rec, err := iter.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Next: %v", err)
		}
		// rec.Start is 0-based
		positions = append(positions, rec.Start)
	}

	// VCF pos 500 → 0-based 499
	if len(positions) != 1 || positions[0] != 499 {
		t.Errorf("expected [499], got %v", positions)
	}
}

func TestTabixReaderVCFQueryMultiple(t *testing.T) {
	tr, err := NewTabixReader(testVCFPath)
	if err != nil {
		t.Fatalf("NewTabixReader: %v", err)
	}
	defer tr.Close()

	// Query [0, 600) should match VCF pos 100 and 500 (0-based 99 and 499).
	iter, err := tr.Query("chr1", 0, 600)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var positions []int
	for {
		rec, err := iter.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Next: %v", err)
		}
		positions = append(positions, rec.Start)
	}

	if len(positions) != 2 || positions[0] != 99 || positions[1] != 499 {
		t.Errorf("expected [99, 499], got %v", positions)
	}
}

func TestTabixReaderUnknownRef(t *testing.T) {
	tr, err := NewTabixReader(testBEDPath)
	if err != nil {
		t.Fatalf("NewTabixReader: %v", err)
	}
	defer tr.Close()

	_, err = tr.Query("chrX", 0, 100)
	if err == nil {
		t.Error("expected error for unknown reference")
	}
}

func TestTabixReaderMultipleQueries(t *testing.T) {
	tr, err := NewTabixReader(testBEDPath)
	if err != nil {
		t.Fatalf("NewTabixReader: %v", err)
	}
	defer tr.Close()

	// First query
	iter, err := tr.Query("chr1", 0, 200)
	if err != nil {
		t.Fatalf("Query 1: %v", err)
	}
	rec, err := iter.Next()
	if err != nil {
		t.Fatalf("Next 1: %v", err)
	}
	fields := splitTabs(rec.Line)
	if len(fields) < 4 || fields[3] != "feature1" {
		t.Errorf("query 1: got %q", rec.Line)
	}

	// Second query (tests re-seeking)
	iter, err = tr.Query("chr2", 0, 1000)
	if err != nil {
		t.Fatalf("Query 2: %v", err)
	}
	rec, err = iter.Next()
	if err != nil {
		t.Fatalf("Next 2: %v", err)
	}
	fields = splitTabs(rec.Line)
	if len(fields) < 4 || fields[3] != "feature5" {
		t.Errorf("query 2: got %q", rec.Line)
	}
}

func TestTabixReaderCoordinates(t *testing.T) {
	tr, err := NewTabixReader(testBEDPath)
	if err != nil {
		t.Fatalf("NewTabixReader: %v", err)
	}
	defer tr.Close()

	// BED is 0-based. feature1 is chr1:50-150.
	iter, err := tr.Query("chr1", 50, 150)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	rec, err := iter.Next()
	if err != nil {
		t.Fatalf("Next: %v", err)
	}
	if rec.Start != 50 || rec.End != 150 {
		t.Errorf("coordinates: got [%d, %d), want [50, 150)", rec.Start, rec.End)
	}
}

// CSI index tests — same queries as TBI, using CSI-indexed files.

func TestTabixReaderCSIBEDQuery(t *testing.T) {
	tr, err := NewTabixReader(testCSIBEDPath)
	if err != nil {
		t.Fatalf("NewTabixReader: %v", err)
	}
	defer tr.Close()

	iter, err := tr.Query("chr1", 400, 600)
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
		fields := splitTabs(rec.Line)
		if len(fields) >= 4 {
			names = append(names, fields[3])
		}
	}

	if len(names) != 1 || names[0] != "feature2" {
		t.Errorf("expected [feature2], got %v", names)
	}
}

func TestTabixReaderCSIBEDQueryMultiple(t *testing.T) {
	tr, err := NewTabixReader(testCSIBEDPath)
	if err != nil {
		t.Fatalf("NewTabixReader: %v", err)
	}
	defer tr.Close()

	iter, err := tr.Query("chr1", 0, 700)
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
		fields := splitTabs(rec.Line)
		if len(fields) >= 4 {
			names = append(names, fields[3])
		}
	}

	if len(names) != 2 || names[0] != "feature1" || names[1] != "feature2" {
		t.Errorf("expected [feature1, feature2], got %v", names)
	}
}

func TestTabixReaderCSIVCFQuery(t *testing.T) {
	tr, err := NewTabixReader(testCSIVCFPath)
	if err != nil {
		t.Fatalf("NewTabixReader: %v", err)
	}
	defer tr.Close()

	iter, err := tr.Query("chr1", 400, 600)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var positions []int
	for {
		rec, err := iter.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Next: %v", err)
		}
		positions = append(positions, rec.Start)
	}

	if len(positions) != 1 || positions[0] != 499 {
		t.Errorf("expected [499], got %v", positions)
	}
}

func TestTabixReaderCSIChr2(t *testing.T) {
	tr, err := NewTabixReader(testCSIBEDPath)
	if err != nil {
		t.Fatalf("NewTabixReader: %v", err)
	}
	defer tr.Close()

	iter, err := tr.Query("chr2", 0, 1000)
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
		fields := splitTabs(rec.Line)
		if len(fields) >= 4 {
			names = append(names, fields[3])
		}
	}

	if len(names) != 1 || names[0] != "feature5" {
		t.Errorf("expected [feature5], got %v", names)
	}
}

func splitTabs(s string) []string {
	return splitByByte(s, '\t')
}

func splitByByte(s string, sep byte) []string {
	var parts []string
	start := 0
	for i := 0; i < len(s); i++ {
		if s[i] == sep {
			parts = append(parts, s[start:i])
			start = i + 1
		}
	}
	parts = append(parts, s[start:])
	return parts
}
