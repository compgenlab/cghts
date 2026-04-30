package htsio

import (
	"strings"
	"testing"
)

const testBEDPath = "testdata/test.bed.gz"
const testVCFPath = "testdata/test.vcf.gz"
const testCSIBEDPath = "testdata/test_csi.bed.gz"
const testCSIVCFPath = "testdata/test_csi.vcf.gz"

func bedName(line string) string {
	fields := strings.Split(line, "\t")
	if len(fields) >= 4 {
		return fields[3]
	}
	return ""
}

func TestTabixReaderBEDQuery(t *testing.T) {
	tr, err := NewTabixReader(testBEDPath)
	if err != nil {
		t.Fatalf("NewTabixReader: %v", err)
	}
	defer tr.Close()

	records, err := tr.Query("chr1", 400, 600)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var names []string
	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter: %v", err)
		}
		names = append(names, bedName(rec.Line))
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

	records, err := tr.Query("chr1", 0, 700)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var names []string
	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter: %v", err)
		}
		names = append(names, bedName(rec.Line))
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

	records, err := tr.Query("chr1", 2000, 3000)
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

func TestTabixReaderBEDQueryChr2(t *testing.T) {
	tr, err := NewTabixReader(testBEDPath)
	if err != nil {
		t.Fatalf("NewTabixReader: %v", err)
	}
	defer tr.Close()

	records, err := tr.Query("chr2", 0, 1000)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var names []string
	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter: %v", err)
		}
		names = append(names, bedName(rec.Line))
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

	records, err := tr.Query("chr1", 400, 600)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var positions []int
	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter: %v", err)
		}
		positions = append(positions, rec.Start)
	}

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

	records, err := tr.Query("chr1", 0, 600)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var positions []int
	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter: %v", err)
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
	records, err := tr.Query("chr1", 0, 200)
	if err != nil {
		t.Fatalf("Query 1: %v", err)
	}
	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter 1: %v", err)
		}
		if bedName(rec.Line) != "feature1" {
			t.Errorf("query 1: got %q", rec.Line)
		}
		break
	}

	// Second query (tests re-seeking)
	records, err = tr.Query("chr2", 0, 1000)
	if err != nil {
		t.Fatalf("Query 2: %v", err)
	}
	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter 2: %v", err)
		}
		if bedName(rec.Line) != "feature5" {
			t.Errorf("query 2: got %q", rec.Line)
		}
		break
	}
}

func TestTabixReaderCoordinates(t *testing.T) {
	tr, err := NewTabixReader(testBEDPath)
	if err != nil {
		t.Fatalf("NewTabixReader: %v", err)
	}
	defer tr.Close()

	records, err := tr.Query("chr1", 50, 150)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter: %v", err)
		}
		if rec.Start != 50 || rec.End != 150 {
			t.Errorf("coordinates: got [%d, %d), want [50, 150)", rec.Start, rec.End)
		}
		break
	}
}

// CSI index tests

func TestTabixReaderCSIBEDQuery(t *testing.T) {
	tr, err := NewTabixReader(testCSIBEDPath)
	if err != nil {
		t.Fatalf("NewTabixReader: %v", err)
	}
	defer tr.Close()

	records, err := tr.Query("chr1", 400, 600)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var names []string
	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter: %v", err)
		}
		names = append(names, bedName(rec.Line))
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

	records, err := tr.Query("chr1", 0, 700)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var names []string
	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter: %v", err)
		}
		names = append(names, bedName(rec.Line))
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

	records, err := tr.Query("chr1", 400, 600)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var positions []int
	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter: %v", err)
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

	records, err := tr.Query("chr2", 0, 1000)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var names []string
	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter: %v", err)
		}
		names = append(names, bedName(rec.Line))
	}

	if len(names) != 1 || names[0] != "feature5" {
		t.Errorf("expected [feature5], got %v", names)
	}
}
