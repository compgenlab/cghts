package htsio

import (
	"io"
	"os"
	"testing"
)

func TestSortedBamWriterCoordSort(t *testing.T) {
	dir := t.TempDir()
	outPath := dir + "/sorted.bam"

	header := NewSamHeader()
	header.AddLine("@HD\tVN:1.6\tSO:coordinate")
	header.AddLine("@SQ\tSN:chr1\tLN:100000")
	header.AddLine("@SQ\tSN:chr2\tLN:50000")

	// Write records out of order.
	records := []*SamRecord{
		{ReadName: "read3", Flag: 0, RefName: "chr1", Pos: 1000, MapQ: 60, Cigar: "50M", RefNext: "*", PosNext: 0, InsertLen: 0, Seq: "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC", Qual: "*", Tags: map[string]SamTag{}},
		{ReadName: "read1", Flag: 0, RefName: "chr1", Pos: 100, MapQ: 60, Cigar: "50M", RefNext: "*", PosNext: 0, InsertLen: 0, Seq: "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC", Qual: "*", Tags: map[string]SamTag{}},
		{ReadName: "read5", Flag: 0, RefName: "chr2", Pos: 200, MapQ: 60, Cigar: "50M", RefNext: "*", PosNext: 0, InsertLen: 0, Seq: "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC", Qual: "*", Tags: map[string]SamTag{}},
		{ReadName: "read2", Flag: 0, RefName: "chr1", Pos: 500, MapQ: 60, Cigar: "50M", RefNext: "*", PosNext: 0, InsertLen: 0, Seq: "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC", Qual: "*", Tags: map[string]SamTag{}},
		{ReadName: "unmapped", Flag: 4, RefName: "*", Pos: 0, MapQ: 0, Cigar: "*", RefNext: "*", PosNext: 0, InsertLen: 0, Seq: "ACGT", Qual: "*", Tags: map[string]SamTag{}},
	}

	opts := SamWriterOptions(header).BAM().SortCoord()
	writer, err := NewSamWriter(outPath, opts)
	if err != nil {
		t.Fatalf("NewSamWriter: %v", err)
	}
	for _, rec := range records {
		if err := writer.Write(rec); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := writer.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	// Read back and verify order.
	f, err := os.Open(outPath)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	reader, err := NewBamReader(f)
	if err != nil {
		t.Fatalf("NewBamReader: %v", err)
	}
	defer reader.Close()

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

	// Expected order: chr1:100, chr1:500, chr1:1000, chr2:200, unmapped
	expected := []string{"read1", "read2", "read3", "read5", "unmapped"}
	if len(names) != len(expected) {
		t.Fatalf("expected %d records, got %d: %v", len(expected), len(names), names)
	}
	for i, name := range names {
		if name != expected[i] {
			t.Errorf("record[%d]: got %q, want %q", i, name, expected[i])
		}
	}
}

func TestSortedBamWriterNameSort(t *testing.T) {
	dir := t.TempDir()
	outPath := dir + "/namesorted.bam"

	header := NewSamHeader()
	header.AddLine("@HD\tVN:1.6\tSO:queryname")
	header.AddLine("@SQ\tSN:chr1\tLN:100000")

	records := []*SamRecord{
		{ReadName: "charlie", Flag: 0, RefName: "chr1", Pos: 300, MapQ: 60, Cigar: "50M", RefNext: "*", PosNext: 0, InsertLen: 0, Seq: "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC", Qual: "*", Tags: map[string]SamTag{}},
		{ReadName: "alpha", Flag: 0, RefName: "chr1", Pos: 100, MapQ: 60, Cigar: "50M", RefNext: "*", PosNext: 0, InsertLen: 0, Seq: "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC", Qual: "*", Tags: map[string]SamTag{}},
		{ReadName: "bravo", Flag: 0, RefName: "chr1", Pos: 200, MapQ: 60, Cigar: "50M", RefNext: "*", PosNext: 0, InsertLen: 0, Seq: "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC", Qual: "*", Tags: map[string]SamTag{}},
	}

	opts := SamWriterOptions(header).BAM().SortName()
	writer, err := NewSamWriter(outPath, opts)
	if err != nil {
		t.Fatalf("NewSamWriter: %v", err)
	}
	for _, rec := range records {
		if err := writer.Write(rec); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := writer.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	f, err := os.Open(outPath)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	reader, err := NewBamReader(f)
	if err != nil {
		t.Fatalf("NewBamReader: %v", err)
	}
	defer reader.Close()

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

	expected := []string{"alpha", "bravo", "charlie"}
	if len(names) != len(expected) {
		t.Fatalf("expected %d records, got %d: %v", len(expected), len(names), names)
	}
	for i, name := range names {
		if name != expected[i] {
			t.Errorf("record[%d]: got %q, want %q", i, name, expected[i])
		}
	}
}

func TestSortedBamWriterMerge(t *testing.T) {
	dir := t.TempDir()
	outPath := dir + "/merged.bam"

	header := NewSamHeader()
	header.AddLine("@HD\tVN:1.6\tSO:coordinate")
	header.AddLine("@SQ\tSN:chr1\tLN:100000")

	// Use a tiny buffer to force multiple temp files.
	sw, err := newSortedBamWriter(outPath, header, true, "")
	if err != nil {
		t.Fatalf("newSortedBamWriter: %v", err)
	}
	sw.maxMem = 500 // force flush after ~2 records

	records := []*SamRecord{
		{ReadName: "read4", Flag: 0, RefName: "chr1", Pos: 4000, MapQ: 60, Cigar: "50M", RefNext: "*", PosNext: 0, InsertLen: 0, Seq: "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC", Qual: "*", Tags: map[string]SamTag{}},
		{ReadName: "read1", Flag: 0, RefName: "chr1", Pos: 100, MapQ: 60, Cigar: "50M", RefNext: "*", PosNext: 0, InsertLen: 0, Seq: "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC", Qual: "*", Tags: map[string]SamTag{}},
		{ReadName: "read3", Flag: 0, RefName: "chr1", Pos: 3000, MapQ: 60, Cigar: "50M", RefNext: "*", PosNext: 0, InsertLen: 0, Seq: "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC", Qual: "*", Tags: map[string]SamTag{}},
		{ReadName: "read2", Flag: 0, RefName: "chr1", Pos: 2000, MapQ: 60, Cigar: "50M", RefNext: "*", PosNext: 0, InsertLen: 0, Seq: "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC", Qual: "*", Tags: map[string]SamTag{}},
	}

	for _, rec := range records {
		if err := sw.Write(rec); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := sw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	// Verify no temp files remain.
	matches, _ := os.ReadDir(dir)
	for _, m := range matches {
		if m.Name() != "merged.bam" {
			t.Errorf("leftover temp file: %s", m.Name())
		}
	}

	// Read back and verify order.
	f, err := os.Open(outPath)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	reader, err := NewBamReader(f)
	if err != nil {
		t.Fatalf("NewBamReader: %v", err)
	}
	defer reader.Close()

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

	expected := []string{"read1", "read2", "read3", "read4"}
	if len(names) != len(expected) {
		t.Fatalf("expected %d records, got %d: %v", len(expected), len(names), names)
	}
	for i, name := range names {
		if name != expected[i] {
			t.Errorf("record[%d]: got %q, want %q", i, name, expected[i])
		}
	}
}

func TestSortedBamWriterEmpty(t *testing.T) {
	dir := t.TempDir()
	outPath := dir + "/empty.bam"

	header := NewSamHeader()
	header.AddLine("@HD\tVN:1.6\tSO:coordinate")
	header.AddLine("@SQ\tSN:chr1\tLN:100000")

	opts := SamWriterOptions(header).BAM().SortCoord()
	writer, err := NewSamWriter(outPath, opts)
	if err != nil {
		t.Fatalf("NewSamWriter: %v", err)
	}
	if err := writer.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	// Should produce a valid empty BAM.
	f, err := os.Open(outPath)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	reader, err := NewBamReader(f)
	if err != nil {
		t.Fatalf("NewBamReader: %v", err)
	}
	defer reader.Close()

	_, err = reader.Next()
	if err != io.EOF {
		t.Errorf("expected io.EOF, got %v", err)
	}
}
