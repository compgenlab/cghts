package htsio

import (
	"bufio"
	"os"
	"strings"
	"testing"

	"github.com/compgen-io/cgltk/htsio/bgzf"
)

// readBGZFLines reads all non-empty lines from a BGZF-compressed file.
func readBGZFLines(t *testing.T, path string) []string {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer f.Close()

	r := bgzf.NewReader(f)
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 0, 1024*1024), 10*1024*1024)

	var lines []string
	for scanner.Scan() {
		line := scanner.Text()
		if line != "" {
			lines = append(lines, line)
		}
	}
	return lines
}

func TestTabixWriterBEDSort(t *testing.T) {
	dir := t.TempDir()
	outPath := dir + "/sorted.bed.gz"

	opts := NewTabixWriterOpts().BED()
	tw := NewTabixWriter(outPath, opts)

	tw.Write("chr1\t500\t600\tfeature3")
	tw.Write("chr1\t100\t200\tfeature1")
	tw.Write("chr1\t300\t400\tfeature2")

	if err := tw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	lines := readBGZFLines(t, outPath)
	if len(lines) != 3 {
		t.Fatalf("expected 3 lines, got %d", len(lines))
	}

	expected := []string{
		"chr1\t100\t200\tfeature1",
		"chr1\t300\t400\tfeature2",
		"chr1\t500\t600\tfeature3",
	}
	for i, line := range lines {
		if line != expected[i] {
			t.Errorf("line[%d]: got %q, want %q", i, line, expected[i])
		}
	}
}

func TestTabixWriterVCFSort(t *testing.T) {
	dir := t.TempDir()
	outPath := dir + "/sorted.vcf.gz"

	opts := NewTabixWriterOpts().VCF()
	tw := NewTabixWriter(outPath, opts)

	tw.WriteHeader("##fileformat=VCFv4.2")
	tw.WriteHeader("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

	tw.Write("chr1\t500\t.\tA\tT\t30\tPASS\t.")
	tw.Write("chr1\t100\t.\tG\tC\t40\tPASS\t.")
	tw.Write("chr2\t200\t.\tT\tA\t50\tPASS\t.")

	if err := tw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	lines := readBGZFLines(t, outPath)
	if len(lines) != 5 {
		t.Fatalf("expected 5 lines, got %d: %v", len(lines), lines)
	}
	if lines[0] != "##fileformat=VCFv4.2" {
		t.Errorf("header 0: %q", lines[0])
	}
	if lines[2] != "chr1\t100\t.\tG\tC\t40\tPASS\t." {
		t.Errorf("data 0: %q", lines[2])
	}
	if lines[3] != "chr1\t500\t.\tA\tT\t30\tPASS\t." {
		t.Errorf("data 1: %q", lines[3])
	}
	if lines[4] != "chr2\t200\t.\tT\tA\t50\tPASS\t." {
		t.Errorf("data 2: %q", lines[4])
	}
}

func TestTabixWriterMultiRef(t *testing.T) {
	dir := t.TempDir()
	outPath := dir + "/multi.bed.gz"

	opts := NewTabixWriterOpts().BED()
	tw := NewTabixWriter(outPath, opts)

	tw.Write("chr2\t100\t200\tb1")
	tw.Write("chr1\t500\t600\ta2")
	tw.Write("chr1\t100\t200\ta1")
	tw.Write("chr2\t500\t600\tb2")

	if err := tw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	lines := readBGZFLines(t, outPath)
	if len(lines) != 4 {
		t.Fatalf("expected 4 lines, got %d", len(lines))
	}

	// Reference order preserved from first appearance: chr2, chr1.
	expected := []string{
		"chr2\t100\t200\tb1",
		"chr2\t500\t600\tb2",
		"chr1\t100\t200\ta1",
		"chr1\t500\t600\ta2",
	}
	for i, line := range lines {
		if line != expected[i] {
			t.Errorf("line[%d]: got %q, want %q", i, line, expected[i])
		}
	}
}

func TestTabixWriterMerge(t *testing.T) {
	dir := t.TempDir()
	outPath := dir + "/merged.bed.gz"

	opts := NewTabixWriterOpts().BED()
	tw := NewTabixWriter(outPath, opts)
	tw.maxMem = 200 // force multiple temp files

	tw.Write("chr1\t400\t500\tfeature4")
	tw.Write("chr1\t100\t200\tfeature1")
	tw.Write("chr1\t300\t400\tfeature3")
	tw.Write("chr1\t200\t300\tfeature2")

	if err := tw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	// Verify no temp files remain.
	matches, _ := os.ReadDir(dir)
	for _, m := range matches {
		if m.Name() != "merged.bed.gz" {
			t.Errorf("leftover temp file: %s", m.Name())
		}
	}

	lines := readBGZFLines(t, outPath)
	expected := []string{
		"chr1\t100\t200\tfeature1",
		"chr1\t200\t300\tfeature2",
		"chr1\t300\t400\tfeature3",
		"chr1\t400\t500\tfeature4",
	}
	if len(lines) != len(expected) {
		t.Fatalf("expected %d lines, got %d: %v", len(expected), len(lines), lines)
	}
	for i, line := range lines {
		if line != expected[i] {
			t.Errorf("line[%d]: got %q, want %q", i, line, expected[i])
		}
	}
}

func TestTabixWriterAutoIndex(t *testing.T) {
	dir := t.TempDir()
	outPath := dir + "/indexed.bed.gz"

	opts := NewTabixWriterOpts().BED().AutoIndex()
	tw := NewTabixWriter(outPath, opts)

	tw.Write("chr1\t100\t200\tfeature1")
	tw.Write("chr1\t500\t600\tfeature2")
	tw.Write("chr2\t100\t200\tfeature3")

	if err := tw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	// Verify .tbi was created.
	tbiPath := outPath + ".tbi"
	if _, err := os.Stat(tbiPath); err != nil {
		t.Fatalf("TBI not created: %v", err)
	}

	// Read back using TabixReader to verify the index works.
	tr, err := NewTabixReader(outPath)
	if err != nil {
		t.Fatalf("NewTabixReader: %v", err)
	}
	defer tr.Close()

	records, err := tr.Query("chr1", 400, 700)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var names []string
	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter: %v", err)
		}
		fields := strings.Split(rec.Line, "\t")
		if len(fields) >= 4 {
			names = append(names, fields[3])
		}
	}

	if len(names) != 1 || names[0] != "feature2" {
		t.Errorf("expected [feature2], got %v", names)
	}
}

func TestTabixWriterEmpty(t *testing.T) {
	dir := t.TempDir()
	outPath := dir + "/empty.bed.gz"

	opts := NewTabixWriterOpts().BED()
	tw := NewTabixWriter(outPath, opts)

	if err := tw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	lines := readBGZFLines(t, outPath)
	if len(lines) != 0 {
		t.Errorf("expected 0 lines, got %d", len(lines))
	}
}
