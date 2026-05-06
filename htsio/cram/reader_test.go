package cram

import (
	"fmt"
	"os/exec"
	"strings"
	"testing"
)

func TestReadCRAMRaw(t *testing.T) {
	cramFile := "testdata/test_raw.cram"
	refFile := "testdata/ref.fa"

	cmd := exec.Command("samtools", "view", "-T", refFile, cramFile)
	out, err := cmd.Output()
	if err != nil {
		t.Fatalf("samtools view failed: %v", err)
	}
	expectedLines := strings.Split(strings.TrimSpace(string(out)), "\n")

	reader, err := NewReader(cramFile, refFile)
	if err != nil {
		t.Fatalf("NewReader failed: %v", err)
	}
	defer reader.Close()

	var gotRecords []string
	for rec, err := range reader.Records() {
		if err != nil {
			t.Fatalf("Records() error: %v", err)
		}
		line := fmt.Sprintf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
			rec.ReadName, rec.Flag, rec.RefName, rec.Pos, rec.MapQ,
			rec.Cigar, rec.RefNext, rec.PosNext, rec.InsertLen, rec.Seq, rec.Qual)
		gotRecords = append(gotRecords, line)
	}

	if len(gotRecords) != len(expectedLines) {
		t.Fatalf("record count mismatch: got %d, expected %d", len(gotRecords), len(expectedLines))
	}

	for i, expected := range expectedLines {
		expFields := strings.SplitN(expected, "\t", 12)
		gotFields := strings.SplitN(gotRecords[i], "\t", 12)
		expCore := strings.Join(expFields[:11], "\t")
		gotCore := strings.Join(gotFields[:11], "\t")
		if expCore != gotCore {
			t.Errorf("record %d mismatch:\n  got:  %s\n  want: %s", i, gotCore, expCore)
		}
	}
}

func TestReadCRAM(t *testing.T) {
	cramFile := "testdata/test.cram"
	refFile := "testdata/ref.fa"

	// Get expected output from samtools.
	cmd := exec.Command("samtools", "view", "-T", refFile, cramFile)
	out, err := cmd.Output()
	if err != nil {
		t.Fatalf("samtools view failed: %v", err)
	}
	expectedLines := strings.Split(strings.TrimSpace(string(out)), "\n")

	// Read with our CRAM reader.
	reader, err := NewReader(cramFile, refFile)
	if err != nil {
		t.Fatalf("NewReader failed: %v", err)
	}
	defer reader.Close()

	hdr, err := reader.Header()
	if err != nil {
		t.Fatalf("Header failed: %v", err)
	}
	if hdr == nil {
		t.Fatal("Header is nil")
	}

	// Check we got @SQ lines.
	refs := hdr.References()
	if len(refs) != 2 {
		t.Fatalf("expected 2 references, got %d", len(refs))
	}
	if refs[0].Name != "chr1" || refs[1].Name != "chr2" {
		t.Fatalf("unexpected references: %v", refs)
	}

	var gotRecords []string
	for rec, err := range reader.Records() {
		if err != nil {
			t.Fatalf("Records() error: %v", err)
		}
		// Format as SAM line (core fields only, no tags for now).
		line := fmt.Sprintf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
			rec.ReadName, rec.Flag, rec.RefName, rec.Pos, rec.MapQ,
			rec.Cigar, rec.RefNext, rec.PosNext, rec.InsertLen, rec.Seq, rec.Qual)
		gotRecords = append(gotRecords, line)
	}

	if len(gotRecords) != len(expectedLines) {
		t.Fatalf("record count mismatch: got %d, expected %d", len(gotRecords), len(expectedLines))
	}

	// Compare core fields (first 11 tab-separated fields).
	for i, expected := range expectedLines {
		expFields := strings.SplitN(expected, "\t", 12)
		gotFields := strings.SplitN(gotRecords[i], "\t", 12)

		expCore := strings.Join(expFields[:11], "\t")
		gotCore := strings.Join(gotFields[:11], "\t")

		if expCore != gotCore {
			t.Errorf("record %d mismatch:\n  got:  %s\n  want: %s", i, gotCore, expCore)
		}
	}
}

func TestReadCRAMSimple(t *testing.T) {
	cramFile := "testdata/simple.cram"
	refFile := "testdata/ref.fa"

	reader, err := NewReader(cramFile, refFile)
	if err != nil {
		t.Fatalf("NewReader failed: %v", err)
	}
	defer reader.Close()

	for rec, err := range reader.Records() {
		if err != nil {
			t.Fatalf("Records() error: %v", err)
		}
		t.Logf("Got: %s %d %s %d %s %s", rec.ReadName, rec.Flag, rec.RefName, rec.Pos, rec.Cigar, rec.Seq)
		if rec.Seq != "ACGTACGTAC" {
			t.Errorf("Sequence mismatch: got %q, want %q", rec.Seq, "ACGTACGTAC")
		}
		if rec.Cigar != "10M" {
			t.Errorf("CIGAR mismatch: got %q, want %q", rec.Cigar, "10M")
		}
	}
}

func TestReadCRAMHeader(t *testing.T) {
	cramFile := "testdata/test.cram"
	refFile := "testdata/ref.fa"

	reader, err := NewReader(cramFile, refFile)
	if err != nil {
		t.Fatalf("NewReader failed: %v", err)
	}
	defer reader.Close()

	hdr, err := reader.Header()
	if err != nil {
		t.Fatalf("Header failed: %v", err)
	}

	// Check header has expected lines.
	foundHD := false
	for _, line := range hdr.Lines {
		if strings.HasPrefix(line, "@HD") {
			foundHD = true
		}
	}
	if !foundHD {
		t.Error("missing @HD header line")
	}
}
