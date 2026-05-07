package cram

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgen-io/cgltk/htsio"
)

func TestWriterRoundtrip(t *testing.T) {
	refFile := "testdata/ref.fa"

	// Build a header.
	header := htsio.NewSamHeader()
	header.AddLine("@HD\tVN:1.6\tSO:coordinate")
	header.AddLine("@SQ\tSN:chr1\tLN:100000")
	header.AddLine("@SQ\tSN:chr2\tLN:50000")
	header.AddLine("@RG\tID:sample1\tSM:sample1")

	records := []*htsio.SamRecord{
		{
			ReadName: "read1", Flag: 0, RefName: "chr1", Pos: 100, MapQ: 60,
			Cigar: "10M", RefNext: "*", PosNext: 0, InsertLen: 0,
			Seq: "ACGTACGTAC", Qual: "IIIIIIIIII",
			Tags:     map[string]htsio.SamTag{"NM": {Type: 'i', Value: "0"}},
			TagOrder: []string{"NM"},
		},
		{
			ReadName: "read2", Flag: 16, RefName: "chr1", Pos: 500, MapQ: 30,
			Cigar: "5M1I4M", RefNext: "*", PosNext: 0, InsertLen: 0,
			Seq: "ACGTAACGTA", Qual: "IIIIIIIIII",
			Tags:     map[string]htsio.SamTag{"NM": {Type: 'i', Value: "1"}},
			TagOrder: []string{"NM"},
		},
		{
			ReadName: "read3", Flag: 4, RefName: "*", Pos: 0, MapQ: 0,
			Cigar: "*", RefNext: "*", PosNext: 0, InsertLen: 0,
			Seq: "NNNNNNNNNN", Qual: "IIIIIIIIII",
			Tags:     map[string]htsio.SamTag{},
			TagOrder: []string{},
		},
	}

	for _, version := range []struct {
		name string
		ver  Version
	}{
		{"v2.1", V2},
		{"v3.0", V3},
		{"v3.1", V31},
	} {
		t.Run(version.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			cramFile := filepath.Join(tmpDir, "test.cram")

			opts := NewWriterOpts().SetVersion(version.ver).Reference(refFile)
			w, err := NewWriter(cramFile, header, opts)
			if err != nil {
				t.Fatalf("NewWriter: %v", err)
			}

			for _, rec := range records {
				if err := w.Write(rec); err != nil {
					t.Fatalf("Write: %v", err)
				}
			}
			if err := w.Close(); err != nil {
				t.Fatalf("Close: %v", err)
			}

			// Read back with our reader.
			reader, err := NewReader(cramFile, refFile)
			if err != nil {
				t.Fatalf("NewReader: %v", err)
			}
			defer reader.Close()

			var gotRecords []string
			for rec, err := range reader.Records() {
				if err != nil {
					t.Fatalf("Records: %v", err)
				}
				line := fmt.Sprintf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
					rec.ReadName, rec.Flag, rec.RefName, rec.Pos, rec.MapQ,
					rec.Cigar, rec.RefNext, rec.PosNext, rec.InsertLen, rec.Seq, rec.Qual)
				gotRecords = append(gotRecords, line)
			}

			if len(gotRecords) != len(records) {
				t.Fatalf("record count: got %d, want %d", len(gotRecords), len(records))
			}

			// Check core fields match.
			for i, rec := range records {
				expected := fmt.Sprintf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
					rec.ReadName, rec.Flag, rec.RefName, rec.Pos, rec.MapQ,
					rec.Cigar, rec.RefNext, rec.PosNext, rec.InsertLen, rec.Seq, rec.Qual)
				if gotRecords[i] != expected {
					t.Errorf("record %d:\n  got:  %s\n  want: %s", i, gotRecords[i], expected)
				}
			}
		})
	}
}

func TestWriterFileRoundtrip(t *testing.T) {
	// Read records from existing CRAM, write as each version, read back, compare.
	// Uses a small recordsPerSlice to force multiple containers.
	refFile := "testdata/ref.fa"
	inputFile := "testdata/test.cram"

	// Read all records from the source file.
	reader, err := NewReader(inputFile, refFile)
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}

	var origRecords []*htsio.SamRecord
	for rec, err := range reader.Records() {
		if err != nil {
			t.Fatalf("reading source: %v", err)
		}
		origRecords = append(origRecords, rec)
	}
	reader.Close()

	if len(origRecords) == 0 {
		t.Fatal("no records in source file")
	}
	t.Logf("Source file has %d records", len(origRecords))

	origHeader, err := reader.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}

	// Also read from the 500-record test file for a more thorough test.
	reader500, err := NewReader("testdata/test_v31_500.cram", refFile)
	if err != nil {
		t.Fatalf("NewReader(500): %v", err)
	}
	var manyRecords []*htsio.SamRecord
	for rec, err := range reader500.Records() {
		if err != nil {
			t.Fatalf("reading 500-record source: %v", err)
		}
		manyRecords = append(manyRecords, rec)
	}
	reader500.Close()
	t.Logf("Large source has %d records", len(manyRecords))

	for _, version := range []struct {
		name string
		ver  Version
	}{
		{"v2.1", V2},
		{"v3.0", V3},
		{"v3.1", V31},
	} {
		t.Run(version.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			cramFile := filepath.Join(tmpDir, "roundtrip.cram")

			// Write all records with small slice size to force multiple containers.
			opts := NewWriterOpts().SetVersion(version.ver).Reference(refFile).RecordsPerSlice(50)
			w, err := NewWriter(cramFile, origHeader, opts)
			if err != nil {
				t.Fatalf("NewWriter: %v", err)
			}
			for _, rec := range manyRecords {
				if err := w.Write(rec); err != nil {
					t.Fatalf("Write: %v", err)
				}
			}
			if err := w.Close(); err != nil {
				t.Fatalf("Close: %v", err)
			}

			stat, _ := os.Stat(cramFile)
			t.Logf("wrote %d records, file size: %d bytes", len(manyRecords), stat.Size())

			// Read back.
			reader2, err := NewReader(cramFile, refFile)
			if err != nil {
				t.Fatalf("NewReader: %v", err)
			}
			defer reader2.Close()

			var gotRecords []*htsio.SamRecord
			for rec, err := range reader2.Records() {
				if err != nil {
					t.Fatalf("reading roundtrip: %v", err)
				}
				gotRecords = append(gotRecords, rec)
			}

			if len(gotRecords) != len(manyRecords) {
				t.Fatalf("record count: got %d, want %d", len(gotRecords), len(manyRecords))
			}

			// Compare each record's core fields.
			for i, orig := range manyRecords {
				got := gotRecords[i]
				if got.ReadName != orig.ReadName {
					t.Errorf("record %d ReadName: got %q, want %q", i, got.ReadName, orig.ReadName)
				}
				if got.Flag != orig.Flag {
					t.Errorf("record %d Flag: got %d, want %d", i, got.Flag, orig.Flag)
				}
				if got.RefName != orig.RefName {
					t.Errorf("record %d RefName: got %q, want %q", i, got.RefName, orig.RefName)
				}
				if got.Pos != orig.Pos {
					t.Errorf("record %d Pos: got %d, want %d", i, got.Pos, orig.Pos)
				}
				if got.MapQ != orig.MapQ {
					t.Errorf("record %d MapQ: got %d, want %d", i, got.MapQ, orig.MapQ)
				}
				if got.Cigar != orig.Cigar {
					t.Errorf("record %d Cigar: got %q, want %q", i, got.Cigar, orig.Cigar)
				}
				if got.Seq != orig.Seq {
					t.Errorf("record %d Seq: got %q, want %q", i, got.Seq, orig.Seq)
				}
			}

			// Also verify samtools can read it (if available).
			if _, err := exec.LookPath("samtools"); err == nil {
				absRefFile, _ := filepath.Abs(refFile)
				cmd := exec.Command("samtools", "view", "-T", absRefFile, cramFile)
				out, err := cmd.CombinedOutput()
				if err != nil {
					t.Errorf("samtools view failed: %v\noutput: %s", err, string(out))
				} else {
					lines := strings.Split(strings.TrimSpace(string(out)), "\n")
					if len(lines) != len(manyRecords) {
						t.Errorf("samtools got %d records, want %d", len(lines), len(manyRecords))
					}
				}
			}
		})
	}
}

func TestWriterSamtoolsReadable(t *testing.T) {
	// Skip if samtools not available.
	if _, err := exec.LookPath("samtools"); err != nil {
		t.Skip("samtools not found")
	}

	refFile := "testdata/ref.fa"

	header := htsio.NewSamHeader()
	header.AddLine("@HD\tVN:1.6\tSO:coordinate")
	header.AddLine("@SQ\tSN:chr1\tLN:10000")
	header.AddLine("@RG\tID:sample1\tSM:sample1")

	records := []*htsio.SamRecord{
		{
			ReadName: "read1", Flag: 0, RefName: "chr1", Pos: 100, MapQ: 60,
			Cigar: "10M", RefNext: "*", PosNext: 0, InsertLen: 0,
			Seq: "ACGTACGTAC", Qual: "IIIIIIIIII",
			Tags:     map[string]htsio.SamTag{},
			TagOrder: []string{},
		},
	}

	for _, version := range []struct {
		name string
		ver  Version
	}{
		{"v2.1", V2},
		{"v3.0", V3},
		{"v3.1", V31},
	} {
		t.Run(version.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			cramFile := filepath.Join(tmpDir, "test.cram")

			// Write absolute reference path in header so samtools can find it.
			absRefFile, _ := filepath.Abs(refFile)
			hdr := htsio.NewSamHeader()
			hdr.AddLine("@HD\tVN:1.6\tSO:coordinate")
			hdr.AddLine(fmt.Sprintf("@SQ\tSN:chr1\tLN:100000\tUR:%s", absRefFile))

			opts := NewWriterOpts().SetVersion(version.ver).Reference(refFile)
			w, err := NewWriter(cramFile, hdr, opts)
			if err != nil {
				t.Fatalf("NewWriter: %v", err)
			}

			for _, rec := range records {
				if err := w.Write(rec); err != nil {
					t.Fatalf("Write: %v", err)
				}
			}
			if err := w.Close(); err != nil {
				t.Fatalf("Close: %v", err)
			}

			// Verify samtools can read it.
			cmd := exec.Command("samtools", "view", "-T", absRefFile, cramFile)
			out, err := cmd.CombinedOutput()
			if err != nil {
				// Log the file for debugging.
				stat, _ := os.Stat(cramFile)
				t.Logf("CRAM file size: %d", stat.Size())
				t.Fatalf("samtools view failed: %v\noutput: %s", err, string(out))
			}

			lines := strings.Split(strings.TrimSpace(string(out)), "\n")
			if len(lines) != len(records) {
				t.Errorf("samtools got %d records, want %d\noutput: %s", len(lines), len(records), string(out))
			}
		})
	}
}
