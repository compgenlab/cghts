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
