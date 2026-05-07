package htsio_test

import (
	"crypto/sha1"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgen-io/cgkit/htsio"
	"github.com/compgen-io/cgkit/htsio/bam"
	_ "github.com/compgen-io/cgkit/htsio/cram"
	_ "github.com/compgen-io/cgkit/htsio/sam"
)

// samtoolsViewSHA1 returns the SHA1 of `samtools view <file>` output (records only).
func samtoolsViewSHA1(t *testing.T, filename string) string {
	t.Helper()
	cmd := exec.Command("samtools", "view", filename)
	out, err := cmd.Output()
	if err != nil {
		t.Fatalf("samtools view %s: %v", filename, err)
	}
	return fmt.Sprintf("%x", sha1.Sum(out))
}

// readAllRecords reads all records from a htsio.SamReader.
func readAllRecords(t *testing.T, reader htsio.SamReader) []*htsio.SamRecord {
	t.Helper()
	var records []*htsio.SamRecord
	for rec, err := range reader.Records() {
		if err != nil {
			t.Fatalf("Records: %v", err)
		}
		records = append(records, rec)
	}
	return records
}

// writeSAMFile writes records as SAM text (with header) to a file.
func writeSAMFile(t *testing.T, filename string, header *htsio.SamHeader, records []*htsio.SamRecord) {
	t.Helper()
	f, err := os.Create(filename)
	if err != nil {
		t.Fatalf("create %s: %v", filename, err)
	}
	defer f.Close()

	if header != nil {
		text := header.Text()
		if text != "" {
			if _, err := f.WriteString(text); err != nil {
				t.Fatalf("write header: %v", err)
			}
		}
	}
	for _, rec := range records {
		if _, err := fmt.Fprintln(f, rec.String()); err != nil {
			t.Fatalf("write record: %v", err)
		}
	}
}

// writeBAMFile writes records as BAM to a file.
func writeBAMFile(t *testing.T, filename string, header *htsio.SamHeader, records []*htsio.SamRecord) {
	t.Helper()
	writer, err := bam.NewWriter(filename, header)
	if err != nil {
		t.Fatalf("bam.NewWriter: %v", err)
	}
	for _, rec := range records {
		if err := writer.Write(rec); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := writer.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
}

func TestRoundTripBAMToSAM(t *testing.T) {
	inputs := []string{
		"testdata/test.bam",
		"testdata/test_tags.bam",
	}
	for _, input := range inputs {
		t.Run(filepath.Base(input), func(t *testing.T) {
			expected := samtoolsViewSHA1(t, input)

			reader, err := htsio.NewSamReader(input)
			if err != nil {
				t.Fatalf("htsio.NewSamReader: %v", err)
			}
			header, _ := reader.Header()
			records := readAllRecords(t, reader)
			reader.Close()

			// Write as SAM and compare
			samFile := filepath.Join(t.TempDir(), "output.sam")
			writeSAMFile(t, samFile, header, records)
			got := samtoolsViewSHA1(t, samFile)
			if got != expected {
				t.Errorf("SAM output mismatch: got %s, want %s", got, expected)
			}
		})
	}
}

func TestRoundTripBAMToBAM(t *testing.T) {
	inputs := []string{
		"testdata/test.bam",
		"testdata/test_tags.bam",
	}
	for _, input := range inputs {
		t.Run(filepath.Base(input), func(t *testing.T) {
			expected := samtoolsViewSHA1(t, input)

			reader, err := htsio.NewSamReader(input)
			if err != nil {
				t.Fatalf("htsio.NewSamReader: %v", err)
			}
			header, _ := reader.Header()
			records := readAllRecords(t, reader)
			reader.Close()

			// Write as BAM and compare
			bamFile := filepath.Join(t.TempDir(), "output.bam")
			writeBAMFile(t, bamFile, header, records)
			got := samtoolsViewSHA1(t, bamFile)
			if got != expected {
				t.Errorf("BAM output mismatch: got %s, want %s", got, expected)
			}
		})
	}
}

func TestRoundTripSAMToSAM(t *testing.T) {
	inputs := []string{
		"testdata/test.sam",
		"testdata/test_tags.sam",
	}
	for _, input := range inputs {
		t.Run(filepath.Base(input), func(t *testing.T) {
			expected := samtoolsViewSHA1(t, input)

			reader, err := htsio.NewSamReader(input)
			if err != nil {
				t.Fatalf("htsio.NewSamReader: %v", err)
			}
			header, _ := reader.Header()
			records := readAllRecords(t, reader)
			reader.Close()

			// Write as SAM and compare
			samFile := filepath.Join(t.TempDir(), "output.sam")
			writeSAMFile(t, samFile, header, records)
			got := samtoolsViewSHA1(t, samFile)
			if got != expected {
				t.Errorf("SAM output mismatch: got %s, want %s", got, expected)
			}
		})
	}
}

func TestRoundTripSAMToBAM(t *testing.T) {
	inputs := []string{
		"testdata/test.sam",
		"testdata/test_tags.sam",
	}
	for _, input := range inputs {
		t.Run(filepath.Base(input), func(t *testing.T) {
			expected := samtoolsViewSHA1(t, input)

			reader, err := htsio.NewSamReader(input)
			if err != nil {
				t.Fatalf("htsio.NewSamReader: %v", err)
			}
			header, _ := reader.Header()
			records := readAllRecords(t, reader)
			reader.Close()

			// Write as BAM and compare
			bamFile := filepath.Join(t.TempDir(), "output.bam")
			writeBAMFile(t, bamFile, header, records)
			got := samtoolsViewSHA1(t, bamFile)
			if got != expected {
				t.Errorf("BAM output mismatch: got %s, want %s", got, expected)
			}
		})
	}
}

func TestRoundTripPreservesTagOrder(t *testing.T) {
	// Read test_tags.bam, write BAM, read back — TagOrder should be identical.
	reader, err := htsio.NewSamReader("testdata/test_tags.bam")
	if err != nil {
		t.Fatalf("htsio.NewSamReader: %v", err)
	}
	header, _ := reader.Header()
	records := readAllRecords(t, reader)
	reader.Close()

	bamFile := filepath.Join(t.TempDir(), "output.bam")
	writeBAMFile(t, bamFile, header, records)

	reader2, err := htsio.NewSamReader(bamFile)
	if err != nil {
		t.Fatalf("htsio.NewSamReader (re-read): %v", err)
	}
	records2 := readAllRecords(t, reader2)
	reader2.Close()

	if len(records) != len(records2) {
		t.Fatalf("record count mismatch: %d vs %d", len(records), len(records2))
	}

	for i := range records {
		r1 := records[i]
		r2 := records2[i]
		if r1.ReadName != r2.ReadName {
			t.Errorf("record %d: name mismatch %q vs %q", i, r1.ReadName, r2.ReadName)
			continue
		}
		order1 := strings.Join(r1.TagOrder, ",")
		order2 := strings.Join(r2.TagOrder, ",")
		if order1 != order2 {
			t.Errorf("record %d (%s): TagOrder mismatch %q vs %q", i, r1.ReadName, order1, order2)
		}
	}
}

func TestRoundTripFieldConsistency(t *testing.T) {
	// Read test_tags.bam, write BAM, read back — all fields should match.
	reader, err := htsio.NewSamReader("testdata/test_tags.bam")
	if err != nil {
		t.Fatalf("htsio.NewSamReader: %v", err)
	}
	header, _ := reader.Header()
	records := readAllRecords(t, reader)
	reader.Close()

	bamFile := filepath.Join(t.TempDir(), "output.bam")
	writeBAMFile(t, bamFile, header, records)

	reader2, err := htsio.NewSamReader(bamFile)
	if err != nil {
		t.Fatalf("htsio.NewSamReader (re-read): %v", err)
	}
	records2 := readAllRecords(t, reader2)
	reader2.Close()

	if len(records) != len(records2) {
		t.Fatalf("record count mismatch: %d vs %d", len(records), len(records2))
	}

	for i := range records {
		r1 := records[i]
		r2 := records2[i]

		if r1.String() != r2.String() {
			t.Errorf("record %d mismatch:\n  got:  %s\n  want: %s", i, r2.String(), r1.String())
		}
	}
}
