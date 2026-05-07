package htsio_test

import (
	"testing"

	"github.com/compgen-io/cgkit/htsio"
	_ "github.com/compgen-io/cgkit/htsio/bam"
	_ "github.com/compgen-io/cgkit/htsio/cram"
	_ "github.com/compgen-io/cgkit/htsio/sam"
)

func TestBamReaderQuery(t *testing.T) {
	reader, err := htsio.NewSamReader("testdata/test.bam")
	if err != nil {
		t.Fatalf("NewSamReader: %v", err)
	}
	defer reader.Close()

	records, err := reader.Query("chr1", 400, 600)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}

	var names []string
	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter: %v", err)
		}
		names = append(names, rec.ReadName)
	}

	if len(names) != 1 || names[0] != "read2" {
		t.Errorf("expected [read2], got %v", names)
	}

	records, err = reader.Query("chr2", 0, 1000)
	if err != nil {
		t.Fatalf("Query chr2: %v", err)
	}

	names = nil
	for rec, err := range records {
		if err != nil {
			t.Fatalf("iter: %v", err)
		}
		names = append(names, rec.ReadName)
	}

	if len(names) != 1 || names[0] != "read5" {
		t.Errorf("expected [read5], got %v", names)
	}
}
