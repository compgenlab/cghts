package htsio

import (
	"io"
	"os"
	"testing"
)

func TestBamWriterRoundTrip(t *testing.T) {
	tmpFile := t.TempDir() + "/test.bam"

	header := NewSamHeader()
	header.AddLine("@HD\tVN:1.6\tSO:coordinate")
	header.AddLine("@SQ\tSN:chr1\tLN:100000")
	header.AddLine("@SQ\tSN:chr2\tLN:50000")

	writer, err := newBamWriter(tmpFile, header)
	if err != nil {
		t.Fatalf("newBamWriter: %v", err)
	}

	records := []*SamRecord{
		{
			ReadName:  "read1",
			Flag:      0,
			RefName:   "chr1",
			Pos:       100,
			MapQ:      60,
			Cigar:     "10M",
			RefNext:   "*",
			PosNext:   0,
			InsertLen: 0,
			Seq:       "ACGTACGTAC",
			Qual:      "??????????",
			Tags: map[string]SamTag{
				"RG": {Type: 'Z', Value: "sample1"},
				"NM": {Type: 'i', Value: "2"},
			},
		},
		{
			ReadName:  "read2",
			Flag:      16,
			RefName:   "chr2",
			Pos:       500,
			MapQ:      30,
			Cigar:     "5S10M2I3M",
			RefNext:   "*",
			PosNext:   0,
			InsertLen: 0,
			Seq:       "ACGTAACGTAACGTAACGTA",
			Qual:      "IIIIIIIIIIIIIIIIIIII",
			Tags:      map[string]SamTag{},
		},
		{
			ReadName:  "unmapped",
			Flag:      4,
			RefName:   "*",
			Pos:       0,
			MapQ:      0,
			Cigar:     "*",
			RefNext:   "*",
			PosNext:   0,
			InsertLen: 0,
			Seq:       "ACGT",
			Qual:      "*",
			Tags:      map[string]SamTag{},
		},
	}

	for _, rec := range records {
		if err := writer.Write(rec); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := writer.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	// Read back
	f, err := os.Open(tmpFile)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	reader, err := NewBamReader(f)
	if err != nil {
		t.Fatalf("NewBamReader: %v", err)
	}
	defer reader.Close()

	hdr, err := reader.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	refs := hdr.References()
	if len(refs) != 2 {
		t.Fatalf("expected 2 refs, got %d", len(refs))
	}
	if refs[0].Name != "chr1" || refs[1].Name != "chr2" {
		t.Errorf("refs: %v", refs)
	}

	// Read record 1
	r, err := reader.Next()
	if err != nil {
		t.Fatalf("Next[0]: %v", err)
	}
	if r.ReadName != "read1" {
		t.Errorf("ReadName: got %q, want %q", r.ReadName, "read1")
	}
	if r.RefName != "chr1" {
		t.Errorf("RefName: got %q, want %q", r.RefName, "chr1")
	}
	if r.Pos != 100 {
		t.Errorf("Pos: got %d, want 100", r.Pos)
	}
	if r.Cigar != "10M" {
		t.Errorf("Cigar: got %q, want %q", r.Cigar, "10M")
	}
	if r.Seq != "ACGTACGTAC" {
		t.Errorf("Seq: got %q, want %q", r.Seq, "ACGTACGTAC")
	}
	if r.Qual != "??????????" {
		t.Errorf("Qual: got %q, want %q", r.Qual, "??????????")
	}
	rg, ok := r.Tags["RG"]
	if !ok || rg.Value != "sample1" {
		t.Errorf("RG tag: %v", r.Tags["RG"])
	}
	nm, ok := r.Tags["NM"]
	if !ok || nm.Value != "2" {
		t.Errorf("NM tag: %v", r.Tags["NM"])
	}

	// Read record 2
	r, err = reader.Next()
	if err != nil {
		t.Fatalf("Next[1]: %v", err)
	}
	if r.ReadName != "read2" {
		t.Errorf("ReadName: got %q, want %q", r.ReadName, "read2")
	}
	if r.Cigar != "5S10M2I3M" {
		t.Errorf("Cigar: got %q, want %q", r.Cigar, "5S10M2I3M")
	}
	if !r.IsReverse() {
		t.Error("expected reverse strand")
	}

	// Read record 3 (unmapped)
	r, err = reader.Next()
	if err != nil {
		t.Fatalf("Next[2]: %v", err)
	}
	if r.ReadName != "unmapped" {
		t.Errorf("ReadName: got %q, want %q", r.ReadName, "unmapped")
	}
	if r.RefName != "*" {
		t.Errorf("RefName: got %q, want *", r.RefName)
	}
	if !r.IsUnmapped() {
		t.Error("expected unmapped flag")
	}

	// EOF
	_, err = reader.Next()
	if err != io.EOF {
		t.Errorf("expected io.EOF, got %v", err)
	}
}

func TestBamWriterMateOnSameRef(t *testing.T) {
	tmpFile := t.TempDir() + "/test_mate.bam"

	header := NewSamHeader()
	header.AddLine("@HD\tVN:1.6")
	header.AddLine("@SQ\tSN:chr1\tLN:100000")

	writer, err := newBamWriter(tmpFile, header)
	if err != nil {
		t.Fatalf("newBamWriter: %v", err)
	}

	rec := &SamRecord{
		ReadName:  "paired",
		Flag:      99, // paired, proper pair, mate reverse, first in pair
		RefName:   "chr1",
		Pos:       100,
		MapQ:      60,
		Cigar:     "50M",
		RefNext:   "=",
		PosNext:   300,
		InsertLen: 250,
		Seq:       "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC",
		Qual:      "*",
		Tags:      map[string]SamTag{},
	}

	if err := writer.Write(rec); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := writer.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	f, err := os.Open(tmpFile)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	reader, err := NewBamReader(f)
	if err != nil {
		t.Fatalf("NewBamReader: %v", err)
	}
	defer reader.Close()

	r, err := reader.Next()
	if err != nil {
		t.Fatalf("Next: %v", err)
	}
	if r.RefNext != "=" {
		t.Errorf("RefNext: got %q, want =", r.RefNext)
	}
	if r.PosNext != 300 {
		t.Errorf("PosNext: got %d, want 300", r.PosNext)
	}
	if r.InsertLen != 250 {
		t.Errorf("InsertLen: got %d, want 250", r.InsertLen)
	}
}
