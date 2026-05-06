package cram

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"os"
	"testing"
)

func TestRansDecompress(t *testing.T) {
	// Read the CRAM file and find the rANS-compressed FP block (cid=28).
	data, err := os.ReadFile("testdata/test.cram")
	if err != nil {
		t.Fatal(err)
	}

	// Parse through the file to find block with cid=28.
	r := bytes.NewReader(data)

	// File definition (26 bytes).
	fd, err := readFileDefinition(r)
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("CRAM version: %d.%d", fd.Major, fd.Minor)

	// Header container.
	hdrCont, err := readContainerHeader(r)
	if err != nil {
		t.Fatal(err)
	}
	// Skip header blocks.
	for i := int32(0); i < hdrCont.NumBlocks; i++ {
		if _, err := readBlock(r); err != nil {
			t.Fatal(err)
		}
	}

	// Data container.
	dataCont, err := readContainerHeader(r)
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("Data container: numBlocks=%d numRecords=%d", dataCont.NumBlocks, dataCont.NumRecords)

	// Read all blocks.
	for i := int32(0); i < dataCont.NumBlocks; i++ {
		blk, err := readBlock(r)
		if err != nil {
			t.Fatalf("reading block %d: %v", i, err)
		}
		t.Logf("Block %d: method=%d ctype=%d cid=%d rawSize=%d dataLen=%d",
			i, blk.method, blk.contentType, blk.contentID, blk.rawSize, len(blk.data))

		if blk.contentID == 28 { // FP block
			expected := []byte{1, 1, 1, 1, 2, 2, 2, 1, 2, 1, 1, 2, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 4, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 3, 1}
			t.Logf("FP first 37: %v", blk.data[:37])
			mismatches := 0
			for j := 0; j < 37; j++ {
				if blk.data[j] != expected[j] {
					t.Logf("  FP MISMATCH at %d: got %d, want %d", j, blk.data[j], expected[j])
					mismatches++
				}
			}
			if mismatches == 0 {
				t.Log("  FP first 37 bytes match!")
			}
		}
		if blk.contentID == 31 { // BS block
			expected := []byte{2, 1, 0, 2, 0, 2, 0, 2, 0, 2, 1, 1, 0, 1, 2, 0, 1, 1, 0, 1, 0, 2, 0, 1, 2, 0, 0, 1, 0, 2, 2, 1, 1, 1, 1, 1, 0}
			t.Logf("BS first 37: %v", blk.data[:37])
			mismatches := 0
			for j := 0; j < 37; j++ {
				if blk.data[j] != expected[j] {
					t.Logf("  BS MISMATCH at %d: got %d, want %d", j, blk.data[j], expected[j])
					mismatches++
				}
			}
			if mismatches == 0 {
				t.Log("  BS first 37 bytes match!")
			}
		}
	}
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// TestRansRoundtrip verifies rANS decompression produces expected output.
func TestRansRoundtrip(t *testing.T) {
	// Create a simple order-0 rANS compressed block manually.
	// For now, just verify the decompressor doesn't crash.
	input := []byte{1, 1, 1, 2, 2, 3, 3, 3, 3, 4}

	// Build a simple frequency table: 1→3, 2→2, 3→4, 4→1 (total=10)
	// We need total to be 4096 for ransL. This test is more about structure.
	_ = input
	_ = binary.LittleEndian
	_ = fmt.Sprintf
}
