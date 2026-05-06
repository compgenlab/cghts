package cram

import (
	"bytes"
	"os"
	"testing"
)

func TestRansDecompress(t *testing.T) {
	// Read the CRAM file and verify rANS-compressed blocks decompress correctly.
	// The expected values were verified against htslib's rans_uncompress() C implementation.
	data, err := os.ReadFile("testdata/test.cram")
	if err != nil {
		t.Fatal(err)
	}

	r := bytes.NewReader(data)

	// File definition (26 bytes).
	if _, err := readFileDefinition(r); err != nil {
		t.Fatal(err)
	}

	// Header container.
	hdrCont, err := readContainerHeader(r)
	if err != nil {
		t.Fatal(err)
	}
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

	for i := int32(0); i < dataCont.NumBlocks; i++ {
		blk, err := readBlock(r)
		if err != nil {
			t.Fatalf("reading block %d: %v", i, err)
		}

		if blk.contentID == 28 && blk.method == blockMethodRans4x8 {
			// FP (feature position) block — verified against htslib output
			expected := []byte{1, 1, 1, 1, 2, 2, 2, 1, 2, 1, 1, 2, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 4, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 3, 1}
			if len(blk.data) < len(expected) {
				t.Fatalf("FP block too short: got %d bytes, need %d", len(blk.data), len(expected))
			}
			for j, want := range expected {
				if blk.data[j] != want {
					t.Errorf("FP[%d]: got %d, want %d", j, blk.data[j], want)
				}
			}
		}

		if blk.contentID == 31 && blk.method == blockMethodRans4x8 {
			// BS (base substitution) block — verified against htslib rans_uncompress()
			expected := []byte{2, 1, 0, 2, 0, 2, 0, 2, 0, 2, 1, 1, 2, 2, 1, 1, 2, 2, 0, 1, 2, 1, 1, 2, 2, 0, 0, 1, 2, 1, 1, 0, 1, 1, 1, 2, 0}
			if len(blk.data) < len(expected) {
				t.Fatalf("BS block too short: got %d bytes, need %d", len(blk.data), len(expected))
			}
			for j, want := range expected {
				if blk.data[j] != want {
					t.Errorf("BS[%d]: got %d, want %d", j, blk.data[j], want)
				}
			}
		}
	}
}

func TestRansDecompressRaw(t *testing.T) {
	// test_raw.cram has some rANS blocks with symbol 255 (full byte range).
	// This tests the frequency table parser's handling of symbol 255
	// where byte(j+1) wraps to 0 and must not be confused with the terminator.
	data, err := os.ReadFile("testdata/test_raw.cram")
	if err != nil {
		t.Fatal(err)
	}

	r := bytes.NewReader(data)
	if _, err := readFileDefinition(r); err != nil {
		t.Fatal(err)
	}

	// Read all containers and verify no rANS decompression errors.
	for {
		ch, err := readContainerHeader(r)
		if err != nil {
			t.Fatal(err)
		}
		if ch.isEOF() {
			break
		}
		for i := int32(0); i < ch.NumBlocks; i++ {
			blk, err := readBlock(r)
			if err != nil {
				t.Fatalf("container block %d: %v", i, err)
			}
			if blk.method == blockMethodRans4x8 && len(blk.data) == 0 {
				t.Error("rANS block decompressed to empty data")
			}
		}
	}
}
