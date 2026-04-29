package htsio

import (
	"bytes"
	"encoding/binary"
	"os"
	"testing"

	"github.com/compgen-io/cgltk/htsio/bgzf"
)

// writeSyntheticBAI creates a minimal BAI file with one reference sequence
// containing the given bins and linear index.
func writeSyntheticBAI(bins map[uint32][]Chunk, linearIdx []bgzf.VirtualOffset) []byte {
	var buf bytes.Buffer

	// Magic
	buf.Write([]byte("BAI\x01"))

	// n_ref = 1
	binary.Write(&buf, binary.LittleEndian, int32(1))

	// Bins
	binary.Write(&buf, binary.LittleEndian, int32(len(bins)))
	for binNum, chunks := range bins {
		binary.Write(&buf, binary.LittleEndian, binNum)
		binary.Write(&buf, binary.LittleEndian, int32(len(chunks)))
		for _, c := range chunks {
			binary.Write(&buf, binary.LittleEndian, uint64(c.Begin))
			binary.Write(&buf, binary.LittleEndian, uint64(c.End))
		}
	}

	// Linear index
	binary.Write(&buf, binary.LittleEndian, int32(len(linearIdx)))
	for _, off := range linearIdx {
		binary.Write(&buf, binary.LittleEndian, uint64(off))
	}

	return buf.Bytes()
}

func TestLoadBAI(t *testing.T) {
	// A record at position 10000 with ref length 100 → bin 4681 + (10000>>14) = 4681
	binNum := uint32(4681) // bin for 0-16383
	chunk := Chunk{
		Begin: bgzf.NewVirtualOffset(1000, 0),
		End:   bgzf.NewVirtualOffset(2000, 0),
	}
	linear := []bgzf.VirtualOffset{
		bgzf.NewVirtualOffset(1000, 0), // window 0 (0-16383)
	}

	data := writeSyntheticBAI(
		map[uint32][]Chunk{binNum: {chunk}},
		linear,
	)

	tmpFile := t.TempDir() + "/test.bai"
	if err := os.WriteFile(tmpFile, data, 0644); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}

	idx, err := LoadBAI(tmpFile)
	if err != nil {
		t.Fatalf("LoadBAI: %v", err)
	}

	if len(idx.refs) != 1 {
		t.Fatalf("expected 1 ref, got %d", len(idx.refs))
	}
	if len(idx.refs[0].bins) != 1 {
		t.Errorf("expected 1 bin, got %d", len(idx.refs[0].bins))
	}
	if len(idx.refs[0].linearIdx) != 1 {
		t.Errorf("expected 1 linear entry, got %d", len(idx.refs[0].linearIdx))
	}
}

func TestReg2bins(t *testing.T) {
	// Region [0, 16384) should include bin 0 and bin 4681
	bins := reg2bins(0, 16384)
	found := make(map[uint32]bool)
	for _, b := range bins {
		found[b] = true
	}
	if !found[0] {
		t.Error("expected bin 0")
	}
	if !found[4681] {
		t.Error("expected bin 4681")
	}
}

func TestBinIndexQuery(t *testing.T) {
	// Set up an index with chunks in bins that overlap region [5000, 10000).
	// Bin 4681 covers [0, 16384).
	bin4681 := uint32(4681)
	chunk1 := Chunk{
		Begin: bgzf.NewVirtualOffset(100, 0),
		End:   bgzf.NewVirtualOffset(200, 0),
	}
	chunk2 := Chunk{
		Begin: bgzf.NewVirtualOffset(300, 0),
		End:   bgzf.NewVirtualOffset(400, 0),
	}

	linear := []bgzf.VirtualOffset{
		bgzf.NewVirtualOffset(50, 0), // window 0
	}

	data := writeSyntheticBAI(
		map[uint32][]Chunk{bin4681: {chunk1, chunk2}},
		linear,
	)

	tmpFile := t.TempDir() + "/test.bai"
	os.WriteFile(tmpFile, data, 0644)

	idx, err := LoadBAI(tmpFile)
	if err != nil {
		t.Fatalf("LoadBAI: %v", err)
	}

	chunks := idx.Query(0, 5000, 10000)
	if len(chunks) != 2 {
		t.Fatalf("expected 2 chunks, got %d", len(chunks))
	}

	// Chunks should be sorted by Begin.
	if chunks[0].Begin > chunks[1].Begin {
		t.Error("chunks not sorted")
	}
}

func TestBinIndexQueryMerge(t *testing.T) {
	// Two overlapping chunks should be merged.
	bin4681 := uint32(4681)
	chunk1 := Chunk{
		Begin: bgzf.NewVirtualOffset(100, 0),
		End:   bgzf.NewVirtualOffset(300, 0),
	}
	chunk2 := Chunk{
		Begin: bgzf.NewVirtualOffset(200, 0), // overlaps with chunk1
		End:   bgzf.NewVirtualOffset(400, 0),
	}

	linear := []bgzf.VirtualOffset{
		bgzf.NewVirtualOffset(0, 0),
	}

	data := writeSyntheticBAI(
		map[uint32][]Chunk{bin4681: {chunk1, chunk2}},
		linear,
	)

	tmpFile := t.TempDir() + "/test.bai"
	os.WriteFile(tmpFile, data, 0644)

	idx, _ := LoadBAI(tmpFile)
	chunks := idx.Query(0, 0, 16384)

	if len(chunks) != 1 {
		t.Fatalf("expected 1 merged chunk, got %d", len(chunks))
	}
	if chunks[0].Begin != bgzf.NewVirtualOffset(100, 0) {
		t.Errorf("merged begin: got %d, want %d", chunks[0].Begin, bgzf.NewVirtualOffset(100, 0))
	}
	if chunks[0].End != bgzf.NewVirtualOffset(400, 0) {
		t.Errorf("merged end: got %d, want %d", chunks[0].End, bgzf.NewVirtualOffset(400, 0))
	}
}

func TestBinIndexQueryLinearFilter(t *testing.T) {
	// A chunk that ends before the linear index minimum should be filtered out.
	bin4681 := uint32(4681)
	earlyChunk := Chunk{
		Begin: bgzf.NewVirtualOffset(10, 0),
		End:   bgzf.NewVirtualOffset(50, 0),
	}
	laterChunk := Chunk{
		Begin: bgzf.NewVirtualOffset(200, 0),
		End:   bgzf.NewVirtualOffset(400, 0),
	}

	linear := []bgzf.VirtualOffset{
		bgzf.NewVirtualOffset(100, 0), // minimum for window 0
	}

	data := writeSyntheticBAI(
		map[uint32][]Chunk{bin4681: {earlyChunk, laterChunk}},
		linear,
	)

	tmpFile := t.TempDir() + "/test.bai"
	os.WriteFile(tmpFile, data, 0644)

	idx, _ := LoadBAI(tmpFile)
	chunks := idx.Query(0, 0, 16384)

	// earlyChunk should be filtered out (ends at 50, before linear min 100).
	if len(chunks) != 1 {
		t.Fatalf("expected 1 chunk after linear filter, got %d", len(chunks))
	}
	if chunks[0].Begin != bgzf.NewVirtualOffset(200, 0) {
		t.Errorf("expected later chunk, got begin=%d", chunks[0].Begin)
	}
}

func TestBinIndexQueryOutOfRange(t *testing.T) {
	data := writeSyntheticBAI(
		map[uint32][]Chunk{},
		nil,
	)

	tmpFile := t.TempDir() + "/test.bai"
	os.WriteFile(tmpFile, data, 0644)

	idx, _ := LoadBAI(tmpFile)

	// Query non-existent ref.
	chunks := idx.Query(5, 0, 100)
	if chunks != nil {
		t.Errorf("expected nil for out-of-range refID, got %v", chunks)
	}

	// Query existing ref but no bins.
	chunks = idx.Query(0, 0, 100)
	if chunks != nil {
		t.Errorf("expected nil for empty bins, got %v", chunks)
	}
}

func TestSplitNulTerminated(t *testing.T) {
	data := []byte("chr1\x00chr2\x00chr3\x00")
	names := splitNulTerminated(data)
	if len(names) != 3 {
		t.Fatalf("expected 3 names, got %d", len(names))
	}
	if names[0] != "chr1" || names[1] != "chr2" || names[2] != "chr3" {
		t.Errorf("names: %v", names)
	}
}
