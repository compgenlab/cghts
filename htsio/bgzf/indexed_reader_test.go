package bgzf

import (
	"bytes"
	"encoding/binary"
	"io"
	"os"
	"testing"
)

// writeBgzfBlocks writes multiple BGZF blocks with explicit Flush boundaries
// and returns the data plus the compressed offset of each block.
func writeBgzfBlocks(blocks [][]byte) ([]byte, []int64) {
	var buf bytes.Buffer
	w := NewWriter(&buf)
	var offsets []int64

	for _, block := range blocks {
		offsets = append(offsets, int64(buf.Len()))
		w.Write(block)
		w.Flush()
	}
	w.Close()

	return buf.Bytes(), offsets
}

func TestIndexedReaderSeek(t *testing.T) {
	block1 := []byte("Hello from block one!")
	block2 := []byte("Greetings from block two!")
	block3 := []byte("Welcome to block three!")

	data, offsets := writeBgzfBlocks([][]byte{block1, block2, block3})

	rs := bytes.NewReader(data)
	ir := NewIndexedReader(rs)

	// Seek to start of block2 and read it.
	vo := NewVirtualOffset(offsets[1], 0)
	if err := ir.SeekToVirtualOffset(vo); err != nil {
		t.Fatalf("SeekToVirtualOffset: %v", err)
	}

	got := make([]byte, len(block2))
	if _, err := io.ReadFull(ir, got); err != nil {
		t.Fatalf("ReadFull: %v", err)
	}
	if string(got) != string(block2) {
		t.Errorf("got %q, want %q", got, block2)
	}

	// Seek to middle of block1.
	vo = NewVirtualOffset(offsets[0], 6)
	if err := ir.SeekToVirtualOffset(vo); err != nil {
		t.Fatalf("SeekToVirtualOffset: %v", err)
	}
	got = make([]byte, len("from block one!"))
	if _, err := io.ReadFull(ir, got); err != nil {
		t.Fatalf("ReadFull: %v", err)
	}
	if string(got) != "from block one!" {
		t.Errorf("got %q, want %q", got, "from block one!")
	}
}

func TestIndexedReaderSequentialAcrossBlocks(t *testing.T) {
	block1 := []byte("AAAA")
	block2 := []byte("BBBB")
	block3 := []byte("CCCC")

	data, offsets := writeBgzfBlocks([][]byte{block1, block2, block3})

	rs := bytes.NewReader(data)
	ir := NewIndexedReader(rs)

	// Seek to start of block1, then read across all three blocks.
	if err := ir.SeekToVirtualOffset(NewVirtualOffset(offsets[0], 0)); err != nil {
		t.Fatalf("Seek: %v", err)
	}

	got, err := io.ReadAll(ir)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if string(got) != "AAAABBBBCCCC" {
		t.Errorf("got %q, want %q", got, "AAAABBBBCCCC")
	}
}

func TestIndexedReaderCacheHit(t *testing.T) {
	block1 := []byte("first block data")
	block2 := []byte("second block data")

	data, offsets := writeBgzfBlocks([][]byte{block1, block2})

	rs := bytes.NewReader(data)
	ir := NewIndexedReaderSize(rs, 4)

	// Read block1.
	if err := ir.SeekToVirtualOffset(NewVirtualOffset(offsets[0], 0)); err != nil {
		t.Fatalf("Seek: %v", err)
	}
	got := make([]byte, len(block1))
	io.ReadFull(ir, got)

	// Seek to block2.
	if err := ir.SeekToVirtualOffset(NewVirtualOffset(offsets[1], 0)); err != nil {
		t.Fatalf("Seek: %v", err)
	}
	got = make([]byte, len(block2))
	io.ReadFull(ir, got)

	// Seek back to block1 — should be a cache hit (no file seek needed,
	// though we can't directly observe that, we verify correctness).
	if err := ir.SeekToVirtualOffset(NewVirtualOffset(offsets[0], 0)); err != nil {
		t.Fatalf("Seek: %v", err)
	}
	got = make([]byte, len(block1))
	if _, err := io.ReadFull(ir, got); err != nil {
		t.Fatalf("ReadFull: %v", err)
	}
	if string(got) != string(block1) {
		t.Errorf("got %q, want %q", got, block1)
	}
}

func TestIndexedReaderVirtualTell(t *testing.T) {
	block1 := []byte("ABCDEF")
	block2 := []byte("GHIJKL")

	data, offsets := writeBgzfBlocks([][]byte{block1, block2})

	rs := bytes.NewReader(data)
	ir := NewIndexedReader(rs)

	// Seek to block1, offset 3.
	vo := NewVirtualOffset(offsets[0], 3)
	if err := ir.SeekToVirtualOffset(vo); err != nil {
		t.Fatalf("Seek: %v", err)
	}

	tell := ir.VirtualTell()
	if tell != vo {
		t.Errorf("VirtualTell: got %d, want %d", tell, vo)
	}

	// Read one byte, tell should advance.
	ir.ReadByte()
	tell = ir.VirtualTell()
	expected := NewVirtualOffset(offsets[0], 4)
	if tell != expected {
		t.Errorf("VirtualTell after read: got %d, want %d", tell, expected)
	}
}

func TestIndexedReaderCacheEviction(t *testing.T) {
	// Create more blocks than the cache can hold.
	blocks := make([][]byte, 5)
	for i := range blocks {
		blocks[i] = []byte{byte('A' + i), byte('A' + i), byte('A' + i), byte('A' + i)}
	}

	data, offsets := writeBgzfBlocks(blocks)

	// Cache size of 2 — blocks 0,1 get evicted when 2,3 are loaded.
	rs := bytes.NewReader(data)
	ir := NewIndexedReaderSize(rs, 2)

	// Read through all blocks sequentially.
	for i := 0; i < len(blocks); i++ {
		if err := ir.SeekToVirtualOffset(NewVirtualOffset(offsets[i], 0)); err != nil {
			t.Fatalf("Seek[%d]: %v", i, err)
		}
		got := make([]byte, 4)
		if _, err := io.ReadFull(ir, got); err != nil {
			t.Fatalf("ReadFull[%d]: %v", i, err)
		}
		if string(got) != string(blocks[i]) {
			t.Errorf("block[%d]: got %q, want %q", i, got, blocks[i])
		}
	}

	// Seek back to block 0 — should have been evicted but still works
	// (just requires a re-read from file).
	if err := ir.SeekToVirtualOffset(NewVirtualOffset(offsets[0], 0)); err != nil {
		t.Fatalf("Seek back to 0: %v", err)
	}
	got := make([]byte, 4)
	if _, err := io.ReadFull(ir, got); err != nil {
		t.Fatalf("ReadFull after eviction: %v", err)
	}
	if string(got) != string(blocks[0]) {
		t.Errorf("after eviction: got %q, want %q", got, blocks[0])
	}
}

func TestIndexedReaderSeekWithGZI(t *testing.T) {
	block1 := []byte("AAAAAAAAAA") // 10 bytes
	block2 := []byte("BBBBBBBBBB") // 10 bytes
	block3 := []byte("CCCCCCCCCC") // 10 bytes

	data, offsets := writeBgzfBlocks([][]byte{block1, block2, block3})

	// Build a GZIndex manually.
	idx := &GZIndex{
		entries: []gziEntry{
			{compressedOffset: offsets[0], uncompressedOffset: 0},
			{compressedOffset: offsets[1], uncompressedOffset: 10},
			{compressedOffset: offsets[2], uncompressedOffset: 20},
		},
	}

	rs := bytes.NewReader(data)
	ir := NewIndexedReader(rs)
	ir.setGZIndex(idx)

	// Seek to uncompressed offset 15 → block2, within-block 5 → "BBBBB"
	pos, err := ir.Seek(15, io.SeekStart)
	if err != nil {
		t.Fatalf("Seek: %v", err)
	}
	if pos != 15 {
		t.Errorf("Seek returned %d, want 15", pos)
	}

	got := make([]byte, 5)
	if _, err := io.ReadFull(ir, got); err != nil {
		t.Fatalf("ReadFull: %v", err)
	}
	if string(got) != "BBBBB" {
		t.Errorf("got %q, want %q", got, "BBBBB")
	}

	// Seek to uncompressed offset 0 → start
	_, err = ir.Seek(0, io.SeekStart)
	if err != nil {
		t.Fatalf("Seek(0): %v", err)
	}
	got = make([]byte, 10)
	if _, err := io.ReadFull(ir, got); err != nil {
		t.Fatalf("ReadFull: %v", err)
	}
	if string(got) != "AAAAAAAAAA" {
		t.Errorf("got %q, want %q", got, "AAAAAAAAAA")
	}

	// Seek across blocks: read from offset 8 through block boundary
	_, err = ir.Seek(8, io.SeekStart)
	if err != nil {
		t.Fatalf("Seek(8): %v", err)
	}
	got = make([]byte, 4)
	if _, err := io.ReadFull(ir, got); err != nil {
		t.Fatalf("ReadFull: %v", err)
	}
	if string(got) != "AABB" {
		t.Errorf("got %q, want %q", got, "AABB")
	}
}

func TestIndexedReaderSeekWithoutGZI(t *testing.T) {
	block := []byte("test data")
	data, _ := writeBgzfBlocks([][]byte{block})

	rs := bytes.NewReader(data)
	ir := NewIndexedReader(rs)

	_, err := ir.Seek(0, io.SeekStart)
	if err == nil {
		t.Fatal("expected error when seeking without GZI")
	}
}

func TestIndexedReaderGZIAutoLoad(t *testing.T) {
	block1 := []byte("first block!")
	block2 := []byte("second block")

	data, offsets := writeBgzfBlocks([][]byte{block1, block2})

	dir := t.TempDir()
	bgzfPath := dir + "/test.bgz"
	gziPath := bgzfPath + ".gzi"

	// Write the BGZF file.
	if err := os.WriteFile(bgzfPath, data, 0644); err != nil {
		t.Fatalf("WriteFile bgzf: %v", err)
	}

	// Write a .gzi index.
	var gziBuf bytes.Buffer
	binary.Write(&gziBuf, binary.LittleEndian, uint64(1)) // count (excluding implicit block 0)
	binary.Write(&gziBuf, binary.LittleEndian, uint64(offsets[1]))
	binary.Write(&gziBuf, binary.LittleEndian, uint64(12)) // uncompressed offset of block2
	if err := os.WriteFile(gziPath, gziBuf.Bytes(), 0644); err != nil {
		t.Fatalf("WriteFile gzi: %v", err)
	}

	// Open with auto-loading.
	ir, err := OpenIndexedReader(bgzfPath)
	if err != nil {
		t.Fatalf("OpenIndexedReader: %v", err)
	}

	// Seek by uncompressed offset should work (gzi was auto-loaded).
	pos, err := ir.Seek(12, io.SeekStart)
	if err != nil {
		t.Fatalf("Seek: %v", err)
	}
	if pos != 12 {
		t.Errorf("Seek returned %d, want 12", pos)
	}

	got := make([]byte, len(block2))
	if _, err := io.ReadFull(ir, got); err != nil {
		t.Fatalf("ReadFull: %v", err)
	}
	if string(got) != string(block2) {
		t.Errorf("got %q, want %q", got, block2)
	}
}

func TestIndexedReaderReadByte(t *testing.T) {
	block := []byte{0xDE, 0xAD, 0xBE, 0xEF}
	data, offsets := writeBgzfBlocks([][]byte{block})

	rs := bytes.NewReader(data)
	ir := NewIndexedReader(rs)

	if err := ir.SeekToVirtualOffset(NewVirtualOffset(offsets[0], 0)); err != nil {
		t.Fatalf("Seek: %v", err)
	}

	for i, want := range block {
		got, err := ir.ReadByte()
		if err != nil {
			t.Fatalf("ReadByte[%d]: %v", i, err)
		}
		if got != want {
			t.Errorf("ReadByte[%d]: got %#x, want %#x", i, got, want)
		}
	}
}
