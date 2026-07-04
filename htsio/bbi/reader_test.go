package bbi

import (
	"encoding/binary"
	"math"
	"os"
	"path/filepath"
	"testing"
)

// --- minimal in-test BBI writers (single chrom, single uncompressed data block,
// single R-tree leaf). Enough to exercise the reader without external kent tools.

type wigItem struct {
	start, end uint32
	val        float32
}
type bedItem struct {
	start, end uint32
	rest       string // tab-separated remaining columns (name, score, …)
}

func buildBBI(magic uint32, chrom string, chromSize uint32, dataBlock []byte, spanStart, spanEnd uint32) []byte {
	le := binary.LittleEndian
	keySize := uint32(len(chrom))

	// chrom B+ tree: 32-byte header + one leaf node (4 + key+8).
	chromTree := make([]byte, 0, 48)
	{
		h := make([]byte, 32)
		le.PutUint32(h[0:], magicChromBP)
		le.PutUint32(h[4:], 1) // blockSize
		le.PutUint32(h[8:], keySize)
		le.PutUint32(h[12:], 8) // valSize
		le.PutUint64(h[16:], 1) // itemCount
		chromTree = append(chromTree, h...)
		node := make([]byte, 4)
		node[0] = 1 // isLeaf
		le.PutUint16(node[2:], 1)
		key := make([]byte, keySize)
		copy(key, chrom)
		item := make([]byte, keySize+8)
		copy(item, key)
		le.PutUint32(item[keySize:], 0) // chromId 0
		le.PutUint32(item[keySize+4:], chromSize)
		node = append(node, item...)
		chromTree = append(chromTree, node...)
	}

	const chromTreeOff = 64
	dataOff := uint64(chromTreeOff + len(chromTree))
	indexOff := dataOff + uint64(len(dataBlock))

	// R-tree: 48-byte header + one leaf node (4 + 32).
	rtree := make([]byte, 0, 84)
	{
		h := make([]byte, 48)
		le.PutUint32(h[0:], magicCIRTree)
		le.PutUint32(h[4:], 1)  // blockSize
		le.PutUint64(h[8:], 1)  // itemCount
		le.PutUint32(h[16:], 0) // startChromIx
		le.PutUint32(h[20:], spanStart)
		le.PutUint32(h[24:], 0) // endChromIx
		le.PutUint32(h[28:], spanEnd)
		le.PutUint64(h[32:], indexOff) // endFileOffset
		le.PutUint32(h[40:], 1)        // itemsPerSlot
		rtree = append(rtree, h...)
		node := make([]byte, 4)
		node[0] = 1 // isLeaf
		le.PutUint16(node[2:], 1)
		it := make([]byte, 32)
		le.PutUint32(it[0:], 0)
		le.PutUint32(it[4:], spanStart)
		le.PutUint32(it[8:], 0)
		le.PutUint32(it[12:], spanEnd)
		le.PutUint64(it[16:], dataOff)
		le.PutUint64(it[24:], uint64(len(dataBlock)))
		node = append(node, it...)
		rtree = append(rtree, node...)
	}

	hdr := make([]byte, 64)
	le.PutUint32(hdr[0:], magic)
	le.PutUint16(hdr[4:], 4) // version
	le.PutUint64(hdr[8:], chromTreeOff)
	le.PutUint64(hdr[16:], dataOff)
	le.PutUint64(hdr[24:], indexOff)
	// uncompressBufSize = 0 (data stored uncompressed)

	out := make([]byte, 0, indexOff+uint64(len(rtree)))
	out = append(out, hdr...)
	out = append(out, chromTree...)
	out = append(out, dataBlock...)
	out = append(out, rtree...)
	return out
}

func writeBigWig(t *testing.T, path, chrom string, chromSize uint32, items []wigItem) {
	t.Helper()
	le := binary.LittleEndian
	var minS, maxE uint32 = math.MaxUint32, 0
	for _, it := range items {
		if it.start < minS {
			minS = it.start
		}
		if it.end > maxE {
			maxE = it.end
		}
	}
	// one bedGraph section (type 1): 24-byte header + N×(start,end,val).
	blk := make([]byte, 24)
	le.PutUint32(blk[0:], 0) // chromId
	le.PutUint32(blk[4:], minS)
	le.PutUint32(blk[8:], maxE)
	blk[20] = 1 // bedGraph
	le.PutUint16(blk[22:], uint16(len(items)))
	for _, it := range items {
		row := make([]byte, 12)
		le.PutUint32(row[0:], it.start)
		le.PutUint32(row[4:], it.end)
		le.PutUint32(row[8:], math.Float32bits(it.val))
		blk = append(blk, row...)
	}
	if err := os.WriteFile(path, buildBBI(magicBigWig, chrom, chromSize, blk, minS, maxE), 0o644); err != nil {
		t.Fatal(err)
	}
}

func writeBigBed(t *testing.T, path, chrom string, chromSize uint32, items []bedItem) {
	t.Helper()
	le := binary.LittleEndian
	var minS, maxE uint32 = math.MaxUint32, 0
	var blk []byte
	for _, it := range items {
		if it.start < minS {
			minS = it.start
		}
		if it.end > maxE {
			maxE = it.end
		}
		row := make([]byte, 12)
		le.PutUint32(row[0:], 0) // chromId
		le.PutUint32(row[4:], it.start)
		le.PutUint32(row[8:], it.end)
		blk = append(blk, row...)
		blk = append(blk, []byte(it.rest)...)
		blk = append(blk, 0)
	}
	if err := os.WriteFile(path, buildBBI(magicBigBed, chrom, chromSize, blk, minS, maxE), 0o644); err != nil {
		t.Fatal(err)
	}
}

func collect(t *testing.T, r *Reader, ref string, start, end int) []*Record {
	t.Helper()
	seq, err := r.Query(ref, start, end)
	if err != nil {
		t.Fatalf("Query: %v", err)
	}
	var out []*Record
	for rec, err := range seq {
		if err != nil {
			t.Fatalf("iterate: %v", err)
		}
		out = append(out, rec)
	}
	return out
}

func TestBigWigQuery(t *testing.T) {
	path := filepath.Join(t.TempDir(), "s.bw")
	writeBigWig(t, path, "chr1", 1000, []wigItem{
		{100, 101, 0.5},
		{101, 102, 1.25},
		{200, 201, 9.0},
	})
	r, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	if r.Kind() != BigWig {
		t.Errorf("kind = %v, want BigWig", r.Kind())
	}
	if !r.HasRef("chr1") || r.HasRef("chr2") {
		t.Errorf("HasRef wrong: %v", r.RefNames())
	}
	// point query at base 100 (0-based [100,101)).
	got := collect(t, r, "chr1", 100, 101)
	if len(got) != 1 || got[0].Value != 0.5 {
		t.Fatalf("at 100: %+v", got)
	}
	// range covering the first two.
	got = collect(t, r, "chr1", 100, 102)
	if len(got) != 2 || got[1].Value != 1.25 {
		t.Fatalf("range: %+v", got)
	}
	// gap returns nothing.
	if got = collect(t, r, "chr1", 150, 160); len(got) != 0 {
		t.Fatalf("gap: %+v", got)
	}
}

func TestBigBedQuery(t *testing.T) {
	path := filepath.Join(t.TempDir(), "f.bb")
	writeBigBed(t, path, "chr1", 1000, []bedItem{
		{100, 150, "geneA\t42"},
		{300, 400, "geneB\t7"},
	})
	r, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	if r.Kind() != BigBed {
		t.Errorf("kind = %v, want BigBed", r.Kind())
	}
	got := collect(t, r, "chr1", 120, 130) // overlaps geneA only
	if len(got) != 1 || got[0].Line != "chr1\t100\t150\tgeneA\t42" {
		t.Fatalf("overlap: %+v", got)
	}
	if got = collect(t, r, "chr1", 200, 250); len(got) != 0 {
		t.Fatalf("gap: %+v", got)
	}
	got = collect(t, r, "chr1", 0, 1000) // both
	if len(got) != 2 {
		t.Fatalf("all: %+v", got)
	}
}
