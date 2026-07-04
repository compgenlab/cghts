package annotate

import (
	"encoding/binary"
	"math"
	"os"
	"path/filepath"
	"testing"
)

// Minimal BBI writers for annotator tests: one chrom, one uncompressed data
// block, one R-tree leaf. (The reader has its own coverage in htsio/bbi; these
// exist so the annotators can be tested against real .bw/.bb files.)

const (
	tMagicBigWig  = 0x888FFC26
	tMagicBigBed  = 0x8789F2EB
	tMagicChromBP = 0x78CA8C91
	tMagicCIRTree = 0x2468ACE0
)

func buildTestBBI(magic uint32, chrom string, chromSize uint32, dataBlock []byte, spanStart, spanEnd uint32) []byte {
	le := binary.LittleEndian
	keySize := uint32(len(chrom))

	chromTree := make([]byte, 32)
	le.PutUint32(chromTree[0:], tMagicChromBP)
	le.PutUint32(chromTree[4:], 1)
	le.PutUint32(chromTree[8:], keySize)
	le.PutUint32(chromTree[12:], 8)
	le.PutUint64(chromTree[16:], 1)
	node := make([]byte, 4)
	node[0] = 1
	le.PutUint16(node[2:], 1)
	item := make([]byte, keySize+8)
	copy(item, chrom)
	le.PutUint32(item[keySize:], 0)
	le.PutUint32(item[keySize+4:], chromSize)
	chromTree = append(chromTree, append(node, item...)...)

	const chromTreeOff = 64
	dataOff := uint64(chromTreeOff + len(chromTree))
	indexOff := dataOff + uint64(len(dataBlock))

	rtree := make([]byte, 48)
	le.PutUint32(rtree[0:], tMagicCIRTree)
	le.PutUint32(rtree[4:], 1)
	le.PutUint64(rtree[8:], 1)
	le.PutUint32(rtree[20:], spanStart)
	le.PutUint32(rtree[28:], spanEnd)
	le.PutUint64(rtree[32:], indexOff)
	le.PutUint32(rtree[40:], 1)
	rn := make([]byte, 4)
	rn[0] = 1
	le.PutUint16(rn[2:], 1)
	rit := make([]byte, 32)
	le.PutUint32(rit[4:], spanStart)
	le.PutUint32(rit[12:], spanEnd)
	le.PutUint64(rit[16:], dataOff)
	le.PutUint64(rit[24:], uint64(len(dataBlock)))
	rtree = append(rtree, append(rn, rit...)...)

	hdr := make([]byte, 64)
	le.PutUint32(hdr[0:], magic)
	le.PutUint16(hdr[4:], 4)
	le.PutUint64(hdr[8:], chromTreeOff)
	le.PutUint64(hdr[16:], dataOff)
	le.PutUint64(hdr[24:], indexOff)

	out := append([]byte{}, hdr...)
	out = append(out, chromTree...)
	out = append(out, dataBlock...)
	out = append(out, rtree...)
	return out
}

type twig struct {
	start, end uint32
	val        float32
}

func writeTestBigWig(t *testing.T, chrom string, items []twig) string {
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
	blk := make([]byte, 24)
	le.PutUint32(blk[4:], minS)
	le.PutUint32(blk[8:], maxE)
	blk[20] = 1
	le.PutUint16(blk[22:], uint16(len(items)))
	for _, it := range items {
		row := make([]byte, 12)
		le.PutUint32(row[0:], it.start)
		le.PutUint32(row[4:], it.end)
		le.PutUint32(row[8:], math.Float32bits(it.val))
		blk = append(blk, row...)
	}
	path := filepath.Join(t.TempDir(), chrom+".bw")
	if err := os.WriteFile(path, buildTestBBI(tMagicBigWig, chrom, 100000, blk, minS, maxE), 0o644); err != nil {
		t.Fatal(err)
	}
	return path
}

type tbed struct {
	start, end uint32
	rest       string
}

func writeTestBigBed(t *testing.T, chrom string, items []tbed) string {
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
		le.PutUint32(row[4:], it.start)
		le.PutUint32(row[8:], it.end)
		blk = append(blk, row...)
		blk = append(blk, []byte(it.rest)...)
		blk = append(blk, 0)
	}
	path := filepath.Join(t.TempDir(), chrom+".bb")
	if err := os.WriteFile(path, buildTestBBI(tMagicBigBed, chrom, 100000, blk, minS, maxE), 0o644); err != nil {
		t.Fatal(err)
	}
	return path
}
