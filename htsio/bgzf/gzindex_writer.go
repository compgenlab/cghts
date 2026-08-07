package bgzf

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"encoding/binary"
	"fmt"
	"io"
	"os"
)

// BlockOffset is one BGZF block's position: where it begins in the compressed
// file, and what uncompressed offset its first byte carries.
type BlockOffset struct {
	Compressed   int64
	Uncompressed int64
}

// WalkBlocks decompresses a BGZF stream, returning the uncompressed bytes and
// the position of every block after the first.
//
// The boundaries are the point: a .gzi is exactly that list, and neither Reader
// nor IndexedReader exposes it — they consume boundaries rather than reporting
// them. The first block is omitted because (0,0) is where every reader already
// starts, which is also how LoadGZIndex reads the format back.
//
// The EOF marker — an empty terminating block — indexes nothing and is skipped;
// counting it yields a .gzi one entry longer than htslib writes.
func WalkBlocks(r io.Reader) ([]byte, []BlockOffset, error) {
	var (
		blocks  []BlockOffset
		out     bytes.Buffer
		coffset int64
		uoffset int64
		hdr     = make([]byte, 18)
		first   = true
	)
	br := bufio.NewReaderSize(r, 1<<20)
	for {
		n, err := io.ReadFull(br, hdr)
		if err == io.EOF || (err == io.ErrUnexpectedEOF && n == 0) {
			break
		}
		if err != nil {
			return nil, nil, fmt.Errorf("bgzf: read block header at %d: %w", coffset, err)
		}
		if hdr[0] != 0x1f || hdr[1] != 0x8b {
			return nil, nil, fmt.Errorf("bgzf: block at %d has no gzip magic", coffset)
		}
		if hdr[3]&0x04 == 0 {
			return nil, nil, fmt.Errorf("bgzf: plain gzip, not BGZF")
		}
		// The BC subfield carries BSIZE-1, the total block length, and is the
		// first extra subfield in every BGZF block — which is what the fixed
		// 18-byte header read above assumes.
		if hdr[12] != 'B' || hdr[13] != 'C' {
			return nil, nil, fmt.Errorf("bgzf: block at %d has no BC subfield", coffset)
		}
		bsize := int(binary.LittleEndian.Uint16(hdr[16:18])) + 1
		if bsize <= 18 {
			return nil, nil, fmt.Errorf("bgzf: block at %d declares length %d", coffset, bsize)
		}

		rest := make([]byte, bsize-18)
		if _, err := io.ReadFull(br, rest); err != nil {
			return nil, nil, fmt.Errorf("bgzf: read block at %d: %w", coffset, err)
		}
		whole := make([]byte, 0, bsize)
		whole = append(append(whole, hdr...), rest...)

		zr, err := gzip.NewReader(bytes.NewReader(whole))
		if err != nil {
			return nil, nil, fmt.Errorf("bgzf: block at %d: %w", coffset, err)
		}
		data, err := io.ReadAll(zr)
		zr.Close()
		if err != nil {
			return nil, nil, fmt.Errorf("bgzf: block at %d: %w", coffset, err)
		}
		if len(data) == 0 {
			break // EOF marker
		}
		if !first {
			blocks = append(blocks, BlockOffset{Compressed: coffset, Uncompressed: uoffset})
		}
		first = false

		out.Write(data)
		coffset += int64(bsize)
		uoffset += int64(len(data))
	}
	return out.Bytes(), blocks, nil
}

// WriteGZIndex writes a .gzi: uint64 count, then count pairs of (compressed,
// uncompressed) uint64, little-endian. The inverse of LoadGZIndex.
func WriteGZIndex(filename string, blocks []BlockOffset) error {
	buf := make([]byte, 8+16*len(blocks))
	binary.LittleEndian.PutUint64(buf[:8], uint64(len(blocks)))
	for i, b := range blocks {
		off := 8 + 16*i
		binary.LittleEndian.PutUint64(buf[off:off+8], uint64(b.Compressed))
		binary.LittleEndian.PutUint64(buf[off+8:off+16], uint64(b.Uncompressed))
	}
	return os.WriteFile(filename, buf, 0o644)
}
