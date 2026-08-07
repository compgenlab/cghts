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

// IndexingReader decompresses a BGZF stream and records block boundaries as it
// goes.
//
// A reader rather than a function returning bytes, because the whole point is
// indexing files too large to hold: a human genome is ~3 GB uncompressed, and
// buffering that to find block boundaries trades a disk read for an OOM kill.
// Memory here is one block at a time.
//
// Blocks reports the position of every block after the first — the first is
// implicit at (0,0), which is where every reader already starts and how
// LoadGZIndex reads the format back. The EOF marker, an empty terminating
// block, indexes nothing and is skipped; counting it yields a .gzi one entry
// longer than htslib writes.
type IndexingReader struct {
	br      *bufio.Reader
	blocks  []BlockOffset
	buf     []byte // current block's uncompressed bytes, unread portion
	coffset int64
	uoffset int64
	first   bool
	done    bool
	err     error
}

// NewIndexingReader wraps a BGZF stream.
func NewIndexingReader(r io.Reader) *IndexingReader {
	return &IndexingReader{br: bufio.NewReaderSize(r, 1<<20), first: true}
}

// Blocks returns the boundaries seen so far; complete once Read has returned
// io.EOF.
func (r *IndexingReader) Blocks() []BlockOffset { return r.blocks }

func (r *IndexingReader) Read(p []byte) (int, error) {
	for len(r.buf) == 0 {
		if r.done {
			if r.err != nil {
				return 0, r.err
			}
			return 0, io.EOF
		}
		if err := r.next(); err != nil {
			r.done = true
			if err != io.EOF {
				r.err = err
				return 0, err
			}
			return 0, io.EOF
		}
	}
	n := copy(p, r.buf)
	r.buf = r.buf[n:]
	return n, nil
}

// next decodes one block into buf, advancing the offsets.
func (r *IndexingReader) next() error {
	hdr := make([]byte, 18)
	n, err := io.ReadFull(r.br, hdr)
	if err == io.EOF || (err == io.ErrUnexpectedEOF && n == 0) {
		return io.EOF
	}
	if err != nil {
		return fmt.Errorf("bgzf: read block header at %d: %w", r.coffset, err)
	}
	if hdr[0] != 0x1f || hdr[1] != 0x8b {
		return fmt.Errorf("bgzf: block at %d has no gzip magic", r.coffset)
	}
	if hdr[3]&0x04 == 0 {
		return fmt.Errorf("bgzf: plain gzip, not BGZF")
	}
	// The BC subfield carries BSIZE-1, the total block length, and is the first
	// extra subfield in every BGZF block — which is what the fixed 18-byte
	// header read above assumes.
	if hdr[12] != 'B' || hdr[13] != 'C' {
		return fmt.Errorf("bgzf: block at %d has no BC subfield", r.coffset)
	}
	bsize := int(binary.LittleEndian.Uint16(hdr[16:18])) + 1
	if bsize <= 18 {
		return fmt.Errorf("bgzf: block at %d declares length %d", r.coffset, bsize)
	}

	rest := make([]byte, bsize-18)
	if _, err := io.ReadFull(r.br, rest); err != nil {
		return fmt.Errorf("bgzf: read block at %d: %w", r.coffset, err)
	}
	whole := make([]byte, 0, bsize)
	whole = append(append(whole, hdr...), rest...)

	zr, err := gzip.NewReader(bytes.NewReader(whole))
	if err != nil {
		return fmt.Errorf("bgzf: block at %d: %w", r.coffset, err)
	}
	data, err := io.ReadAll(zr)
	zr.Close()
	if err != nil {
		return fmt.Errorf("bgzf: block at %d: %w", r.coffset, err)
	}
	if len(data) == 0 {
		return io.EOF // the EOF marker
	}
	if !r.first {
		r.blocks = append(r.blocks, BlockOffset{Compressed: r.coffset, Uncompressed: r.uoffset})
	}
	r.first = false

	r.buf = data
	r.coffset += int64(bsize)
	r.uoffset += int64(len(data))
	return nil
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
