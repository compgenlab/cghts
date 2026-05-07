package cram

import (
	"bytes"
	"compress/bzip2"
	"compress/gzip"
	"fmt"
	"hash/crc32"
	"io"

	"github.com/compgen-io/cgkit/htsio/codec"
	"github.com/ulikunitz/xz/lzma"
)

// Block compression methods.
const (
	blockMethodRaw    = 0
	blockMethodGzip   = 1
	blockMethodBzip2  = 2
	blockMethodLzma   = 3
	blockMethodRans4x8  = 4
	blockMethodRans4x16 = 5
	blockMethodAdaptive = 6
	blockMethodFqzcomp  = 7
	blockMethodNameTok  = 8
)

// Block content types.
const (
	blockContentFileHeader       = 0
	blockContentCompressionHeader = 1
	blockContentSliceHeader      = 2
	blockContentReserved         = 3
	blockContentExternalData     = 4
	blockContentCoreData         = 5
)

// block represents a single CRAM block.
type block struct {
	method     byte
	contentType byte
	contentID  int32
	rawSize    int32
	data       []byte // decompressed data
}

// readBlock reads and decompresses a single CRAM block.
func readBlock(r io.Reader, majorVersion byte) (*block, error) {
	// v3+ uses CRC32 on blocks.
	var tr io.Reader
	var h hash32
	if majorVersion >= 3 {
		h = crc32.NewIEEE()
		tr = io.TeeReader(r, h)
	} else {
		tr = r
	}

	// Read method byte.
	var buf [1]byte
	if _, err := io.ReadFull(tr, buf[:]); err != nil {
		return nil, err
	}
	method := buf[0]

	// Content type.
	if _, err := io.ReadFull(tr, buf[:]); err != nil {
		return nil, fmt.Errorf("reading content type: %w", err)
	}
	contentType := buf[0]

	// Content ID.
	contentID, err := readITF8(tr)
	if err != nil {
		return nil, fmt.Errorf("reading content ID: %w", err)
	}

	// Compressed size.
	compSize, err := readITF8(tr)
	if err != nil {
		return nil, fmt.Errorf("reading compressed size: %w", err)
	}

	// Raw size.
	rawSize, err := readITF8(tr)
	if err != nil {
		return nil, fmt.Errorf("reading raw size: %w", err)
	}

	// Read compressed data.
	compData := make([]byte, compSize)
	if _, err := io.ReadFull(tr, compData); err != nil {
		return nil, fmt.Errorf("reading block data (%d bytes): %w", compSize, err)
	}

	// v3+ has CRC32 after the block data.
	if majorVersion >= 3 {
		var crcBuf [4]byte
		if _, err := io.ReadFull(r, crcBuf[:]); err != nil {
			return nil, fmt.Errorf("reading block CRC32: %w", err)
		}
		stored := uint32(crcBuf[0]) | uint32(crcBuf[1])<<8 | uint32(crcBuf[2])<<16 | uint32(crcBuf[3])<<24
		if computed := h.Sum32(); computed != stored {
			return nil, fmt.Errorf("block CRC32 mismatch: computed %08x, stored %08x", computed, stored)
		}
	}

	// Decompress.
	data, err := decompressBlock(method, compData, rawSize)
	if err != nil {
		return nil, fmt.Errorf("decompressing block (method=%d): %w", method, err)
	}

	return &block{
		method:      method,
		contentType: contentType,
		contentID:   contentID,
		rawSize:     rawSize,
		data:        data,
	}, nil
}

var methodNames = map[byte]string{
	blockMethodRaw:      "raw",
	blockMethodGzip:     "gzip",
	blockMethodBzip2:    "bzip2",
	blockMethodLzma:     "lzma",
	blockMethodRans4x8:  "rANS 4x8",
	blockMethodRans4x16: "rANS 4x16",
	blockMethodAdaptive: "adaptive arithmetic",
	blockMethodFqzcomp:  "fqzcomp",
	blockMethodNameTok:  "name tokenizer",
}

func decompressBlock(method byte, data []byte, rawSize int32) ([]byte, error) {
	switch method {
	case blockMethodRaw:
		return data, nil
	case blockMethodGzip:
		return decompressGzip(data)
	case blockMethodBzip2:
		return decompressBzip2(data)
	case blockMethodLzma:
		return decompressLzma(data)
	case blockMethodRans4x8:
		return codec.DecodeRans4x8(data)
	case blockMethodRans4x16:
		return codec.DecodeRansNx16(data)
	case blockMethodAdaptive:
		return codec.DecodeArithDynamic(data)
	case blockMethodFqzcomp:
		return codec.DecodeFqzcomp(data)
	case blockMethodNameTok:
		return codec.DecodeNameTokenizer(data)
	default:
		name, ok := methodNames[method]
		if !ok {
			name = fmt.Sprintf("unknown(%d)", method)
		}
		return nil, fmt.Errorf("unsupported compression method: %s", name)
	}
}

func decompressGzip(data []byte) ([]byte, error) {
	gr, err := gzip.NewReader(bytes.NewReader(data))
	if err != nil {
		return nil, fmt.Errorf("gzip reader: %w", err)
	}
	defer gr.Close()
	return io.ReadAll(gr)
}

func decompressBzip2(data []byte) ([]byte, error) {
	return io.ReadAll(bzip2.NewReader(bytes.NewReader(data)))
}

func decompressLzma(data []byte) ([]byte, error) {
	r, err := lzma.NewReader(bytes.NewReader(data))
	if err != nil {
		return nil, fmt.Errorf("lzma reader: %w", err)
	}
	return io.ReadAll(r)
}

