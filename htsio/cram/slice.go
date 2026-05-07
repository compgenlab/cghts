package cram

import (
	"fmt"
	"io"
)

// sliceHeader holds metadata for a CRAM slice.
type sliceHeader struct {
	refSeqID       int32
	alignmentStart int32
	alignmentSpan  int32
	numRecords     int32
	recordCounter  int64
	numBlocks      int32
	blockContentIDs []int32
	embeddedRefID  int32
	referenceMD5   [16]byte
}

func readSliceHeader(data []byte) (*sliceHeader, error) {
	r := newByteReader(data)

	refSeqID, err := readITF8(r)
	if err != nil {
		return nil, fmt.Errorf("reading ref seq ID: %w", err)
	}

	alignmentStart, err := readITF8(r)
	if err != nil {
		return nil, fmt.Errorf("reading alignment start: %w", err)
	}

	alignmentSpan, err := readITF8(r)
	if err != nil {
		return nil, fmt.Errorf("reading alignment span: %w", err)
	}

	numRecords, err := readITF8(r)
	if err != nil {
		return nil, fmt.Errorf("reading num records: %w", err)
	}

	recordCounter, err := readLTF8(r)
	if err != nil {
		return nil, fmt.Errorf("reading record counter: %w", err)
	}

	numBlocks, err := readITF8(r)
	if err != nil {
		return nil, fmt.Errorf("reading num blocks: %w", err)
	}

	blockContentIDs, err := readITF8Array(r)
	if err != nil {
		return nil, fmt.Errorf("reading block content IDs: %w", err)
	}

	embeddedRefID, err := readITF8(r)
	if err != nil {
		return nil, fmt.Errorf("reading embedded ref ID: %w", err)
	}

	var referenceMD5 [16]byte
	if _, err := io.ReadFull(r, referenceMD5[:]); err != nil {
		return nil, fmt.Errorf("reading reference MD5: %w", err)
	}

	// Optional tags may follow but we skip them.

	return &sliceHeader{
		refSeqID:        refSeqID,
		alignmentStart:  alignmentStart,
		alignmentSpan:   alignmentSpan,
		numRecords:      numRecords,
		recordCounter:   recordCounter,
		numBlocks:       numBlocks,
		blockContentIDs: blockContentIDs,
		embeddedRefID:   embeddedRefID,
		referenceMD5:    referenceMD5,
	}, nil
}

// byteReader wraps a byte slice to implement io.Reader.
type byteReader struct {
	data []byte
	pos  int
}

func newByteReader(data []byte) *byteReader {
	return &byteReader{data: data}
}

func (r *byteReader) Read(p []byte) (int, error) {
	if r.pos >= len(r.data) {
		return 0, io.EOF
	}
	n := copy(p, r.data[r.pos:])
	r.pos += n
	return n, nil
}
