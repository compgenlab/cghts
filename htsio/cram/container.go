package cram

import (
	"encoding/binary"
	"fmt"
	"io"
)

// fileDefinition is the 26-byte header at the start of a CRAM file.
type fileDefinition struct {
	Major  byte
	Minor  byte
	FileID [20]byte
}

func readFileDefinition(r io.Reader) (*fileDefinition, error) {
	var magic [4]byte
	if _, err := io.ReadFull(r, magic[:]); err != nil {
		return nil, fmt.Errorf("reading magic: %w", err)
	}
	if magic != [4]byte{'C', 'R', 'A', 'M'} {
		return nil, fmt.Errorf("invalid CRAM magic: %x", magic)
	}

	var ver [2]byte
	if _, err := io.ReadFull(r, ver[:]); err != nil {
		return nil, fmt.Errorf("reading version: %w", err)
	}

	fd := &fileDefinition{
		Major: ver[0],
		Minor: ver[1],
	}

	if _, err := io.ReadFull(r, fd.FileID[:]); err != nil {
		return nil, fmt.Errorf("reading file ID: %w", err)
	}

	if fd.Major != 3 {
		return nil, fmt.Errorf("unsupported CRAM version: %d.%d (only v3.x supported)", fd.Major, fd.Minor)
	}

	return fd, nil
}

// containerHeader holds metadata for a CRAM container.
type containerHeader struct {
	Length         int32
	RefSeqID       int32
	StartPos       int32
	AlignmentSpan  int32
	NumRecords     int32
	RecordCounter  int64
	Bases          int64
	NumBlocks      int32
	Landmarks      []int32
}

// isEOF returns true if this is the EOF container.
func (c *containerHeader) isEOF() bool {
	return c.Length == 15 && c.RefSeqID == -1 && c.NumRecords == 0 && c.NumBlocks == 1
}

func readContainerHeader(r io.Reader) (*containerHeader, error) {
	// Length: int32 little-endian
	var length int32
	if err := binary.Read(r, binary.LittleEndian, &length); err != nil {
		return nil, err
	}

	// Ref seq ID: itf8
	refSeqID, err := readITF8(r)
	if err != nil {
		return nil, fmt.Errorf("reading ref seq ID: %w", err)
	}

	// Start position: itf8
	startPos, err := readITF8(r)
	if err != nil {
		return nil, fmt.Errorf("reading start pos: %w", err)
	}

	// Alignment span: itf8
	alignmentSpan, err := readITF8(r)
	if err != nil {
		return nil, fmt.Errorf("reading alignment span: %w", err)
	}

	// Number of records: itf8
	numRecords, err := readITF8(r)
	if err != nil {
		return nil, fmt.Errorf("reading num records: %w", err)
	}

	// Record counter: ltf8
	recordCounter, err := readLTF8(r)
	if err != nil {
		return nil, fmt.Errorf("reading record counter: %w", err)
	}

	// Bases: ltf8
	bases, err := readLTF8(r)
	if err != nil {
		return nil, fmt.Errorf("reading bases: %w", err)
	}

	// Number of blocks: itf8
	numBlocks, err := readITF8(r)
	if err != nil {
		return nil, fmt.Errorf("reading num blocks: %w", err)
	}

	// Landmarks: array<itf8>
	landmarks, err := readITF8Array(r)
	if err != nil {
		return nil, fmt.Errorf("reading landmarks: %w", err)
	}

	// CRC32: 4 bytes
	var crc [4]byte
	if _, err := io.ReadFull(r, crc[:]); err != nil {
		return nil, fmt.Errorf("reading container CRC32: %w", err)
	}

	return &containerHeader{
		Length:        length,
		RefSeqID:      refSeqID,
		StartPos:      startPos,
		AlignmentSpan: alignmentSpan,
		NumRecords:    numRecords,
		RecordCounter: recordCounter,
		Bases:         bases,
		NumBlocks:     numBlocks,
		Landmarks:     landmarks,
	}, nil
}
