package cram

import (
	"encoding/binary"
	"fmt"
	"hash/crc32"
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

	if fd.Major != 2 && fd.Major != 3 {
		return nil, fmt.Errorf("unsupported CRAM version: %d.%d (only v2.x and v3.x supported)", fd.Major, fd.Minor)
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
// v3: length=15, refSeqID=-1, numRecords=0, numBlocks=1
// v2: length=11, refSeqID=-1, numRecords=0, numBlocks=1
func (c *containerHeader) isEOF() bool {
	return (c.Length == 15 || c.Length == 11) && c.RefSeqID == -1 && c.NumRecords == 0 && c.NumBlocks == 1
}

func readContainerHeader(r io.Reader, majorVersion byte) (*containerHeader, error) {
	// v3+ uses CRC32 on container headers.
	var tr io.Reader
	var h hash32
	if majorVersion >= 3 {
		h = crc32.NewIEEE()
		tr = io.TeeReader(r, h)
	} else {
		tr = r
	}

	// Length: int32 little-endian
	var length int32
	if err := binary.Read(tr, binary.LittleEndian, &length); err != nil {
		return nil, err
	}

	// Ref seq ID: itf8
	refSeqID, err := readITF8(tr)
	if err != nil {
		return nil, fmt.Errorf("reading ref seq ID: %w", err)
	}

	// Start position: itf8
	startPos, err := readITF8(tr)
	if err != nil {
		return nil, fmt.Errorf("reading start pos: %w", err)
	}

	// Alignment span: itf8
	alignmentSpan, err := readITF8(tr)
	if err != nil {
		return nil, fmt.Errorf("reading alignment span: %w", err)
	}

	// Number of records: itf8
	numRecords, err := readITF8(tr)
	if err != nil {
		return nil, fmt.Errorf("reading num records: %w", err)
	}

	// Record counter: itf8 in v2, ltf8 in v3+
	var recordCounter int64
	if majorVersion >= 3 {
		recordCounter, err = readLTF8(tr)
	} else {
		var rc int32
		rc, err = readITF8(tr)
		recordCounter = int64(rc)
	}
	if err != nil {
		return nil, fmt.Errorf("reading record counter: %w", err)
	}

	// Bases: ltf8
	bases, err := readLTF8(tr)
	if err != nil {
		return nil, fmt.Errorf("reading bases: %w", err)
	}

	// Number of blocks: itf8
	numBlocks, err := readITF8(tr)
	if err != nil {
		return nil, fmt.Errorf("reading num blocks: %w", err)
	}

	// Landmarks: array<itf8>
	landmarks, err := readITF8Array(tr)
	if err != nil {
		return nil, fmt.Errorf("reading landmarks: %w", err)
	}

	// v3+ has CRC32 after the header.
	if majorVersion >= 3 {
		var crc [4]byte
		if _, err := io.ReadFull(r, crc[:]); err != nil {
			return nil, fmt.Errorf("reading container CRC32: %w", err)
		}
		stored := uint32(crc[0]) | uint32(crc[1])<<8 | uint32(crc[2])<<16 | uint32(crc[3])<<24
		if computed := h.Sum32(); computed != stored {
			return nil, fmt.Errorf("container header CRC32 mismatch: computed %08x, stored %08x", computed, stored)
		}
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

// hash32 is the interface satisfied by crc32.Hash.
type hash32 interface {
	io.Writer
	Sum32() uint32
}
