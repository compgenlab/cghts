package cram

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// craiEntry represents a single line from a CRAI index file.
type craiEntry struct {
	seqID           int   // reference sequence ID (0-based, -1=unmapped, -2=multi-ref)
	alignmentStart  int   // 1-based start position
	alignmentSpan   int   // span in bases
	containerOffset int64 // byte offset of the container in the CRAM file
	sliceOffset     int64 // byte offset of the slice within the container
	sliceSize       int64 // size of the slice in bytes
}

// craiIndex holds parsed CRAI index entries.
type craiIndex struct {
	entries []craiEntry
}

// loadCRAI reads and parses a gzip-compressed CRAI index file.
func loadCRAI(path string) (*craiIndex, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	gz, err := gzip.NewReader(f)
	if err != nil {
		return nil, fmt.Errorf("crai: gzip reader: %w", err)
	}
	defer gz.Close()

	var entries []craiEntry
	scanner := bufio.NewScanner(gz)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 6 {
			return nil, fmt.Errorf("crai: expected 6 fields, got %d: %q", len(fields), line)
		}

		seqID, err := strconv.Atoi(fields[0])
		if err != nil {
			return nil, fmt.Errorf("crai: parsing seqID: %w", err)
		}
		start, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("crai: parsing start: %w", err)
		}
		span, err := strconv.Atoi(fields[2])
		if err != nil {
			return nil, fmt.Errorf("crai: parsing span: %w", err)
		}
		contOff, err := strconv.ParseInt(fields[3], 10, 64)
		if err != nil {
			return nil, fmt.Errorf("crai: parsing containerOffset: %w", err)
		}
		sliceOff, err := strconv.ParseInt(fields[4], 10, 64)
		if err != nil {
			return nil, fmt.Errorf("crai: parsing sliceOffset: %w", err)
		}
		sliceSize, err := strconv.ParseInt(fields[5], 10, 64)
		if err != nil {
			return nil, fmt.Errorf("crai: parsing sliceSize: %w", err)
		}

		entries = append(entries, craiEntry{
			seqID:           seqID,
			alignmentStart:  start,
			alignmentSpan:   span,
			containerOffset: contOff,
			sliceOffset:     sliceOff,
			sliceSize:       sliceSize,
		})
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("crai: scanning: %w", err)
	}

	return &craiIndex{entries: entries}, nil
}

// query returns CRAI entries whose slices overlap [start, end) on the given reference.
// start and end are 0-based half-open coordinates.
func (idx *craiIndex) query(seqID int, start, end int) []craiEntry {
	var result []craiEntry
	for _, e := range idx.entries {
		if e.seqID != seqID {
			continue
		}
		// CRAI alignmentStart is 1-based; convert to 0-based for comparison.
		eStart := e.alignmentStart - 1
		eEnd := eStart + e.alignmentSpan
		// Check overlap: [eStart, eEnd) ∩ [start, end)
		if eStart < end && eEnd > start {
			result = append(result, e)
		}
	}
	return result
}
