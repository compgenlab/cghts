package htsio

import (
	"fmt"
	"strconv"
	"strings"
)

// ParseRegion parses a samtools-style region string into 0-based half-open
// coordinates. Supported formats:
//   - "chr1"           → ref="chr1", start=0, end=-1 (whole chromosome)
//   - "chr1:1000-2000" → ref="chr1", start=999, end=2000 (1-based inclusive → 0-based half-open)
//   - "chr1:1000"      → ref="chr1", start=999, end=-1 (to end of chromosome)
//
// An end of -1 means "to end of reference."
func ParseRegion(region string) (ref string, start, end int, err error) {
	colonIdx := strings.Index(region, ":")
	if colonIdx < 0 {
		return region, 0, -1, nil
	}

	ref = region[:colonIdx]
	coords := region[colonIdx+1:]

	// Remove commas (samtools allows "1,000" style).
	coords = strings.ReplaceAll(coords, ",", "")

	dashIdx := strings.Index(coords, "-")
	if dashIdx < 0 {
		s, err := strconv.Atoi(coords)
		if err != nil {
			return "", 0, 0, fmt.Errorf("invalid region %q: %w", region, err)
		}
		return ref, s - 1, -1, nil
	}

	startStr := coords[:dashIdx]
	endStr := coords[dashIdx+1:]

	s, err := strconv.Atoi(startStr)
	if err != nil {
		return "", 0, 0, fmt.Errorf("invalid region start %q: %w", region, err)
	}
	e, err := strconv.Atoi(endStr)
	if err != nil {
		return "", 0, 0, fmt.Errorf("invalid region end %q: %w", region, err)
	}

	return ref, s - 1, e, nil
}
