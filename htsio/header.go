package htsio

import (
	"fmt"
	"strconv"
	"strings"
)

// SamHeader holds the parsed header from a SAM/BAM/CRAM file.
type SamHeader struct {
	Lines []string // raw header lines (including the @ prefix)
}

// SamHeaderRef represents a reference sequence from an @SQ header line.
type SamHeaderRef struct {
	Name   string
	Length int
}

// NewSamHeader creates an empty header.
func NewSamHeader() *SamHeader {
	return &SamHeader{}
}

// AddLine appends a raw header line.
func (h *SamHeader) AddLine(line string) {
	h.Lines = append(h.Lines, line)
}

// References returns the reference sequences parsed from @SQ lines.
func (h *SamHeader) References() []SamHeaderRef {
	var refs []SamHeaderRef
	for _, line := range h.Lines {
		if !strings.HasPrefix(line, "@SQ\t") {
			continue
		}
		ref := SamHeaderRef{}
		for _, field := range strings.Split(line, "\t")[1:] {
			if strings.HasPrefix(field, "SN:") {
				ref.Name = field[3:]
			} else if strings.HasPrefix(field, "LN:") {
				ref.Length, _ = strconv.Atoi(field[3:])
			}
		}
		if ref.Name != "" {
			refs = append(refs, ref)
		}
	}
	return refs
}

// ReadGroups returns the read group IDs from @RG lines.
func (h *SamHeader) ReadGroups() []string {
	var rgs []string
	for _, line := range h.Lines {
		if !strings.HasPrefix(line, "@RG\t") {
			continue
		}
		for _, field := range strings.Split(line, "\t")[1:] {
			if strings.HasPrefix(field, "ID:") {
				rgs = append(rgs, field[3:])
				break
			}
		}
	}
	return rgs
}

// Text returns the full header as a SAM-formatted string (with trailing newline).
func (h *SamHeader) Text() string {
	if len(h.Lines) == 0 {
		return ""
	}
	var sb strings.Builder
	for _, line := range h.Lines {
		fmt.Fprintln(&sb, line)
	}
	return sb.String()
}
