package htsio

import (
	"fmt"
	"os"
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

// AddPGLine appends a @PG header line with PN:{pn} and the current command line.
// The id is used for the ID field. If a @PG line with the same ID already exists
// in the header, a numeric suffix (.1, .2, etc.) is appended to make it unique.
// Extra tab-delimited fields (e.g. "DS:...") can be appended via extras.
func (h *SamHeader) AddPGLine(id string, pn string, extras ...string) {
	uniqueID := h.uniquePGID(id)
	cl := strings.Join(os.Args, " ")
	var line strings.Builder
	line.WriteString(fmt.Sprintf("@PG\tID:%s\tPN:%s\tCL:%s", uniqueID, pn, cl))
	for _, extra := range extras {
		line.WriteString("\t" + extra)
	}
	h.AddLine(line.String())
}

// uniquePGID returns id if no @PG line uses it, otherwise id.1, id.2, etc.
func (h *SamHeader) uniquePGID(id string) string {
	existing := make(map[string]bool)
	for _, line := range h.Lines {
		if !strings.HasPrefix(line, "@PG\t") {
			continue
		}
		for _, field := range strings.Split(line, "\t")[1:] {
			if strings.HasPrefix(field, "ID:") {
				existing[field[3:]] = true
				break
			}
		}
	}
	if !existing[id] {
		return id
	}
	for i := 1; ; i++ {
		candidate := fmt.Sprintf("%s.%d", id, i)
		if !existing[candidate] {
			return candidate
		}
	}
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
