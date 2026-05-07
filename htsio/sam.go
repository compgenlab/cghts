package htsio

import (
	"fmt"
	"iter"
	"strconv"
	"strings"
)

// SamTag represents a parsed optional field from a SAM record.
type SamTag struct {
	Type  byte // SAM type character: A, i, f, Z, H, B
	Value string
}

// IntValue returns the tag value as an int. Returns 0, false if the tag is not an integer type.
func (t SamTag) IntValue() (int, bool) {
	if t.Type != 'i' {
		return 0, false
	}
	v, err := strconv.Atoi(t.Value)
	if err != nil {
		return 0, false
	}
	return v, true
}

// StringValue returns the tag value as a string. Returns "", false if the tag is not a string type (Z or H).
func (t SamTag) StringValue() (string, bool) {
	if t.Type != 'Z' && t.Type != 'H' {
		return "", false
	}
	return t.Value, true
}

// FloatValue returns the tag value as a float64. Returns 0, false if the tag is not a float type.
func (t SamTag) FloatValue() (float64, bool) {
	if t.Type != 'f' {
		return 0, false
	}
	v, err := strconv.ParseFloat(t.Value, 64)
	if err != nil {
		return 0, false
	}
	return v, true
}

// SamRecord represents a single alignment record from a SAM/BAM/CRAM file.
type SamRecord struct {
	ReadName  string
	Flag      int
	RefName   string
	Pos       int // 1-based
	MapQ      int
	Cigar     string
	RefNext   string
	PosNext   int
	InsertLen int
	Seq       string
	Qual      string
	Tags      map[string]SamTag
	TagOrder  []string // insertion order of tag keys, for consistent output
}

// IsUnmapped returns true if the read is unmapped (flag 0x4).
func (r *SamRecord) IsUnmapped() bool {
	return r.Flag&0x4 != 0
}

// IsReverse returns true if the read is on the reverse strand (flag 0x10).
func (r *SamRecord) IsReverse() bool {
	return r.Flag&0x10 != 0
}

// IsSecondary returns true if the alignment is secondary (flag 0x100).
func (r *SamRecord) IsSecondary() bool {
	return r.Flag&0x100 != 0
}

// IsSupplementary returns true if the alignment is supplementary (flag 0x800).
func (r *SamRecord) IsSupplementary() bool {
	return r.Flag&0x800 != 0
}

// IsDuplicate returns true if the read is a PCR/optical duplicate (flag 0x400).
func (r *SamRecord) IsDuplicate() bool {
	return r.Flag&0x400 != 0
}

// String returns the record as a SAM-formatted line (no trailing newline).
// Tags are output in TagOrder if set, otherwise in map iteration order.
func (r *SamRecord) String() string {
	var sb strings.Builder
	fmt.Fprintf(&sb, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
		r.ReadName, r.Flag, r.RefName, r.Pos, r.MapQ,
		r.Cigar, r.RefNext, r.PosNext, r.InsertLen, r.Seq, r.Qual)
	if len(r.TagOrder) > 0 {
		for _, tag := range r.TagOrder {
			if val, ok := r.Tags[tag]; ok {
				fmt.Fprintf(&sb, "\t%s:%c:%s", tag, val.Type, val.Value)
			}
		}
	} else {
		for tag, val := range r.Tags {
			fmt.Fprintf(&sb, "\t%s:%c:%s", tag, val.Type, val.Value)
		}
	}
	return sb.String()
}

// SamReader is the interface for reading SAM/BAM/CRAM records.
type SamReader interface {
	// Records returns an iterator over all records in the reader.
	Records() iter.Seq2[*SamRecord, error]
	// Header returns the parsed SAM header. May return nil before the first Records() call.
	Header() (*SamHeader, error)
	// Query returns an iterator over records overlapping the 0-based
	// half-open region [start, end) on the given reference. Returns an
	// error if the file is not indexed.
	Query(ref string, start, end int) (iter.Seq2[*SamRecord, error], error)
	// Close releases resources.
	Close() error
}
