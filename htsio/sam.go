package htsio

import (
	"fmt"
	"io"
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

// TagFilterOp specifies the comparison operation for a tag filter.
type TagFilterOp int

const (
	TagEq          TagFilterOp = iota // equals
	TagNotEq                          // not equals
	TagContains                       // substring match
	TagNotContains                    // no substring match
	TagLt                             // less than (numeric)
	TagGt                             // greater than (numeric)
	TagLte                            // less than or equal (numeric)
	TagGte                            // greater than or equal (numeric)
	TagInSet                          // value is in a set
	TagNotInSet                       // value is not in a set
)

// TagFilter represents a single tag-based filter condition.
type TagFilter struct {
	Tag    string
	Op     TagFilterOp
	Val    string            // single value for eq/not-eq/contains/numeric ops
	Values map[string]bool   // value set for TagInSet/TagNotInSet ops
}

// matchesRecord returns true if the SAM record passes this tag filter.
func (f *TagFilter) matchesRecord(rec *SamRecord) bool {
	t, ok := rec.Tags[f.Tag]
	if !ok {
		// Missing tag: treat as empty string for equality checks.
		switch f.Op {
		case TagEq:
			return f.Val == ""
		case TagNotEq:
			return f.Val != ""
		default:
			return false
		}
	}

	switch f.Op {
	case TagEq:
		return t.Value == f.Val
	case TagNotEq:
		return t.Value != f.Val
	case TagContains:
		return strings.Contains(t.Value, f.Val)
	case TagNotContains:
		return !strings.Contains(t.Value, f.Val)
	case TagLt, TagGt, TagLte, TagGte:
		return f.numericCompare(t)
	case TagInSet:
		return f.Values[t.Value]
	case TagNotInSet:
		return !f.Values[t.Value]
	}
	return false
}

func (f *TagFilter) numericCompare(t SamTag) bool {
	switch t.Type {
	case 'i':
		tv, ok := t.IntValue()
		if !ok {
			return false
		}
		fv, err := strconv.Atoi(f.Val)
		if err != nil {
			return false
		}
		switch f.Op {
		case TagLt:
			return tv < fv
		case TagGt:
			return tv > fv
		case TagLte:
			return tv <= fv
		case TagGte:
			return tv >= fv
		}
	case 'f':
		tv, ok := t.FloatValue()
		if !ok {
			return false
		}
		fv, err := strconv.ParseFloat(f.Val, 64)
		if err != nil {
			return false
		}
		switch f.Op {
		case TagLt:
			return tv < fv
		case TagGt:
			return tv > fv
		case TagLte:
			return tv <= fv
		case TagGte:
			return tv >= fv
		}
	}
	return false
}

// ParseTagFilter parses a "TAG:VALUE" string into a TagFilter with the given op.
func ParseTagFilter(s string, op TagFilterOp) (*TagFilter, error) {
	idx := strings.Index(s, ":")
	if idx < 1 {
		return nil, fmt.Errorf("invalid tag filter %q: expected TAG:VALUE", s)
	}
	return &TagFilter{
		Tag: s[:idx],
		Op:  op,
		Val: s[idx+1:],
	}, nil
}

type SamReaderOpts struct {
	flagReq    int
	flagFilter int
	minMapQ    int
	threads    int
	tagFilters []*TagFilter
}

// NewSamReader opens a SAM/BAM/CRAM file by auto-detecting the format
// from magic bytes. The file is peeked to detect format, then the matched
// reader opens its own file handle(s) directly — no nested readers.
func NewSamReader(filename string, opts ...*SamReaderOpts) (SamReader, error) {
	var o *SamReaderOpts
	if len(opts) > 0 {
		o = opts[0]
	} else {
		o = NewSamReaderOpts()
	}

	reg, err := detectFromFile(filename)
	if err != nil {
		return nil, err
	}
	return reg.NewFromFile(filename, o)
}

// NewSamReaderFromReader creates a SamReader from an io.ReadCloser by
// auto-detecting the format from magic bytes. This is the stream path
// (e.g., stdin) — peeked bytes are prepended back via MultiReader.
// Query() is not supported on stream-based readers.
func NewSamReaderFromReader(r io.ReadCloser, opts ...*SamReaderOpts) (SamReader, error) {
	var o *SamReaderOpts
	if len(opts) > 0 {
		o = opts[0]
	} else {
		o = NewSamReaderOpts()
	}

	reg, fullReader, err := detectFromStream(r)
	if err != nil {
		return nil, err
	}
	return reg.NewFromStream(fullReader, o)
}

func NewSamReaderOpts() *SamReaderOpts {
	return &SamReaderOpts{}
}

// PassesFilters returns true if the record passes all configured filters
// (flag required, flag filter, min mapq, and tag filters).
func (r *SamReaderOpts) PassesFilters(rec *SamRecord) bool {
	if r.flagReq != 0 && rec.Flag&r.flagReq != r.flagReq {
		return false
	}
	if r.flagFilter != 0 && rec.Flag&r.flagFilter != 0 {
		return false
	}
	if r.minMapQ != 0 && rec.MapQ < r.minMapQ {
		return false
	}
	for _, f := range r.tagFilters {
		if !f.matchesRecord(rec) {
			return false
		}
	}
	return true
}

// FlagReqValue returns the required flags value.
func (r *SamReaderOpts) FlagReqValue() int { return r.flagReq }

// FlagFilterValue returns the filter flags value.
func (r *SamReaderOpts) FlagFilterValue() int { return r.flagFilter }

// MinMapQValue returns the minimum mapping quality value.
func (r *SamReaderOpts) MinMapQValue() int { return r.minMapQ }

// ThreadsValue returns the number of threads.
func (r *SamReaderOpts) ThreadsValue() int { return r.threads }

// FlagRequired sets the required flags filter (-f). Only reads with all of these flags set are returned.
func (r *SamReaderOpts) FlagRequired(flag int) *SamReaderOpts {
	r.flagReq = flag
	return r
}

// FlagFilter sets the filtering flags (-F). Reads with any of these flags set are excluded.
func (r *SamReaderOpts) FlagFilter(flag int) *SamReaderOpts {
	r.flagFilter = flag
	return r
}

// MinMapQ sets the minimum mapping quality filter (-q).
func (r *SamReaderOpts) MinMapQ(mapq int) *SamReaderOpts {
	r.minMapQ = mapq
	return r
}

// Threads sets the number of samtools decompression threads (--threads).
func (r *SamReaderOpts) Threads(n int) *SamReaderOpts {
	r.threads = n
	return r
}

// AddTagFilter adds a tag-based filter. Multiple filters are ANDed together.
func (r *SamReaderOpts) AddTagFilter(f *TagFilter) *SamReaderOpts {
	r.tagFilters = append(r.tagFilters, f)
	return r
}


// CigarRefLen returns the number of reference bases consumed by a CIGAR string.
// Operations M, D, N, =, X consume reference; I, S, H, P do not.
func CigarRefLen(cigar string) int {
	if cigar == "*" {
		return 0
	}
	refLen := 0
	num := 0
	for i := 0; i < len(cigar); i++ {
		c := cigar[i]
		if c >= '0' && c <= '9' {
			num = num*10 + int(c-'0')
		} else {
			switch c {
			case 'M', 'D', 'N', '=', 'X':
				refLen += num
			}
			num = 0
		}
	}
	return refLen
}

// ParseRegion parses a samtools-style region string into 0-based half-open
// coordinates. Supported formats:
//   - "chr1"           → ref="chr1", start=0, end=-1 (whole chromosome; caller should use ref length)
//   - "chr1:1000-2000" → ref="chr1", start=999, end=2000 (input is 1-based inclusive)
//   - "chr1:1000"      → ref="chr1", start=999, end=-1 (to end of chromosome)
//
// Returns ref, start, end. An end of -1 means "to end of reference."
func ParseRegion(region string) (ref string, start, end int, err error) {
	colonIdx := strings.Index(region, ":")
	if colonIdx < 0 {
		return region, 0, -1, nil
	}

	ref = region[:colonIdx]
	coords := region[colonIdx+1:]

	// Remove commas (samtools allows "1,000" style)
	coords = strings.ReplaceAll(coords, ",", "")

	dashIdx := strings.Index(coords, "-")
	if dashIdx < 0 {
		// Just start position.
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

	// samtools regions are 1-based inclusive → 0-based half-open
	return ref, s - 1, e, nil
}

// iterReaderState wraps an iter.Seq2 of SamRecords into a SamReader,
// enabling code that expects Records()/Header()/Close() to consume an
// iterator-based query result.
type iterReaderState struct {
	seq  iter.Seq2[*SamRecord, error]
	hdr  *SamHeader
}

// IterReader wraps an iter.Seq2[*SamRecord, error] as a SamReader.
func IterReader(seq iter.Seq2[*SamRecord, error], hdr *SamHeader) SamReader {
	return &iterReaderState{seq: seq, hdr: hdr}
}

func (r *iterReaderState) Records() iter.Seq2[*SamRecord, error] {
	return r.seq
}

func (r *iterReaderState) Header() (*SamHeader, error) {
	return r.hdr, nil
}

func (r *iterReaderState) Query(ref string, start, end int) (iter.Seq2[*SamRecord, error], error) {
	return nil, fmt.Errorf("Query not supported on iterator reader")
}

func (r *iterReaderState) Close() error {
	return nil
}
