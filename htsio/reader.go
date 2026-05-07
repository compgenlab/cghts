package htsio

import (
	"fmt"
	"io"
	"iter"
	"strconv"
	"strings"
)

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
	Val    string          // single value for eq/not-eq/contains/numeric ops
	Values map[string]bool // value set for TagInSet/TagNotInSet ops
}

// matchesRecord returns true if the SAM record passes this tag filter.
func (f *TagFilter) matchesRecord(rec *SamRecord) bool {
	t, ok := rec.Tags[f.Tag]
	if !ok {
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

// SamReaderOpts configures SAM/BAM/CRAM reader behavior.
type SamReaderOpts struct {
	flagReq    int
	flagFilter int
	minMapQ    int
	threads    int
	tagFilters []*TagFilter
}

// NewSamReader opens a SAM/BAM/CRAM file by auto-detecting the format
// from magic bytes.
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
// auto-detecting the format from magic bytes.
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

// NewSamReaderOpts returns default reader options.
func NewSamReaderOpts() *SamReaderOpts {
	return &SamReaderOpts{}
}

// PassesFilters returns true if the record passes all configured filters.
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

func (r *SamReaderOpts) FlagReqValue() int    { return r.flagReq }
func (r *SamReaderOpts) FlagFilterValue() int  { return r.flagFilter }
func (r *SamReaderOpts) MinMapQValue() int     { return r.minMapQ }
func (r *SamReaderOpts) ThreadsValue() int     { return r.threads }

func (r *SamReaderOpts) FlagRequired(flag int) *SamReaderOpts { r.flagReq = flag; return r }
func (r *SamReaderOpts) FlagFilter(flag int) *SamReaderOpts   { r.flagFilter = flag; return r }
func (r *SamReaderOpts) MinMapQ(mapq int) *SamReaderOpts      { r.minMapQ = mapq; return r }
func (r *SamReaderOpts) Threads(n int) *SamReaderOpts         { r.threads = n; return r }
func (r *SamReaderOpts) AddTagFilter(f *TagFilter) *SamReaderOpts {
	r.tagFilters = append(r.tagFilters, f)
	return r
}

// IterReader wraps an iter.Seq2[*SamRecord, error] as a SamReader.
func IterReader(seq iter.Seq2[*SamRecord, error], hdr *SamHeader) SamReader {
	return &iterReaderState{seq: seq, hdr: hdr}
}

type iterReaderState struct {
	seq iter.Seq2[*SamRecord, error]
	hdr *SamHeader
}

func (r *iterReaderState) Records() iter.Seq2[*SamRecord, error] { return r.seq }
func (r *iterReaderState) Header() (*SamHeader, error)           { return r.hdr, nil }
func (r *iterReaderState) Close() error                          { return nil }
func (r *iterReaderState) Query(ref string, start, end int) (iter.Seq2[*SamRecord, error], error) {
	return nil, fmt.Errorf("Query not supported on iterator reader")
}
