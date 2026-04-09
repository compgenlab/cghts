package htsio

import (
	"bufio"
	"fmt"
	"io"
	"os/exec"
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
func (r *SamRecord) String() string {
	var sb strings.Builder
	fmt.Fprintf(&sb, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
		r.ReadName, r.Flag, r.RefName, r.Pos, r.MapQ,
		r.Cigar, r.RefNext, r.PosNext, r.InsertLen, r.Seq, r.Qual)
	for tag, val := range r.Tags {
		fmt.Fprintf(&sb, "\t%s:%c:%s", tag, val.Type, val.Value)
	}
	return sb.String()
}

// SamReader is the interface for reading SAM/BAM/CRAM records.
type SamReader interface {
	// Next returns the next SamRecord. Returns nil, io.EOF when done.
	Next() (*SamRecord, error)
	// Header returns the parsed SAM header. May return nil before the first Next() call.
	Header() (*SamHeader, error)
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
		return false
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
	region     string
	flagReq    int
	flagFilter int
	minMapQ    int
	threads    int
	tagFilters []*TagFilter
}

// SamtoolsSamReader reads SAM/BAM/CRAM files by executing samtools view.
type SamtoolsSamReader struct {
	filename string
	opts     *SamReaderOpts
	header   *SamHeader
	cmd      *exec.Cmd
	stdout   io.ReadCloser
	stderr   io.ReadCloser
	scanner  *bufio.Scanner
	started  bool
	nextLine string
}

func NewSamReader(filename string, opts ...*SamReaderOpts) (SamReader, error) {
	return newSamtoolsReader(filename, opts...)
}

// NewSamtoolsReader creates a SamtoolsSamReader for the given file.
// Returns an error if samtools is not found in PATH.
// Use the builder methods to set options before calling Next().
func newSamtoolsReader(filename string, opts ...*SamReaderOpts) (SamReader, error) {
	if err := checkSamtools(); err != nil {
		return nil, err
	}
	if len(opts) == 0 {
		opts = []*SamReaderOpts{NewSamReaderOpts()}
	}
	return &SamtoolsSamReader{
		filename: filename,
		opts:     opts[0],
	}, nil
}

func NewSamReaderOpts() *SamReaderOpts {
	return &SamReaderOpts{}
}

// Region sets the genomic region to query (e.g. "chr1:1000-2000").
func (r *SamReaderOpts) Region(region string) *SamReaderOpts {
	r.region = region
	return r
}

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

func checkSamtools() error {
	_, err := exec.LookPath("samtools")
	if err != nil {
		return fmt.Errorf("samtools not found in PATH: %w", err)
	}
	return nil
}

func (r *SamtoolsSamReader) start() error {
	if r.started {
		return nil
	}

	args := []string{"view", "-h"}
	if r.opts.threads > 0 {
		args = append(args, "--threads", strconv.Itoa(r.opts.threads))
	}
	if r.opts.flagReq != 0 {
		args = append(args, "-f", strconv.Itoa(r.opts.flagReq))
	}
	if r.opts.flagFilter != 0 {
		args = append(args, "-F", strconv.Itoa(r.opts.flagFilter))
	}
	if r.opts.minMapQ != 0 {
		args = append(args, "-q", strconv.Itoa(r.opts.minMapQ))
	}
	args = append(args, r.filename)
	if r.opts.region != "" {
		args = append(args, r.opts.region)
	}

	r.cmd = exec.Command("samtools", args...)

	var err error
	r.stdout, err = r.cmd.StdoutPipe()
	if err != nil {
		return fmt.Errorf("samtools stdout pipe: %w", err)
	}

	r.stderr, err = r.cmd.StderrPipe()
	if err != nil {
		return fmt.Errorf("samtools stderr pipe: %w", err)
	}

	if err := r.cmd.Start(); err != nil {
		return fmt.Errorf("samtools start: %w", err)
	}

	r.scanner = bufio.NewScanner(r.stdout)
	r.scanner.Buffer(make([]byte, 0, 1024*1024), 10*1024*1024) // up to 10MB lines
	r.started = true
	r.populateHeader()
	return nil
}

// Header returns the parsed SAM header. The header is populated on the first
// call to Next(), so this will return nil before any records have been read.
func (r *SamtoolsSamReader) Header() (*SamHeader, error) {
	if err := r.start(); err != nil {
		return nil, err
	}
	return r.header, nil
}

func (r *SamtoolsSamReader) populateHeader() {
	for r.scanner.Scan() {
		line := r.scanner.Text()
		if strings.HasPrefix(line, "@") {
			if r.header == nil {
				r.header = NewSamHeader()
			}
			r.header.AddLine(line)
			continue
		}
		r.nextLine = line
		return
	}
}

// Next returns the next SamRecord. Returns nil, io.EOF when done.
func (r *SamtoolsSamReader) Next() (*SamRecord, error) {
	if err := r.start(); err != nil {
		return nil, err
	}

	if r.nextLine != "" {
		rec, err := parseSamLine(r.nextLine)
		r.nextLine = ""
		if err != nil {
			return nil, fmt.Errorf("parse SAM: %w", err)
		}
		if r.passesTagFilters(rec) {
			return rec, nil
		}
	}

	for r.scanner.Scan() {
		line := r.scanner.Text()
		if strings.HasPrefix(line, "@") {
			if r.header == nil {
				r.header = NewSamHeader()
			}
			r.header.AddLine(line)
			continue
		}
		rec, err := parseSamLine(line)
		if err != nil {
			return nil, fmt.Errorf("parse SAM: %w", err)
		}
		if r.passesTagFilters(rec) {
			return rec, nil
		}
	}

	if err := r.scanner.Err(); err != nil {
		return nil, fmt.Errorf("samtools read: %w", err)
	}

	return nil, io.EOF
}

// passesTagFilters returns true if the record passes all tag filters.
func (r *SamtoolsSamReader) passesTagFilters(rec *SamRecord) bool {
	for _, f := range r.opts.tagFilters {
		if !f.matchesRecord(rec) {
			return false
		}
	}
	return true
}

// Close waits for the samtools process to finish and releases resources.
func (r *SamtoolsSamReader) Close() error {
	if !r.started {
		return nil
	}
	if r.stdout != nil {
		r.stdout.Close()
	}
	if r.stderr != nil {
		io.ReadAll(r.stderr)
		r.stderr.Close()
	}
	if r.cmd != nil {
		return r.cmd.Wait()
	}
	return nil
}

func parseSamLine(line string) (*SamRecord, error) {
	fields := strings.SplitN(line, "\t", 12)
	if len(fields) < 11 {
		return nil, fmt.Errorf("expected at least 11 fields, got %d", len(fields))
	}

	flag, err := strconv.Atoi(fields[1])
	if err != nil {
		return nil, fmt.Errorf("invalid flag %q: %w", fields[1], err)
	}

	pos, err := strconv.Atoi(fields[3])
	if err != nil {
		return nil, fmt.Errorf("invalid pos %q: %w", fields[3], err)
	}

	mapq, err := strconv.Atoi(fields[4])
	if err != nil {
		return nil, fmt.Errorf("invalid mapq %q: %w", fields[4], err)
	}

	pnext, err := strconv.Atoi(fields[7])
	if err != nil {
		return nil, fmt.Errorf("invalid pnext %q: %w", fields[7], err)
	}

	tlen, err := strconv.Atoi(fields[8])
	if err != nil {
		return nil, fmt.Errorf("invalid tlen %q: %w", fields[8], err)
	}

	rec := &SamRecord{
		ReadName:  fields[0],
		Flag:      flag,
		RefName:   fields[2],
		Pos:       pos,
		MapQ:      mapq,
		Cigar:     fields[5],
		RefNext:   fields[6],
		PosNext:   pnext,
		InsertLen: tlen,
		Seq:       fields[9],
		Qual:      fields[10],
		Tags:      make(map[string]SamTag),
	}

	if len(fields) > 11 {
		for _, raw := range strings.Split(fields[11], "\t") {
			parts := strings.SplitN(raw, ":", 3)
			if len(parts) != 3 {
				continue
			}
			rec.Tags[parts[0]] = SamTag{
				Type:  parts[1][0],
				Value: parts[2],
			}
		}
	}

	return rec, nil
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
