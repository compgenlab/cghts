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
	QName string
	Flag  int
	RName string
	Pos   int // 1-based
	MapQ  int
	Cigar string
	RNext string
	PNext int
	TLen  int
	Seq   string
	Qual  string
	Tags  map[string]SamTag
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
		r.QName, r.Flag, r.RName, r.Pos, r.MapQ,
		r.Cigar, r.RNext, r.PNext, r.TLen, r.Seq, r.Qual)
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
	Header() *SamHeader
	// Close releases resources.
	Close() error
}

// SamtoolsSamReader reads SAM/BAM/CRAM files by executing samtools view.
type SamtoolsSamReader struct {
	filename   string
	region     string
	flagReq    int
	flagFilter int
	minMapQ    int
	header     *SamHeader
	cmd        *exec.Cmd
	stdout     io.ReadCloser
	stderr     io.ReadCloser
	scanner    *bufio.Scanner
	started    bool
}

// NewSamReader creates a SamtoolsSamReader for the given file.
// Use the builder methods to set options before calling Next().
func NewSamReader(filename string) *SamtoolsSamReader {
	return &SamtoolsSamReader{
		filename: filename,
	}
}

// Region sets the genomic region to query (e.g. "chr1:1000-2000").
func (r *SamtoolsSamReader) Region(region string) *SamtoolsSamReader {
	r.region = region
	return r
}

// FlagRequired sets the required flags filter (-f). Only reads with all of these flags set are returned.
func (r *SamtoolsSamReader) FlagRequired(flag int) *SamtoolsSamReader {
	r.flagReq = flag
	return r
}

// FlagFilter sets the filtering flags (-F). Reads with any of these flags set are excluded.
func (r *SamtoolsSamReader) FlagFilter(flag int) *SamtoolsSamReader {
	r.flagFilter = flag
	return r
}

// MinMapQ sets the minimum mapping quality filter (-q).
func (r *SamtoolsSamReader) MinMapQ(mapq int) *SamtoolsSamReader {
	r.minMapQ = mapq
	return r
}

func (r *SamtoolsSamReader) start() error {
	if r.started {
		return nil
	}

	args := []string{"view", "-h"}
	if r.flagReq != 0 {
		args = append(args, "-f", strconv.Itoa(r.flagReq))
	}
	if r.flagFilter != 0 {
		args = append(args, "-F", strconv.Itoa(r.flagFilter))
	}
	if r.minMapQ != 0 {
		args = append(args, "-q", strconv.Itoa(r.minMapQ))
	}
	args = append(args, r.filename)
	if r.region != "" {
		args = append(args, r.region)
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
	return nil
}

// Header returns the parsed SAM header. The header is populated on the first
// call to Next(), so this will return nil before any records have been read.
func (r *SamtoolsSamReader) Header() *SamHeader {
	return r.header
}

// Next returns the next SamRecord. Returns nil, io.EOF when done.
func (r *SamtoolsSamReader) Next() (*SamRecord, error) {
	if err := r.start(); err != nil {
		return nil, err
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
		return rec, nil
	}

	if err := r.scanner.Err(); err != nil {
		return nil, fmt.Errorf("samtools read: %w", err)
	}

	return nil, io.EOF
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
		QName: fields[0],
		Flag:  flag,
		RName: fields[2],
		Pos:   pos,
		MapQ:  mapq,
		Cigar: fields[5],
		RNext: fields[6],
		PNext: pnext,
		TLen:  tlen,
		Seq:   fields[9],
		Qual:  fields[10],
		Tags:  make(map[string]SamTag),
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
