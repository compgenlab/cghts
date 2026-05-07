package sam

import (
	"bufio"
	"fmt"
	"io"
	"iter"
	"os/exec"
	"strconv"
	"strings"

	"github.com/compgen-io/cgltk/htsio"
)

// SamtoolsReader reads SAM/BAM/CRAM files by executing samtools view.
type SamtoolsReader struct {
	filename string
	opts     *htsio.SamReaderOpts
	header   *htsio.SamHeader
	cmd      *exec.Cmd
	stdout   io.ReadCloser
	stderr   io.ReadCloser
	scanner  *bufio.Scanner
	started  bool
	nextLine string
}

// NewSamtoolsReader creates a SamtoolsReader for the given file.
func NewSamtoolsReader(filename string, opts *htsio.SamReaderOpts) (htsio.SamReader, error) {
	if err := checkSamtools(); err != nil {
		return nil, err
	}
	if opts == nil {
		opts = htsio.NewSamReaderOpts()
	}
	return &SamtoolsReader{
		filename: filename,
		opts:     opts,
	}, nil
}

func checkSamtools() error {
	_, err := exec.LookPath("samtools")
	if err != nil {
		return fmt.Errorf("samtools not found in PATH: %w", err)
	}
	return nil
}

func (r *SamtoolsReader) start() error {
	if r.started {
		return nil
	}

	args := []string{"view", "-h"}
	if r.opts.ThreadsValue() > 0 {
		args = append(args, "--threads", strconv.Itoa(r.opts.ThreadsValue()))
	}
	if r.opts.FlagReqValue() != 0 {
		args = append(args, "-f", strconv.Itoa(r.opts.FlagReqValue()))
	}
	if r.opts.FlagFilterValue() != 0 {
		args = append(args, "-F", strconv.Itoa(r.opts.FlagFilterValue()))
	}
	if r.opts.MinMapQValue() != 0 {
		args = append(args, "-q", strconv.Itoa(r.opts.MinMapQValue()))
	}
	args = append(args, r.filename)

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
	r.scanner.Buffer(make([]byte, 0, 1024*1024), 10*1024*1024)
	r.started = true
	r.populateHeader()
	return nil
}

func (r *SamtoolsReader) Header() (*htsio.SamHeader, error) {
	if err := r.start(); err != nil {
		return nil, err
	}
	return r.header, nil
}

func (r *SamtoolsReader) populateHeader() {
	for r.scanner.Scan() {
		line := r.scanner.Text()
		if strings.HasPrefix(line, "@") {
			if r.header == nil {
				r.header = htsio.NewSamHeader()
			}
			r.header.AddLine(line)
			continue
		}
		r.nextLine = line
		return
	}
}

func (r *SamtoolsReader) Records() iter.Seq2[*htsio.SamRecord, error] {
	return func(yield func(*htsio.SamRecord, error) bool) {
		if err := r.start(); err != nil {
			yield(nil, err)
			return
		}

		if r.nextLine != "" {
			rec, err := parseSamLine(r.nextLine)
			r.nextLine = ""
			if err != nil {
				yield(nil, fmt.Errorf("parse SAM: %w", err))
				return
			}
			if r.opts.PassesFilters(rec) {
				if !yield(rec, nil) {
					return
				}
			}
		}

		for r.scanner.Scan() {
			line := r.scanner.Text()
			if strings.HasPrefix(line, "@") {
				if r.header == nil {
					r.header = htsio.NewSamHeader()
				}
				r.header.AddLine(line)
				continue
			}
			rec, err := parseSamLine(line)
			if err != nil {
				yield(nil, fmt.Errorf("parse SAM: %w", err))
				return
			}
			if r.opts.PassesFilters(rec) {
				if !yield(rec, nil) {
					return
				}
			}
		}

		if err := r.scanner.Err(); err != nil {
			yield(nil, fmt.Errorf("samtools read: %w", err))
			return
		}
	}
}

func (r *SamtoolsReader) Close() error {
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

func (r *SamtoolsReader) Query(ref string, start, end int) (iter.Seq2[*htsio.SamRecord, error], error) {
	region := fmt.Sprintf("%s:%d-%d", ref, start+1, end)

	opts := htsio.NewSamReaderOpts()
	if r.opts != nil {
		opts.FlagRequired(r.opts.FlagReqValue())
		opts.FlagFilter(r.opts.FlagFilterValue())
		opts.MinMapQ(r.opts.MinMapQValue())
		opts.Threads(r.opts.ThreadsValue())
	}

	sr := &SamtoolsReader{
		filename: r.filename,
		opts:     opts,
	}

	args := []string{"view", "-h"}
	if opts.ThreadsValue() > 0 {
		args = append(args, "--threads", strconv.Itoa(opts.ThreadsValue()))
	}
	if opts.FlagReqValue() != 0 {
		args = append(args, "-f", strconv.Itoa(opts.FlagReqValue()))
	}
	if opts.FlagFilterValue() != 0 {
		args = append(args, "-F", strconv.Itoa(opts.FlagFilterValue()))
	}
	if opts.MinMapQValue() != 0 {
		args = append(args, "-q", strconv.Itoa(opts.MinMapQValue()))
	}
	args = append(args, r.filename, region)

	sr.cmd = exec.Command("samtools", args...)

	var err error
	sr.stdout, err = sr.cmd.StdoutPipe()
	if err != nil {
		return nil, fmt.Errorf("samtools stdout pipe: %w", err)
	}
	sr.stderr, err = sr.cmd.StderrPipe()
	if err != nil {
		return nil, fmt.Errorf("samtools stderr pipe: %w", err)
	}
	if err := sr.cmd.Start(); err != nil {
		return nil, fmt.Errorf("samtools start: %w", err)
	}

	sr.scanner = bufio.NewScanner(sr.stdout)
	sr.scanner.Buffer(make([]byte, 0, 1024*1024), 10*1024*1024)
	sr.started = true
	sr.populateHeader()

	return func(yield func(*htsio.SamRecord, error) bool) {
		defer sr.Close()
		for rec, err := range sr.Records() {
			if !yield(rec, err) {
				return
			}
			if err != nil {
				return
			}
		}
	}, nil
}
