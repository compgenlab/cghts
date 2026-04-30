package htsio

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"
	"os/exec"
	"strings"
	"sync"
)

// SamOutputFormat specifies the output file format for SamWriter.
type SamOutputFormat int

const (
	FormatSAM  SamOutputFormat = iota // plain text SAM
	FormatBAM                         // compressed BAM
	FormatCRAM                        // CRAM (requires reference)
)

// SamWriter is the interface for writing SAM/BAM/CRAM records.
type SamWriter interface {
	// Write writes a SamRecord to the output.
	Write(rec *SamRecord) error
	// Close flushes remaining data and releases resources.
	Close() error
}

type samWriterOptions struct {
	header        *SamHeader
	format        SamOutputFormat
	threads       int
	reference     string
	sortedCoord   bool
	sortedName    bool
	sortTmpPrefix string
}

// SamtoolsSamWriter writes SAM/BAM/CRAM files by piping SAM text to samtools view.
// It is safe for concurrent use from multiple goroutines.
type SamtoolsSamWriter struct {
	filename string
	opts     *samWriterOptions
	cmd      *exec.Cmd
	outs     io.WriteCloser
	stderr   io.ReadCloser
	started  bool
	writeCh  chan string    // buffered channel for async writes
	writeWg  sync.WaitGroup // tracks the writer goroutine
	writeErr error          // first error from the writer goroutine
	mu       sync.Mutex     // protects start()

	// stderrBuf captures samtools' stderr output across the lifetime of the
	// process. A dedicated goroutine tees samtools stderr into both this
	// buffer (for later diagnostic printing in Close) and cgltk's own
	// stderr (for real-time visibility). stderrWg tracks that goroutine so
	// Close can wait for it to finish before reading the buffer.
	stderrBuf bytes.Buffer
	stderrWg  sync.WaitGroup
}

// NewSamWriter creates a SamWriter for the given output file. For BAM format,
// a native writer is used (no samtools dependency). For CRAM or SAM format,
// samtools is required.
func NewSamWriter(filename string, opts *samWriterOptions) (SamWriter, error) {
	if opts.format == FormatBAM {
		if opts.sortedCoord || opts.sortedName {
			return newSortedBamWriter(filename, opts.header, opts.sortedCoord, opts.sortTmpPrefix)
		}
		return newBamWriter(filename, opts.header)
	}

	if err := checkSamtools(); err != nil {
		return nil, err
	}
	return &SamtoolsSamWriter{
		filename: filename,
		opts:     opts,
	}, nil
}

// Format sets the output format (FormatSAM, FormatBAM, or FormatCRAM).
func SamWriterOptions(h *SamHeader) *samWriterOptions {
	w := &samWriterOptions{}
	w.header = h
	return w
}

// Format sets the output format (FormatSAM, FormatBAM, or FormatCRAM).
func (w *samWriterOptions) BAM() *samWriterOptions {
	w.format = FormatBAM
	return w
}

// Format sets the output format (FormatSAM, FormatBAM, or FormatCRAM).
func (w *samWriterOptions) CRAM(ref string) *samWriterOptions {
	w.format = FormatCRAM
	w.reference = ref
	return w
}

// Threads sets the number of samtools compression threads (--threads).
func (w *samWriterOptions) Threads(n int) *samWriterOptions {
	w.threads = n
	return w
}

// Threads sets the number of samtools compression threads (--threads).
func (w *samWriterOptions) SortCoord() *samWriterOptions {
	w.sortedCoord = true
	return w
}

// Threads sets the number of samtools compression threads (--threads).
func (w *samWriterOptions) SortName() *samWriterOptions {
	w.sortedName = true
	return w
}

// Threads sets the number of samtools compression threads (--threads).
func (w *samWriterOptions) SortTempPrefix(prefix string) *samWriterOptions {
	w.sortTmpPrefix = prefix
	return w
}

func (w *SamtoolsSamWriter) start() error {
	w.mu.Lock()
	defer w.mu.Unlock()

	if w.started {
		return nil
	}
	var err error

	if w.opts.format == FormatSAM {
		w.outs, err = os.Create(w.filename)
		if err != nil {
			return err
		}
	} else {
		args := []string{}

		if w.opts.sortedCoord {
			args = append(args, "sort")
			if w.opts.sortTmpPrefix != "" {
				args = append(args, "-T", w.opts.sortTmpPrefix)
			}
		} else if w.opts.sortedName {
			args = append(args, "sort", "-N")
			if w.opts.sortTmpPrefix != "" {
				args = append(args, "-T", w.opts.sortTmpPrefix)
			}
		} else {
			args = append(args, "view")
		}
		switch w.opts.format {
		case FormatBAM:
			args = append(args, "-O", "bam")
		case FormatCRAM:
			args = append(args, "-O", "cram")
			if w.opts.reference != "" {
				args = append(args, "--reference", w.opts.reference)
			}
		}

		if w.opts.threads > 0 {
			args = append(args, "--threads", fmt.Sprintf("%d", w.opts.threads))
		}
		args = append(args, "--no-PG", "-o", w.filename, "-")
		w.cmd = exec.Command("samtools", args...)

		// Instrumentation point 1: log the full samtools argv once, right
		// before Start(). This pins down exactly what we invoked so the
		// command can be reproduced by hand from a job log.
		fmt.Fprintf(os.Stderr, "samtools cmd: samtools %s\n", strings.Join(args, " "))

		w.outs, err = w.cmd.StdinPipe()
		if err != nil {
			return fmt.Errorf("samtools stdin pipe: %w", err)
		}

		w.stderr, err = w.cmd.StderrPipe()
		if err != nil {
			return fmt.Errorf("samtools stderr pipe: %w", err)
		}

		if err := w.cmd.Start(); err != nil {
			return fmt.Errorf("samtools start: %w", err)
		}

		// Drain samtools stderr concurrently. Every line is tee'd to both
		// cgltk's stderr (with a "samtools: " prefix, so the user sees
		// samtools messages in real time) and an internal bytes.Buffer
		// (so Close can include the full text in any error return).
		//
		// This also eliminates the possibility of samtools blocking on a
		// full stderr pipe buffer — unlikely in practice with samtools
		// sort but a correctness concern worth closing off.
		w.stderrWg.Add(1)
		go func() {
			defer w.stderrWg.Done()
			scanner := bufio.NewScanner(w.stderr)
			scanner.Buffer(make([]byte, 64*1024), 1024*1024)
			for scanner.Scan() {
				line := scanner.Bytes()
				w.stderrBuf.Write(line)
				w.stderrBuf.WriteByte('\n')
				fmt.Fprintf(os.Stderr, "samtools: %s\n", line)
			}
			// Scanner errors (other than EOF) are intentionally ignored
			// here — the pipe may return ErrClosed after cmd.Wait, which
			// isn't something we want to escalate.
		}()
	}

	if w.opts.header != nil {
		headerText := w.opts.header.Text()
		if headerText != "" {
			if _, err := io.WriteString(w.outs, headerText); err != nil {
				return fmt.Errorf("write header: %w", err)
			}
		}
	}

	// Start buffered writer goroutine.
	w.writeCh = make(chan string, 1024)
	w.writeWg.Add(1)
	go func() {
		defer w.writeWg.Done()
		for line := range w.writeCh {
			if _, err := io.WriteString(w.outs, line); err != nil {
				w.writeErr = fmt.Errorf("write record: %w", err)
				// Drain remaining records to unblock senders.
				for range w.writeCh {
				}
				return
			}
		}
	}()

	w.started = true
	return nil
}

// Write serializes a SamRecord and sends it to the buffered write channel.
// It is safe for concurrent use from multiple goroutines.
func (w *SamtoolsSamWriter) Write(rec *SamRecord) error {
	if err := w.start(); err != nil {
		return err
	}
	if w.writeErr != nil {
		return w.writeErr
	}
	w.writeCh <- rec.String() + "\n"
	return nil
}

// Close drains the write buffer, waits for the writer goroutine to finish,
// and waits for the samtools process to exit.
//
// Instrumentation: Close logs the samtools exit status (success or error)
// to cgltk's stderr. The samtools stderr buffer is populated by the
// concurrent drain goroutine started in start(); we only need to wait for
// that goroutine to finish before reading the buffer in the error paths.
func (w *SamtoolsSamWriter) Close() error {
	if !w.started {
		return nil
	}
	close(w.writeCh)
	w.writeWg.Wait()
	if w.outs != nil {
		w.outs.Close()
	}

	// Wait for samtools to exit and its stderr drain goroutine to finish.
	// The drain goroutine returns as soon as the samtools process closes
	// its end of the stderr pipe (i.e., as soon as samtools exits).
	var waitErr error
	if w.cmd != nil {
		waitErr = w.cmd.Wait()
	}
	w.stderrWg.Wait()
	stderrText := w.stderrBuf.String()

	// Instrumentation point 3: log the samtools exit status unconditionally
	// on the success and error paths. If this line is missing from a job's
	// stderr, we know cmd.Wait() never returned — which means cgltk is
	// hanging in Close rather than exiting cleanly.
	if w.cmd != nil {
		if waitErr != nil {
			fmt.Fprintf(os.Stderr, "samtools exit: %v\n", waitErr)
		} else {
			fmt.Fprintf(os.Stderr, "samtools exit: ok\n")
		}
	}

	// If the writer goroutine saw an EPIPE (samtools died before consuming
	// our writes), surface that as the primary error but also include the
	// samtools exit status and any stderr text so we don't lose the real
	// cause.
	if w.writeErr != nil {
		if waitErr != nil {
			return fmt.Errorf("%w (samtools: %v: %s)", w.writeErr, waitErr, stderrText)
		}
		return fmt.Errorf("%w (samtools stderr: %s)", w.writeErr, stderrText)
	}

	// No writer error. If samtools exited non-zero, surface that along
	// with the captured stderr.
	if waitErr != nil {
		return fmt.Errorf("samtools: %w: %s", waitErr, stderrText)
	}
	return nil
}

// StdoutSamWriter writes SAM text directly to stdout without samtools.
type StdoutSamWriter struct {
	header        *SamHeader
	headerWritten bool
}

// NewStdoutSamWriter creates a writer that writes SAM text to stdout.
func NewStdoutSamWriter(header *SamHeader) *StdoutSamWriter {
	return &StdoutSamWriter{
		header: header,
	}
}

// Write writes a SamRecord to stdout as SAM text.
func (w *StdoutSamWriter) Write(rec *SamRecord) error {
	if !w.headerWritten {
		if w.header != nil {
			headerText := w.header.Text()
			if headerText != "" {
				if _, err := io.WriteString(os.Stdout, headerText); err != nil {
					return fmt.Errorf("write header: %w", err)
				}
			}
		}
		w.headerWritten = true
	}
	_, err := fmt.Fprintln(os.Stdout, rec.String())
	if err != nil {
		return fmt.Errorf("write record: %w", err)
	}
	return nil
}

// Close is a no-op for stdout.
func (w *StdoutSamWriter) Close() error {
	return nil
}
