package htsio

import (
	"fmt"
	"io"
	"os"
	"os/exec"
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
}

// NewSamWriter creates a SamtoolsSamWriter for the given output file.
// Returns an error if samtools is not found in PATH.
// Default format is BAM. Use the builder methods to set options before calling Write().
func NewSamWriter(filename string, opts *samWriterOptions) (*SamtoolsSamWriter, error) {
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
func (w *SamtoolsSamWriter) Close() error {
	if !w.started {
		return nil
	}
	close(w.writeCh)
	w.writeWg.Wait()
	if w.outs != nil {
		w.outs.Close()
	}
	if w.writeErr != nil {
		// Still wait for samtools to exit.
		if w.stderr != nil {
			io.ReadAll(w.stderr)
			w.stderr.Close()
		}
		if w.cmd != nil {
			w.cmd.Wait()
		}
		return w.writeErr
	}
	if w.stderr != nil {
		stderrBytes, _ := io.ReadAll(w.stderr)
		w.stderr.Close()
		if w.cmd != nil {
			if err := w.cmd.Wait(); err != nil {
				return fmt.Errorf("samtools: %w: %s", err, string(stderrBytes))
			}
		}
	} else if w.cmd != nil {
		return w.cmd.Wait()
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
