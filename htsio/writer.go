package htsio

import (
	"fmt"
	"io"
	"os"
	"os/exec"
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

// SamtoolsSamWriter writes SAM/BAM/CRAM files by piping SAM text to samtools view.
type SamtoolsSamWriter struct {
	filename string
	format   SamOutputFormat
	header   *SamHeader
	refFile  string // CRAM reference FASTA
	cmd      *exec.Cmd
	stdin    io.WriteCloser
	stderr   io.ReadCloser
	started  bool
}

// NewSamWriter creates a SamtoolsSamWriter for the given output file.
// Returns an error if samtools is not found in PATH.
// Default format is BAM. Use the builder methods to set options before calling Write().
func NewSamWriter(filename string, header *SamHeader) (*SamtoolsSamWriter, error) {
	if err := checkSamtools(); err != nil {
		return nil, err
	}
	return &SamtoolsSamWriter{
		filename: filename,
		header:   header,
		format:   FormatBAM,
	}, nil
}

// Format sets the output format (FormatSAM, FormatBAM, or FormatCRAM).
func (w *SamtoolsSamWriter) Format(f SamOutputFormat) *SamtoolsSamWriter {
	w.format = f
	return w
}

// Reference sets the reference FASTA file, required for CRAM output.
func (w *SamtoolsSamWriter) Reference(ref string) *SamtoolsSamWriter {
	w.refFile = ref
	return w
}

func (w *SamtoolsSamWriter) start() error {
	if w.started {
		return nil
	}

	args := []string{"view", "-S"}

	switch w.format {
	case FormatSAM:
		args = append(args, "-h")
	case FormatBAM:
		args = append(args, "-b")
	case FormatCRAM:
		args = append(args, "-C")
		if w.refFile != "" {
			args = append(args, "-T", w.refFile)
		}
	}

	args = append(args, "-o", w.filename, "-")

	w.cmd = exec.Command("samtools", args...)

	var err error
	w.stdin, err = w.cmd.StdinPipe()
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

	w.started = true

	if w.header != nil {
		headerText := w.header.Text()
		if headerText != "" {
			if _, err := io.WriteString(w.stdin, headerText); err != nil {
				return fmt.Errorf("write header: %w", err)
			}
		}
	}

	return nil
}

// Write writes a SamRecord to the output.
func (w *SamtoolsSamWriter) Write(rec *SamRecord) error {
	if err := w.start(); err != nil {
		return err
	}
	_, err := fmt.Fprintln(w.stdin, rec.String())
	if err != nil {
		return fmt.Errorf("write record: %w", err)
	}
	return nil
}

// Close flushes remaining data and waits for the samtools process to finish.
func (w *SamtoolsSamWriter) Close() error {
	if !w.started {
		return nil
	}
	if w.stdin != nil {
		w.stdin.Close()
	}
	if w.stderr != nil {
		stderrBytes, _ := io.ReadAll(w.stderr)
		w.stderr.Close()
		if err := w.cmd.Wait(); err != nil {
			return fmt.Errorf("samtools: %w: %s", err, string(stderrBytes))
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
