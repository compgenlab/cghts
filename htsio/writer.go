package htsio

import (
	"fmt"
	"io"
	"os"
)

// SamWriter is the interface for writing SAM/BAM/CRAM records.
type SamWriter interface {
	// Write writes a SamRecord to the output.
	Write(rec *SamRecord) error
	// Close flushes remaining data and releases resources.
	Close() error
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
