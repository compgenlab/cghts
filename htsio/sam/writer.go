package sam

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sync"

	"github.com/compgen-io/cgkit/htsio"
)

// Writer writes SAM text files. Implements htsio.SamWriter.
// It is safe for concurrent use by multiple goroutines.
type Writer struct {
	mu            sync.Mutex
	w             *bufio.Writer
	f             *os.File // non-nil if we opened the file
	header        *htsio.SamHeader
	headerWritten bool
}

// NewWriter creates a SAM text writer for the given output file.
// If filename is "-", writes to stdout.
func NewWriter(filename string, header *htsio.SamHeader) (*Writer, error) {
	if filename == "-" {
		return NewWriterFromWriter(os.Stdout, header), nil
	}
	f, err := os.Create(filename)
	if err != nil {
		return nil, err
	}
	return &Writer{
		w:      bufio.NewWriter(f),
		f:      f,
		header: header,
	}, nil
}

// NewWriterFromWriter creates a SAM text writer that writes to w.
func NewWriterFromWriter(w io.Writer, header *htsio.SamHeader) *Writer {
	return &Writer{
		w:      bufio.NewWriter(w),
		header: header,
	}
}

// Write writes a SamRecord as SAM text.
// It is safe for concurrent use by multiple goroutines.
func (w *Writer) Write(rec *htsio.SamRecord) error {
	w.mu.Lock()
	defer w.mu.Unlock()

	if !w.headerWritten {
		if w.header != nil {
			headerText := w.header.Text()
			if headerText != "" {
				if _, err := w.w.WriteString(headerText); err != nil {
					return fmt.Errorf("write header: %w", err)
				}
			}
		}
		w.headerWritten = true
	}
	if _, err := fmt.Fprintln(w.w, rec.String()); err != nil {
		return fmt.Errorf("write record: %w", err)
	}
	return nil
}

// Close flushes the buffer and closes the file (if not stdout).
func (w *Writer) Close() error {
	w.mu.Lock()
	defer w.mu.Unlock()

	if err := w.w.Flush(); err != nil {
		return err
	}
	if w.f != nil {
		return w.f.Close()
	}
	return nil
}
