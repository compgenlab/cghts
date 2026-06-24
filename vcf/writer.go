package vcf

import (
	"bufio"
	"compress/gzip"
	"io"
	"os"
	"strings"
)

// VcfWriter writes a VCF header and records to an io.Writer or a file.
type VcfWriter struct {
	writer *bufio.Writer
	gz     *gzip.Writer
	file   *os.File
}

// NewVcfWriter creates a VcfWriter that writes to w.
func NewVcfWriter(w io.Writer) *VcfWriter {
	return &VcfWriter{writer: bufio.NewWriter(w)}
}

// OpenVcfWriter creates a VcfWriter for the given filename. A filename ending in
// ".gz" is gzip-compressed.
func OpenVcfWriter(filename string) (*VcfWriter, error) {
	f, err := os.Create(filename)
	if err != nil {
		return nil, err
	}
	w := &VcfWriter{file: f}
	if strings.HasSuffix(filename, ".gz") {
		w.gz = gzip.NewWriter(f)
		w.writer = bufio.NewWriter(w.gz)
	} else {
		w.writer = bufio.NewWriter(f)
	}
	return w, nil
}

// WriteHeader writes the header's metadata lines and the #CHROM column line.
func (w *VcfWriter) WriteHeader(h *VcfHeader) error {
	_, err := h.WriteTo(w.writer)
	return err
}

// WriteRecord writes a record's raw line verbatim.
func (w *VcfWriter) WriteRecord(rec *VcfRecord) error {
	return w.WriteLine(rec.Line())
}

// WriteLine writes a single raw line, appending a newline.
func (w *VcfWriter) WriteLine(line string) error {
	if _, err := w.writer.WriteString(line); err != nil {
		return err
	}
	return w.writer.WriteByte('\n')
}

// Close flushes and closes the writer.
func (w *VcfWriter) Close() error {
	if w.writer != nil {
		if err := w.writer.Flush(); err != nil {
			return err
		}
	}
	if w.gz != nil {
		if err := w.gz.Close(); err != nil {
			return err
		}
	}
	if w.file != nil {
		return w.file.Close()
	}
	return nil
}
