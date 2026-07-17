package vcf

import (
	"bufio"
	"io"
	"os"
	"strings"

	"github.com/compgenlab/cghts/htsio/bgzf"
)

// VcfWriter writes a VCF header and records to an io.Writer or a file.
type VcfWriter struct {
	writer *bufio.Writer
	gz     io.WriteCloser // BGZF compression layer, when the file is block-gzipped
	file   *os.File
}

// NewVcfWriter creates a VcfWriter that writes to w.
func NewVcfWriter(w io.Writer) *VcfWriter {
	return &VcfWriter{writer: bufio.NewWriter(w)}
}

// OpenVcfWriter creates a VcfWriter for the given filename. A filename ending in
// ".gz" or ".bgz" is BGZF (block-gzip) compressed, so the output is a valid bgzip
// file that tabix can index — not plain gzip.
func OpenVcfWriter(filename string) (*VcfWriter, error) {
	f, err := os.Create(filename)
	if err != nil {
		return nil, err
	}
	w := &VcfWriter{file: f}
	if strings.HasSuffix(filename, ".gz") || strings.HasSuffix(filename, ".bgz") {
		w.gz = bgzf.NewWriter(f)
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

// WriteRecord writes a record. An unmodified record is emitted verbatim; a
// modified one (see [VcfRecord.Dirty]) is reconstructed from its parsed model.
func (w *VcfWriter) WriteRecord(rec *VcfRecord) error {
	if rec.dirty {
		return w.WriteLine(rec.serialize())
	}
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
