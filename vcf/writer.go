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
	// Every stage runs even after one fails, and the first error is what is
	// reported. Returning early left the layers below open: the file descriptor
	// leaked, and for a .gz output the BGZF writer never ran its own Close, so
	// the EOF block went unwritten and the result read as truncated to every
	// bgzip and tabix consumer -- on precisely the error path where releasing
	// the handle matters most.
	var first error
	fail := func(err error) {
		if err != nil && first == nil {
			first = err
		}
	}
	if w.writer != nil {
		fail(w.writer.Flush())
	}
	if w.gz != nil {
		fail(w.gz.Close())
	}
	if w.file != nil {
		fail(w.file.Close())
	}
	return first
}
