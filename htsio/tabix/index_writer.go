package tabix

import (
	"io"
	"os"

	"github.com/compgenlab/cghts/htsio/bgzf"
)

// IndexWriter builds a tabix (.tbi) index for a file that is already
// BGZF-compressed and sorted (the file is not modified). The column positions,
// meta/comment character, header skip count, and coordinate base come from a
// [WriterOpts] (use its BED/VCF/GFF presets or set fields directly), matching
// the configuration the `tabix` command line accepts.
type IndexWriter struct {
	opts *WriterOpts
}

// NewIndexWriter returns an IndexWriter configured by opts. A nil opts uses the
// defaults from [NewWriterOpts].
func NewIndexWriter(opts *WriterOpts) *IndexWriter {
	if opts == nil {
		opts = NewWriterOpts()
	}
	return &IndexWriter{opts: opts}
}

// WriteIndex reads the BGZF-compressed file and writes a companion ".tbi" index
// (filename + ".tbi"). It walks the file block by block, recording the virtual
// offset at the start of each data line (skipping the configured header lines
// and meta/comment lines).
func (iw *IndexWriter) WriteIndex(filename string) error {
	f, err := os.Open(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	r := bgzf.NewReader(f)
	ib := &tbiIndexBuilder{opts: iw.opts, refs: make(map[string]*tbiRefBuilder)}
	refIdx := make(map[string]int)
	var refOrder []string

	lineNo := 0
	for {
		begin := r.VirtualTell()
		line, rerr := readBGZFLine(r)
		end := r.VirtualTell()

		// A line is available when we read up to a newline (rerr == nil) or hit
		// EOF with a final, newline-less line (line != "").
		if rerr == nil || line != "" {
			switch {
			case lineNo < int(iw.opts.skip):
				// header line, not indexed
			case line == "":
				// blank line
			case iw.opts.meta != 0 && line[0] == byte(iw.opts.meta):
				// comment line
			default:
				l, perr := parseTabixLine(line, iw.opts)
				if perr != nil {
					return perr
				}
				if _, seen := refIdx[l.ref]; !seen {
					refIdx[l.ref] = len(refOrder)
					refOrder = append(refOrder, l.ref)
				}
				ib.addRecord(l, begin, end)
			}
			lineNo++
		}

		if rerr != nil {
			if rerr == io.EOF {
				break
			}
			return rerr
		}
	}

	ib.refOrder = refOrder
	ib.refIdx = refIdx
	return ib.writeTBI(filename + ".tbi")
}

// readBGZFLine reads a single line (without the trailing newline) from r. It
// returns io.EOF when the stream ends; the returned string holds the final line
// when it is not newline-terminated.
func readBGZFLine(r *bgzf.Reader) (string, error) {
	var sb []byte
	for {
		b, err := r.ReadByte()
		if err != nil {
			return string(sb), err
		}
		if b == '\n' {
			return string(sb), nil
		}
		sb = append(sb, b)
	}
}
