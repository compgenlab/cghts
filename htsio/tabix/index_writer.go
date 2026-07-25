package tabix

import (
	"fmt"
	"io"
	"os"

	"github.com/compgenlab/cghts/htsio/bgzf"
)

// UnsortedError reports that a file cannot be indexed because its records are
// not in coordinate order. A tabix index over an unsorted file is not an error
// to build but is wrong to use: its linear index assumes offsets increase with
// position, so queries silently miss records. [IndexWriter.WriteIndex] returns
// this instead of writing such an index.
type UnsortedError struct {
	Line    int    // 1-based line number of the offending record
	Ref     string // its reference/chromosome name
	Pos     int    // its start position, in the input's coordinate base
	PrevRef string // reference of the preceding record
	PrevPos int    // start position of the preceding record
	Revisit bool   // Ref reappeared after records on other references
}

func (e *UnsortedError) Error() string {
	if e.Revisit {
		return fmt.Sprintf("not coordinate sorted: line %d has %s:%d, but %s already gave way to %s",
			e.Line, e.Ref, e.Pos, e.Ref, e.PrevRef)
	}
	return fmt.Sprintf("not coordinate sorted: line %d has %s:%d, which follows %s:%d",
		e.Line, e.Ref, e.Pos, e.PrevRef, e.PrevPos)
}

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
//
// The records must be in coordinate order: positions ascending within a
// reference, and each reference in one contiguous block. Otherwise no index is
// written and the error is an [*UnsortedError].
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

	// Previous record, for the order check. lastRef is "" before the first one.
	lastRef, lastStart := "", 0

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
				_, seen := refIdx[l.ref]
				switch {
				case l.ref != lastRef:
					// A reference that was already left behind cannot resume:
					// the linear index assumes offsets grow with position.
					if seen {
						return iw.unsorted(lineNo+1, l, lastRef, lastStart, true)
					}
					refIdx[l.ref] = len(refOrder)
					refOrder = append(refOrder, l.ref)
				case l.start < lastStart:
					return iw.unsorted(lineNo+1, l, lastRef, lastStart, false)
				}
				lastRef, lastStart = l.ref, l.start
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

// unsorted builds the UnsortedError for a record that breaks coordinate order,
// reporting positions in the input's own coordinate base rather than the
// 0-based one the index uses internally.
func (iw *IndexWriter) unsorted(lineNo int, l tabixLine, prevRef string, prevStart int, revisit bool) error {
	base := 0
	if !iw.opts.zeroBased {
		base = 1
	}
	return &UnsortedError{
		Line:    lineNo,
		Ref:     l.ref,
		Pos:     l.start + base,
		PrevRef: prevRef,
		PrevPos: prevStart + base,
		Revisit: revisit,
	}
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
