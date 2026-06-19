package bed

import (
	"bufio"
	"compress/gzip"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/compgenlab/hts/htsio/tabix"
)

// ColumnMode selects how many columns a [BedWriter] emits per record.
type ColumnMode int

const (
	// ColumnsAuto emits ref/start/end and, when the record HasName, the
	// name/score/strand fields followed by any extra columns. This is the
	// faithful pass-through layout.
	ColumnsAuto ColumnMode = iota
	// Columns3 emits a strict BED3 line (ref/start/end only); extras are dropped.
	Columns3
	// Columns6 emits a strict BED6 line (ref/start/end/name/score/strand);
	// extras are dropped.
	Columns6
)

// BedWriterOpts configures a [BedWriter].
type BedWriterOpts struct {
	forceScoreInt bool
	mode          ColumnMode
	tabix         bool
	bgzip         bool
	rawStrand     bool
}

// NewBedWriterOpts returns a BedWriterOpts with default settings (auto columns,
// float scores, plain output).
func NewBedWriterOpts() *BedWriterOpts {
	return &BedWriterOpts{}
}

// ForceScoreInt coerces the score column to an integer (truncating toward
// zero) on output. It returns o to allow chaining.
func (o *BedWriterOpts) ForceScoreInt(b bool) *BedWriterOpts {
	o.forceScoreInt = b
	return o
}

// Columns selects the output column layout. It returns o to allow chaining.
func (o *BedWriterOpts) Columns(m ColumnMode) *BedWriterOpts {
	o.mode = m
	return o
}

// RawStrand emits the strand field verbatim, including "." for [StrandNone],
// instead of the default behavior of forcing an unspecified strand to "+". It
// returns o to allow chaining.
func (o *BedWriterOpts) RawStrand(b bool) *BedWriterOpts {
	o.rawStrand = b
	return o
}

// Tabix requests tabix-indexed output: a sorted, BGZF-compressed file with a
// companion .tbi index. It is honored only by [OpenBedWriter] (which has a
// filename to index); [NewBedWriter] cannot produce tabix output. It returns o
// to allow chaining.
func (o *BedWriterOpts) Tabix(b bool) *BedWriterOpts {
	o.tabix = b
	return o
}

// Bgzip requests sorted, BGZF-compressed output without a tabix index. Like
// [BedWriterOpts.Tabix] it is honored only by [OpenBedWriter]. [BedWriterOpts.Tabix]
// takes precedence (it already implies BGZF). It returns o to allow chaining.
func (o *BedWriterOpts) Bgzip(b bool) *BedWriterOpts {
	o.bgzip = b
	return o
}

// BedWriter writes BED records one at a time to an io.Writer or a file.
type BedWriter struct {
	writer *bufio.Writer
	gz     *gzip.Writer
	file   *os.File
	tw     *tabix.Writer
	opts   *BedWriterOpts
}

func resolveOpts(opts []*BedWriterOpts) *BedWriterOpts {
	if len(opts) > 0 && opts[0] != nil {
		return opts[0]
	}
	return NewBedWriterOpts()
}

// NewBedWriter creates a BedWriter that writes to w. The Tabix option is not
// supported here (it requires a filename to index); use [OpenBedWriter] for
// tabix output.
func NewBedWriter(w io.Writer, opts ...*BedWriterOpts) *BedWriter {
	return &BedWriter{writer: bufio.NewWriter(w), opts: resolveOpts(opts)}
}

// OpenBedWriter creates a BedWriter for the given filename. If the Tabix option
// is set, output is a sorted, BGZF-compressed file with a companion .tbi index.
// Otherwise, if the filename ends in ".gz" the output is whole-file
// gzip-compressed; otherwise it is plain text.
func OpenBedWriter(filename string, opts ...*BedWriterOpts) (*BedWriter, error) {
	o := resolveOpts(opts)

	if o.tabix {
		tw := tabix.NewWriter(filename, tabix.NewWriterOpts().BED().AutoIndex())
		return &BedWriter{tw: tw, opts: o}, nil
	}
	if o.bgzip {
		// Sorted BGZF without a .tbi index, reusing the tabix writer's sort.
		tw := tabix.NewWriter(filename, tabix.NewWriterOpts().BED())
		return &BedWriter{tw: tw, opts: o}, nil
	}

	f, err := os.Create(filename)
	if err != nil {
		return nil, err
	}
	w := &BedWriter{file: f, opts: o}
	if strings.HasSuffix(filename, ".gz") {
		w.gz = gzip.NewWriter(f)
		w.writer = bufio.NewWriter(w.gz)
	} else {
		w.writer = bufio.NewWriter(f)
	}
	return w, nil
}

// WriteRecord writes a single BED record using the writer's configured layout.
func (w *BedWriter) WriteRecord(rec *BedRecord) error {
	line := w.formatLine(rec)
	if w.tw != nil {
		// The tabix writer appends its own newline.
		return w.tw.Write(line)
	}
	if _, err := w.writer.WriteString(line); err != nil {
		return err
	}
	return w.writer.WriteByte('\n')
}

func (w *BedWriter) formatLine(rec *BedRecord) string {
	fields := []string{rec.Ref, strconv.Itoa(rec.Start), strconv.Itoa(rec.End)}

	emit6 := func() {
		fields = append(fields, rec.Name)
		if w.opts.forceScoreInt {
			fields = append(fields, formatScoreInt(rec.Score))
		} else {
			fields = append(fields, formatScoreFloat(rec.Score))
		}
		strand := rec.Strand
		if strand == StrandNone && !w.opts.rawStrand {
			strand = StrandPlus
		}
		fields = append(fields, string(strand))
	}

	switch w.opts.mode {
	case Columns3:
		// ref/start/end only.
	case Columns6:
		emit6()
	default: // ColumnsAuto
		if rec.HasName {
			emit6()
			fields = append(fields, rec.Extras...)
		}
	}

	return strings.Join(fields, "\t")
}

// Close flushes and closes the writer, finalizing tabix output (sort + .tbi)
// when in tabix mode.
func (w *BedWriter) Close() error {
	if w.tw != nil {
		return w.tw.Close()
	}
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
