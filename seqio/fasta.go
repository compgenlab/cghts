package seqio

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"iter"
	"os"
	"strings"
)

type FastaSeqReader struct {
	buf    string
	reader *bufio.Reader
}

type FastaReader struct {
	filename string
	file     *os.File
	parent   io.Reader
	reader   *bufio.Reader
	closed   bool
	buffer   string
	lastByte byte
}

func (r *FastaReader) Close() {
	if r.file != nil {
		r.file.Close()
	}
	r.closed = true
}

func NewFastaFile(filename string) (*FastaReader, error) {
	f, err := os.Open(filename) // read-only
	if err != nil {
		return nil, err
	}

	// For indexed random access, use NewIndexedFastaReader() instead.

	r := bufio.NewReader(f)
	magic, err := r.Peek(2)
	if err != nil {
		f.Close()
		return nil, err
	}

	if magic[0] == 0x1f && magic[1] == 0x8b {
		// gzip magic number, wrap in a gzip.Reader
		gz, err := gzip.NewReader(r)
		if err != nil {
			f.Close()
			return nil, err
		}
		r = bufio.NewReader(gz)
	}

	return &FastaReader{
		filename: filename,
		file:     f,
		parent:   f,
		closed:   false,
		reader:   r,
		lastByte: '\n',
	}, nil
}

func NewFastaReader(rd io.Reader) (*FastaReader, error) {
	if rd == nil {
		return nil, io.ErrUnexpectedEOF
	}
	return &FastaReader{
		filename: "",
		file:     nil,
		parent:   rd,
		closed:   false,
		reader:   bufio.NewReader(rd),
		lastByte: '\n',
	}, nil
}

func (r *FastaReader) NextSeq() (SeqRecord, error) {
	if r.closed {
		return nil, ClosedSeqReaderError
	}

	for {
		if b, err := r.reader.ReadByte(); err != nil {
			return nil, err
		} else if b == '>' && r.lastByte == '\n' {
			// we only are at the start of a new record if we see a '>' at the start of a line (after a newline)
			// technically a '>' at any other place wouldn't be a valid FASTA file, but we'll ignore that for now.
			line, err := r.reader.ReadString('\n')

			if err != nil && line == "" {
				// if we didn't read anything and got an error, that's a problem.
				// if we get an error, but read something, that means we didn't end with
				// a newline, but that's not a problem, we can still process the record we read.
				return nil, err
			}

			var name, comment string

			line = strings.Trim(line, "\r\n")
			spl := strings.SplitN(line, " ", 2)

			name = spl[0]
			if len(spl) > 1 {
				comment = spl[1]
			}
			r.lastByte = '\n' // we know the last byte is a newline since we just read it
			return &FastaSeqRecord{
				name:     name,
				comment:  comment,
				reader:   r.reader,
				lastByte: '\n',
			}, nil
		} else {
			r.lastByte = b
		}
	}
}

func (r *FastaReader) Names() (iter.Seq[string], error) {
	if r.closed {
		return nil, ClosedSeqReaderError
	}
	return func(yield func(string) bool) {
		for seq, err := r.NextSeq(); err == nil; {
			if !yield(seq.Name()) {
				return
			}
		}
		// this isn't a seek-able reader, so we can't reset it, so close it instead
		r.Close()
	}, nil
}

func (r *FastaReader) FetchRecord(name string) (SeqRecord, error) {
	if r.closed {
		return nil, ClosedSeqReaderError
	}
	for seq, err := r.NextSeq(); err == nil; {
		if seq.Name() == name {
			return seq, nil
		}
	}
	// close the reader since we can't reset it
	r.Close()

	return nil, io.EOF
}

// FASTA records will extract the seq from a reader,
// so they don't have to be fully loaded into memory at once
// (FASTA records can be very large).
type FastaSeqRecord struct {
	name    string
	comment string
	reader  *bufio.Reader
	// lastByte tracks the byte most recently consumed by Chunks so that we can
	// detect the \n> delimiter across Peek boundaries. Initialized to '\n'
	// since the preceding byte is always the header line's terminating newline.
	lastByte byte
}

func (r *FastaSeqRecord) FullSeq() SeqQual {
	// Read the entire sequence into memory
	var buf strings.Builder
	var last byte = '\n'
	for {
		b, err := r.reader.ReadByte()
		if err != nil {
			break
		}
		if last == '\n' && b == '>' {
			r.reader.UnreadByte()
			break
		}
		if b != '\n' && b != '\r' {
			buf.WriteByte(b)
		}
		last = b
	}
	return SeqQual{
		name: r.name,
		seq:  buf.String(),
		pos:  0,
	}
}

// Chunks implements [SeqRecord]. The sequence is streamed from the underlying
// reader in chunks of at most `length` bytes (newlines stripped). On each
// iteration, the next `length` bytes are peeked and scanned for the \n>
// delimiter that marks the start of the next record. If found, only the bytes
// up to that point are consumed, so the reader is left positioned at the next
// record's header.
func (r *FastaSeqRecord) Chunks(length int) iter.Seq[SeqQual] {
	return func(yield func(SeqQual) bool) {
		if length <= 0 {
			return
		}
		curPos := 0
		for {
			peek, _ := r.reader.Peek(length)
			if len(peek) == 0 {
				return
			}

			// Find the first \n> delimiter within the peek window. `prev` is
			// seeded from r.lastByte so that a '>' at peek[0] is detected when
			// the preceding byte (already consumed) was a newline.
			readLen := len(peek)
			prev := r.lastByte
			for i, b := range peek {
				if prev == '\n' && b == '>' {
					readLen = i
					break
				}
				prev = b
			}

			if readLen == 0 {
				// Reader is already at the start of the next record.
				return
			}

			buf := make([]byte, readLen)
			n, err := io.ReadFull(r.reader, buf)
			if n > 0 {
				r.lastByte = buf[n-1]
				chunk := string(buf[:n])
				chunk = strings.ReplaceAll(chunk, "\n", "")
				chunk = strings.ReplaceAll(chunk, "\r", "")
				if len(chunk) > 0 {
					if !yield(SeqQual{
						seq:  chunk,
						name: r.name,
						pos:  curPos,
					}) {
						return
					}
					curPos += len(chunk)
				}
			}
			if err != nil {
				return
			}
		}
	}
}

func (r *FastaSeqRecord) Name() string {
	return r.name
}

func (r *FastaSeqRecord) Comment() string {
	return r.comment
}

// FastaWriterOpts configures a FastaWriter.
type FastaWriterOpts struct {
	wrap int
}

// NewFastaWriterOpts returns a new FastaWriterOpts with default settings.
func NewFastaWriterOpts() *FastaWriterOpts {
	return &FastaWriterOpts{}
}

// Wrap sets the line wrap length for output sequences. 0 means no wrapping.
func (o *FastaWriterOpts) Wrap(n int) *FastaWriterOpts {
	o.wrap = n
	return o
}

// FastaWriter writes FASTA records to a file, optionally gzip-compressed.
type FastaWriter struct {
	writer *bufio.Writer
	gz     *gzip.Writer
	file   *os.File
	opts   *FastaWriterOpts
}

// NewFastaWriter creates a FastaWriter that writes to the given io.Writer.
func NewFastaWriter(w io.Writer, opts ...*FastaWriterOpts) *FastaWriter {
	var o *FastaWriterOpts
	if len(opts) > 0 && opts[0] != nil {
		o = opts[0]
	} else {
		o = NewFastaWriterOpts()
	}
	return &FastaWriter{writer: bufio.NewWriter(w), opts: o}
}

// OpenFastaWriter creates a FastaWriter for the given filename.
// If the filename ends in ".gz", the output will be gzip-compressed.
func OpenFastaWriter(filename string, opts ...*FastaWriterOpts) (*FastaWriter, error) {
	f, err := os.Create(filename)
	if err != nil {
		return nil, err
	}
	var o *FastaWriterOpts
	if len(opts) > 0 && opts[0] != nil {
		o = opts[0]
	} else {
		o = NewFastaWriterOpts()
	}
	w := &FastaWriter{file: f, opts: o}
	if strings.HasSuffix(filename, ".gz") {
		gz := gzip.NewWriter(f)
		w.gz = gz
		w.writer = bufio.NewWriter(gz)
	} else {
		w.writer = bufio.NewWriter(f)
	}
	return w, nil
}

// WriteRecord writes a single FASTA record.
func (w *FastaWriter) WriteRecord(name, comment, seq string) error {
	if comment != "" {
		if _, err := fmt.Fprintf(w.writer, ">%s %s\n", name, comment); err != nil {
			return err
		}
	} else {
		if _, err := fmt.Fprintf(w.writer, ">%s\n", name); err != nil {
			return err
		}
	}
	if w.opts.wrap > 0 {
		for i := 0; i < len(seq); i += w.opts.wrap {
			end := i + w.opts.wrap
			if end > len(seq) {
				end = len(seq)
			}
			if _, err := fmt.Fprintf(w.writer, "%s\n", seq[i:end]); err != nil {
				return err
			}
		}
	} else {
		if _, err := fmt.Fprintf(w.writer, "%s\n", seq); err != nil {
			return err
		}
	}
	return nil
}

// WriteSeq writes a SeqRecord, streaming chunks to avoid loading the full sequence into memory.
func (w *FastaWriter) WriteSeq(rec SeqRecord) error {
	comment := rec.Comment()
	if comment != "" {
		if _, err := fmt.Fprintf(w.writer, ">%s %s\n", rec.Name(), comment); err != nil {
			return err
		}
	} else {
		if _, err := fmt.Fprintf(w.writer, ">%s\n", rec.Name()); err != nil {
			return err
		}
	}

	if w.opts.wrap > 0 {
		buf := ""
		for chunk := range rec.Chunks(w.opts.wrap) {
			buf += chunk.Seq()
			for len(buf) >= w.opts.wrap {
				if _, err := fmt.Fprintf(w.writer, "%s\n", buf[:w.opts.wrap]); err != nil {
					return err
				}
				buf = buf[w.opts.wrap:]
			}
		}
		if len(buf) > 0 {
			if _, err := fmt.Fprintf(w.writer, "%s\n", buf); err != nil {
				return err
			}
		}
	} else {
		for chunk := range rec.Chunks(1024) {
			if _, err := w.writer.WriteString(chunk.Seq()); err != nil {
				return err
			}
		}
		if _, err := w.writer.WriteString("\n"); err != nil {
			return err
		}
	}
	return nil
}

// Close flushes and closes the writer.
func (w *FastaWriter) Close() error {
	if err := w.writer.Flush(); err != nil {
		return err
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
