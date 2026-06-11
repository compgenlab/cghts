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

type FastqReader struct {
	filename string
	file     *os.File
	parent   io.Reader
	reader   *bufio.Reader
	closed   bool
	buffer   string
	lastByte byte
}

func (r *FastqReader) Close() {
	if r.file != nil {
		r.file.Close()
	}
	r.closed = true
}

func NewFastqFile(filename string) (*FastqReader, error) {
	f, err := os.Open(filename) // read-only
	if err != nil {
		return nil, err
	}

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
	return &FastqReader{
		filename: filename,
		file:     f,
		parent:   f,
		closed:   false,
		reader:   r,
		lastByte: '\n',
	}, nil
}

func NewFastqReader(rd io.Reader) (*FastqReader, error) {
	if rd == nil {
		return nil, io.ErrUnexpectedEOF
	}
	return &FastqReader{
		filename: "",
		file:     nil,
		parent:   rd,
		closed:   false,
		reader:   bufio.NewReader(rd),
		lastByte: '\n',
	}, nil
}

func (r *FastqReader) NextFastqSeq() (*FastqSeqRecord, error) {
	if r.closed {
		return nil, ClosedSeqReaderError
	}
	for {
		if b, err := r.reader.ReadByte(); err != nil {
			return nil, err
		} else if b == '@' && r.lastByte == '\n' {
			// fmt.Printf("Byte: %c\n", b)

			// TODO: make this support multi-line FASTQ records.
			//       Right now it only supports single-line sequences
			//       and quality strings, which is the most common format,
			//       but we should support multi-line records.

			// start of a FASTQ record, read each line
			header, err := r.reader.ReadString('\n')
			// fmt.Printf("header: %s\n", header)
			if err != nil {
				return nil, err
			}
			seq, err := r.reader.ReadString('\n')
			// fmt.Printf("seq: %s\n", seq)
			if err != nil {
				return nil, err
			}
			plus, err := r.reader.ReadString('\n')
			// fmt.Printf("plus: %s\n", plus)
			if err != nil {
				return nil, err
			}
			if plus[0] != '+' {
				return nil, io.ErrUnexpectedEOF
			}
			qual, err := r.reader.ReadString('\n')
			// fmt.Printf("qual: %s\n", qual)
			if err != nil && qual == "" {
				// it could be possible to have a quality string
				// without a newline, so we'll allow that, but
				// if we got an error and we didn't get any quality string, then it's an error
				return nil, err
			}

			header = strings.Trim(header, "\r\n")
			seq = strings.Trim(seq, "\r\n")
			qual = strings.Trim(qual, "\r\n")

			spl := strings.SplitN(header, " ", 2)

			name := spl[0]
			comment := ""
			if len(spl) > 1 {
				comment = spl[1]
			}

			r.lastByte = '\n' // we know the last byte is a newline since we just read it
			return &FastqSeqRecord{
				name:    name,
				comment: comment,
				seq:     seq,
				qual:    qual,
			}, nil
		} else {
			r.lastByte = b
		}
	}
}

func (r *FastqReader) NextSeq() (SeqRecord, error) {
	return r.NextFastqSeq()
}

func (r *FastqReader) Names() (iter.Seq[string], error) {
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

func (r *FastqReader) FetchRecord(name string) (SeqRecord, error) {
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

// FASTQ records store the full record in memory
type FastqSeqRecord struct {
	name    string
	comment string
	seq     string
	qual    string
}

func (r *FastqSeqRecord) FullSeq() SeqQual {
	return SeqQual{
		name: r.name,
		seq:  r.seq,
		qual: r.qual,
		pos:  0,
	}
}

// GetChunk implements [SeqRecord].
func (r *FastqSeqRecord) Chunks(length int) iter.Seq[SeqQual] {
	return func(yield func(SeqQual) bool) {
		curPos := 0
		for i := 0; i < len(r.seq); i += length {
			end := min(i+length, len(r.seq))
			if !yield(SeqQual{
				seq:  r.seq[i:end],
				qual: r.qual[i:end],
				name: r.name,
				pos:  curPos,
			}) {
				return
			}
			curPos += (end - i)
		}
	}
}

func (r *FastqSeqRecord) Name() string {
	return r.name
}

func (r *FastqSeqRecord) Comment() string {
	return r.comment
}

// AddCommentTSV appends value to the record's comment, tab-separated.
func (r *FastqSeqRecord) AddCommentTSV(value string) {
	if r.comment == "" {
		r.comment = value
	} else {
		r.comment = r.comment + "\t" + value
	}
}

// FastqWriter writes FASTQ records to a file, optionally gzip-compressed.
type FastqWriter struct {
	writer *bufio.Writer
	gz     *gzip.Writer
	file   *os.File
}

// NewFastqWriter creates a FastqWriter that writes to the given io.Writer.
func NewFastqWriter(w io.Writer) *FastqWriter {
	return &FastqWriter{writer: bufio.NewWriter(w)}
}

// OpenFastqWriter creates a FastqWriter for the given filename.
// If the filename ends in ".gz", the output will be gzip-compressed.
func OpenFastqWriter(filename string) (*FastqWriter, error) {
	f, err := os.Create(filename)
	if err != nil {
		return nil, err
	}
	w := &FastqWriter{file: f}
	if strings.HasSuffix(filename, ".gz") {
		gz := gzip.NewWriter(f)
		w.gz = gz
		w.writer = bufio.NewWriter(gz)
	} else {
		w.writer = bufio.NewWriter(f)
	}
	return w, nil
}

// WriteRecord writes a single FASTQ record to the writer.
func (w *FastqWriter) WriteRecord(name, comment, seq, qual string) error {
	var err error
	if comment != "" {
		_, err = fmt.Fprintf(w.writer, "@%s %s\n%s\n+\n%s\n", name, comment, seq, qual)
	} else {
		_, err = fmt.Fprintf(w.writer, "@%s\n%s\n+\n%s\n", name, seq, qual)
	}
	return err
}

// Write writes a SeqRecord to the file using its Name, Comment, and FullSeq.
func (w *FastqWriter) Write(rec SeqRecord) error {
	sq := rec.FullSeq()
	return w.WriteRecord(rec.Name(), rec.Comment(), sq.Seq(), sq.Qual())
}

// Close flushes and closes the writer.
func (w *FastqWriter) Close() error {
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
