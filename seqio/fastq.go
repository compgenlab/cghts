package seqio

import (
	"bufio"
	"compress/gzip"
	"io"
	"iter"
	"os"
	"strings"
)

type FastqSeqReader struct {
	buf    string
	reader *bufio.Reader
}

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

func (r *FastqReader) NextSeq() (SeqRecord, error) {
	if r.closed {
		return nil, ClosedSeqReaderError
	}
	for {
		if b, err := r.reader.ReadByte(); err != nil {
			return nil, err
		} else if b == '@' && r.lastByte == '\n' {

			// TODO: make this support multi-line FASTQ records.
			//       Right now it only supports single-line sequences
			//       and quality strings, which is the most common format,
			//       but we should support multi-line records.

			// start of a FASTQ record, read each line
			header, err := r.reader.ReadString('\n')
			if err != nil {
				return nil, err
			}
			seq, err := r.reader.ReadString('\n')
			if err != nil {
				return nil, err
			}
			plus, err := r.reader.ReadString('\n')
			if err != nil {
				return nil, err
			}
			if plus[0] != '+' {
				return nil, io.ErrUnexpectedEOF
			}
			qual, err := r.reader.ReadString('\n')
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
		seq:  r.seq,
		qual: r.qual,
	}
}

// GetChunk implements [SeqRecord].
func (r *FastqSeqRecord) Chunks(length int) iter.Seq[SeqQual] {
	return func(yield func(SeqQual) bool) {
		for i := 0; i < len(r.seq); i += length {
			end := i + length
			if end > len(r.seq) {
				end = len(r.seq)
			}
			if !yield(SeqQual{
				seq:  r.seq[i:end],
				qual: r.qual[i:end],
			}) {
				return
			}
		}
	}
}

func (r *FastqSeqRecord) Name() string {
	return r.name
}

func (r *FastqSeqRecord) Comment() string {
	return r.comment
}
