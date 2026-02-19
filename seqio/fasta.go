package seqio

import (
	"bufio"
	"compress/gzip"
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

	// TODO: handle FAIDX indexed files (by returning an IndexedFastaReader)

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
				name:    name,
				comment: comment,
				reader:  r.reader,
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
		buf.WriteByte(b)
		last = b
	}
	return SeqQual{seq: buf.String()}
}

// Chunks implements [SeqRecord].
func (r *FastaSeqRecord) Chunks(length int) iter.Seq[SeqQual] {
	return func(yield func(SeqQual) bool) {
		peek, _ := r.reader.Peek(length)
		if len(peek) == 0 {
			return
		}
		var last byte
		blen := length
		for i, b := range peek {
			if last == '\n' && b == '>' {
				blen = i
				break
			}
			last = b
		}
		buf := make([]byte, blen)
		for {
			n, err := r.reader.Read(buf)
			if n > 0 {
				chunk := string(buf[:n])
				chunk = strings.ReplaceAll(chunk, "\n", "")
				chunk = strings.ReplaceAll(chunk, "\r", "")
				if !yield(SeqQual{seq: chunk}) {
					return
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
