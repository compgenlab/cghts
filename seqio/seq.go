package seqio

import (
	"errors"
	"iter"
)

var DirtySeqReaderError = errors.New("input reader is busy")
var ClosedSeqReaderError = errors.New("input reader is closed")

type SeqRecord interface {
	Name() string
	Comment() string

	// GetChunk returns a chunk of the sequence and quality (if available).
	Chunks(length int) iter.Seq[SeqQual]
	FullSeq() SeqQual
}

type SeqQual struct {
	seq  string
	qual string
}

type SeqReader interface {
	NextSeq() (SeqRecord, error)
	Names() (iter.Seq[string], error)
	FetchRecord(name string) (SeqRecord, error)
}

func (s SeqQual) Len() int {
	return len(s.seq)
}

func (s SeqQual) Seq() string {
	return s.seq
}

func (s SeqQual) Qual() string {
	return s.qual
}

func NewStringSeq(seq string, namecomment ...string) (SeqRecord, error) {
	s := &stringSeq{seq: seq}
	if len(namecomment) > 0 {
		s.name = namecomment[0]
	}
	if len(namecomment) > 1 {
		s.comment = namecomment[1]
	}
	return s, nil
}

func NewStringSeqQual(seq string, qual string, namecomment ...string) (SeqRecord, error) {
	s := &stringSeq{seq: seq, qual: qual}
	if len(namecomment) > 0 {
		s.name = namecomment[0]
	}
	if len(namecomment) > 1 {
		s.comment = namecomment[1]
	}
	return s, nil
}

type stringSeq struct {
	name    string
	comment string
	seq     string
	qual    string
}

func (s *stringSeq) Name() string {
	return s.name
}

func (s *stringSeq) Comment() string {
	return s.comment
}

func (s *stringSeq) FullSeq() SeqQual {
	return SeqQual{
		seq:  s.seq,
		qual: s.qual,
	}
}

func (s *stringSeq) Chunks(n int) iter.Seq[SeqQual] {
	return func(yield func(SeqQual) bool) {
		// Clamp to the shorter of seq/qual to avoid slicing panics.
		total := len(s.seq)
		if len(s.qual) < total {
			total = len(s.qual)
		}
		if total == 0 {
			return
		}

		if n <= 0 || n > total {
			n = total
		}

		for i := 0; i < total; i += n {
			end := i + n
			if end > total {
				end = total
			}

			chunk := SeqQual{
				seq:  s.seq[i:end],
				qual: s.qual[i:end],
			}

			if !yield(chunk) {
				return
			}
		}
	}
}
