package seqio

import (
	"errors"
	"iter"

	"github.com/compgen-io/cgkit/support/sequtils"
	"github.com/compgen-io/cgkit/support/stringutils"
)

var DirtySeqReaderError = errors.New("input reader is busy")
var ClosedSeqReaderError = errors.New("input reader is closed")

type SeqReader interface {
	NextSeq() (SeqRecord, error)
	Names() (iter.Seq[string], error)
	FetchRecord(name string) (SeqRecord, error)
}

type SeqRecord interface {
	Name() string
	Comment() string

	// GetChunk returns a chunk of the sequence and quality (if available).
	Chunks(length int) iter.Seq[SeqQual]
	FullSeq() SeqQual
}

type SeqQual struct {
	seq     string
	qual    string
	name    string
	pos     int
	revcomp bool
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

func (s SeqQual) Name() string {
	return s.name
}

func (s SeqQual) Position() int {
	return s.pos
}

func (s SeqQual) IsRevComp() bool {
	return s.revcomp
}

func (s SeqQual) Strand() string {
	if s.revcomp {
		return "-"
	}
	return "+"
}

func (s SeqQual) RevComp() SeqQual {
	return SeqQual{
		seq:     sequtils.ReverseCompliment(s.seq),
		qual:    stringutils.ReverseString(s.qual),
		name:    s.name,
		pos:     s.pos,
		revcomp: !s.revcomp,
	}
}

func NewStringSeq(seq string, namecomment ...string) SeqRecord {
	s := &stringSeq{seq: seq}
	if len(namecomment) > 0 {
		s.name = namecomment[0]
	}
	if len(namecomment) > 1 {
		s.comment = namecomment[1]
	}
	return s
}

func NewStringSeqQual(seq string, qual string, namecomment ...string) SeqRecord {
	s := &stringSeq{seq: seq, qual: qual}
	if len(namecomment) > 0 {
		s.name = namecomment[0]
	}
	if len(namecomment) > 1 {
		s.comment = namecomment[1]
	}
	return s
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
		name: s.name,
		pos:  0,
	}
}

func (s *stringSeq) Chunks(n int) iter.Seq[SeqQual] {
	return func(yield func(SeqQual) bool) {
		curPos := 0
		// Clamp to the shorter of seq/qual to avoid slicing panics.
		total := min(len(s.qual), len(s.seq))
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
				name: s.name,
				pos:  curPos,
			}
			curPos += (end - i)

			if !yield(chunk) {
				return
			}
		}
	}
}

func (s *SeqQual) Sub(start, end int) SeqQual {
	return SeqQual{
		name:    s.name,
		seq:     s.seq[start:end],
		qual:    s.qual[start:end],
		pos:     s.pos + start,
		revcomp: s.revcomp,
	}
}
