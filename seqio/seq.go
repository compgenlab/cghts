package seqio

import (
	"errors"
	"iter"

	"github.com/compgenlab/hts/support/sequtils"
	"github.com/compgenlab/hts/support/stringutils"
)

// DirtySeqReaderError is returned when an operation is attempted on a reader
// that is still busy emitting a previous record.
var DirtySeqReaderError = errors.New("input reader is busy")

// ClosedSeqReaderError is returned by reader methods after the reader has been
// closed (including the automatic close performed by Names and FetchRecord
// once they exhaust the non-seekable underlying stream).
var ClosedSeqReaderError = errors.New("input reader is closed")

// SeqReader is the common interface for streaming sequence sources such as
// [FastaReader] and [FastqReader]. Readers are forward-only and not seekable.
type SeqReader interface {
	// NextSeq returns the next record, or io.EOF when the stream is exhausted.
	NextSeq() (SeqRecord, error)
	// Names returns an iterator over record names. Because the underlying
	// stream cannot be rewound, iterating consumes it and closes the reader.
	Names() (iter.Seq[string], error)
	// FetchRecord scans forward for a record with the given name. As it cannot
	// rewind, it consumes the stream and closes the reader; returns io.EOF if
	// no matching record is found.
	FetchRecord(name string) (SeqRecord, error)
}

// SeqRecord represents a single sequence entry (one FASTA/FASTQ record).
type SeqRecord interface {
	// Name returns the record identifier (the first whitespace-delimited token
	// of the header).
	Name() string
	// Comment returns any text following the name on the header line.
	Comment() string

	// Chunks yields the sequence (and quality, if present) in pieces of at
	// most length bases. For streaming FASTA records this avoids loading the
	// entire sequence into memory at once.
	Chunks(length int) iter.Seq[SeqQual]
	// FullSeq returns the entire record as a single [SeqQual].
	FullSeq() SeqQual
}

// SeqQual holds a sequence with an optional quality string, plus the record
// name, a position offset, and strand orientation. SeqQual values are
// immutable; methods that transform them return new values.
type SeqQual struct {
	seq     string
	qual    string
	name    string
	pos     int
	revcomp bool
}

// Len returns the length of the sequence in bases.
func (s SeqQual) Len() int {
	return len(s.seq)
}

// Seq returns the sequence string.
func (s SeqQual) Seq() string {
	return s.seq
}

// Qual returns the quality string, or "" if no quality is associated with the
// sequence (e.g. for FASTA records).
func (s SeqQual) Qual() string {
	return s.qual
}

// Name returns the record name.
func (s SeqQual) Name() string {
	return s.name
}

// Position returns the 0-based offset of this sequence within its parent
// record. For a full record it is 0; for chunks it reflects the running offset.
func (s SeqQual) Position() int {
	return s.pos
}

// IsRevComp reports whether this sequence has been reverse-complemented
// relative to the original record.
func (s SeqQual) IsRevComp() bool {
	return s.revcomp
}

// Strand returns "-" if the sequence is reverse-complemented, otherwise "+".
func (s SeqQual) Strand() string {
	if s.revcomp {
		return "-"
	}
	return "+"
}

// RevComp returns a new [SeqQual] that is the reverse complement of the
// sequence, with the quality string reversed and the strand flag toggled. The
// name and position are preserved. The receiver is not modified.
func (s SeqQual) RevComp() SeqQual {
	return SeqQual{
		seq:     sequtils.ReverseComplement(s.seq),
		qual:    stringutils.ReverseString(s.qual),
		name:    s.name,
		pos:     s.pos,
		revcomp: !s.revcomp,
	}
}

// NewStringSeq returns a [SeqRecord] backed by an in-memory sequence string
// with no quality. The optional namecomment arguments set the name (first) and
// comment (second).
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

// NewStringSeqQual returns a [SeqRecord] backed by in-memory sequence and
// quality strings. The optional namecomment arguments set the name (first) and
// comment (second).
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

// Sub returns the subsequence covering [start, end) (0-based, half-open) of the
// receiver, slicing both the sequence and quality strings. The returned
// position is offset by start; the strand flag is preserved. It panics if the
// bounds are out of range. The receiver is not modified.
func (s SeqQual) Sub(start, end int) SeqQual {
	return SeqQual{
		name:    s.name,
		seq:     s.seq[start:end],
		qual:    s.qual[start:end],
		pos:     s.pos + start,
		revcomp: s.revcomp,
	}
}
