package seqio

import (
	"io"
	"strings"
	"testing"
)

// Both readers are used through the SeqReader interface, so both must end it
// the same way. FastqReader.NextSeq used to `return r.NextFastqSeq()` directly,
// and that returns a *FastqSeqRecord -- nil at EOF, but boxed into a SeqRecord
// it becomes a non-nil interface holding a nil pointer. A caller checking
// `rec == nil` to detect the end sailed past it and dereferenced the nil on the
// next method call, which crashed cgkit's fastq-gc on every successful run.
func TestNextSeqReturnsATrueNilAtEOF(t *testing.T) {
	cases := []struct {
		name string
		open func() (SeqReader, error)
	}{
		{"fasta", func() (SeqReader, error) {
			return NewFastaReader(strings.NewReader(">r1\nACGT\n"))
		}},
		{"fastq", func() (SeqReader, error) {
			return NewFastqReader(strings.NewReader("@r1\nACGT\n+\nIIII\n"))
		}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			r, err := c.open()
			if err != nil {
				t.Fatal(err)
			}

			// Drain the one real record.
			if _, err := r.NextSeq(); err != nil {
				t.Fatalf("first record: %v", err)
			}

			rec, err := r.NextSeq()
			if err != io.EOF {
				t.Errorf("err = %v, want io.EOF", err)
			}
			// The interface itself must be nil, not merely hold a nil pointer.
			if rec != nil {
				t.Fatalf("NextSeq at EOF returned a non-nil SeqRecord (%T); a caller "+
					"testing rec == nil would dereference it", rec)
			}
		})
	}
}
