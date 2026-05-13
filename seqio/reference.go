package seqio

import (
	"bufio"
	"bytes"
	"fmt"
	"os"
	"strings"
)

// ReferenceReader provides random access to named reference sequences.
// This is the interface for reference genomes — few long sequences accessed
// by name and position range. Implementations include indexed FASTA (local,
// plain or bgzip) and in-memory FASTA.
type ReferenceReader interface {
	// Names returns the ordered sequence names.
	Names() []string

	// SequenceLength returns the length of the named sequence.
	// Returns 0, false if the sequence is not found.
	SequenceLength(name string) (int, bool)

	// GetSequenceRange returns bases for [start, end) (0-based) of the
	// named sequence. The returned slice is uppercased. Out-of-range
	// coordinates are clamped to the sequence boundaries.
	GetSequenceRange(name string, start, end int) ([]byte, error)

	// GetSequence returns the full sequence for the given name, uppercased.
	// Prefer GetSequenceRange when you know the range you need.
	GetSequence(name string) ([]byte, error)

	// Close releases resources held by the reader.
	Close() error
}

// OpenReference opens a reference sequence file, auto-detecting the format.
//
// Supported formats:
//   - Indexed FASTA (.fa/.fasta with .fai index, plain or bgzip-compressed)
//   - Plain FASTA without index (loaded fully into memory)
//
// For indexed FASTA, sequences are loaded in 10MB chunks with LRU caching.
// For bgzip-compressed FASTA (.fa.gz), a .gzi index is loaded automatically.
func OpenReference(path string) (ReferenceReader, error) {
	if _, err := os.Stat(path); err != nil {
		return nil, fmt.Errorf("reference not found: %s", path)
	}

	// Try indexed mode first (.fai must exist).
	r, err := NewIndexedFastaReader(path)
	if err == nil {
		return r, nil
	}

	// Fall back to in-memory FASTA.
	return newInMemoryFasta(path)
}

// inMemoryFasta loads the entire FASTA into memory. Used as a fallback
// when no .fai index is available (e.g., small test references).
type inMemoryFasta struct {
	path  string
	names []string
	seqs  map[string][]byte
}

func newInMemoryFasta(path string) (*inMemoryFasta, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("opening reference FASTA: %w", err)
	}
	defer f.Close()

	r := &inMemoryFasta{
		path: path,
		seqs: make(map[string][]byte),
	}

	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 0, 1024*1024), 10*1024*1024)

	var name string
	var seq bytes.Buffer

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			if name != "" {
				r.seqs[name] = bytes.ToUpper(seq.Bytes())
			}
			fields := strings.Fields(line[1:])
			if len(fields) > 0 {
				name = fields[0]
			}
			r.names = append(r.names, name)
			seq.Reset()
		} else {
			seq.WriteString(strings.TrimSpace(line))
		}
	}
	if name != "" {
		r.seqs[name] = bytes.ToUpper(seq.Bytes())
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return r, nil
}

func (r *inMemoryFasta) Names() []string { return r.names }

func (r *inMemoryFasta) SequenceLength(name string) (int, bool) {
	seq, ok := r.seqs[name]
	if !ok {
		return 0, false
	}
	return len(seq), true
}

func (r *inMemoryFasta) GetSequenceRange(name string, start, end int) ([]byte, error) {
	seq, ok := r.seqs[name]
	if !ok {
		return nil, fmt.Errorf("reference %q not found in %s", name, r.path)
	}
	if start < 0 {
		start = 0
	}
	if end > len(seq) {
		end = len(seq)
	}
	if start >= end {
		return nil, nil
	}
	// Return a copy to avoid callers mutating internal state.
	out := make([]byte, end-start)
	copy(out, seq[start:end])
	return out, nil
}

func (r *inMemoryFasta) GetSequence(name string) ([]byte, error) {
	seq, ok := r.seqs[name]
	if !ok {
		return nil, fmt.Errorf("reference %q not found in %s", name, r.path)
	}
	out := make([]byte, len(seq))
	copy(out, seq)
	return out, nil
}

func (r *inMemoryFasta) Close() error { return nil }
