package cram

import (
	"bufio"
	"bytes"
	"fmt"
	"os"
	"strings"

	"github.com/compgen-io/cgkit/seqio"
)

// referenceProvider loads reference sequences for CRAM encoding/decoding.
// Uses seqio.IndexedFastaReader when a .fai index is available (chunk-based,
// LRU-cached). Falls back to loading the full FASTA into memory otherwise.
type referenceProvider struct {
	fastaPath string
	indexed   *seqio.IndexedFastaReader // used when .fai exists
	seqs      map[string][]byte         // fallback: full FASTA in memory
}

// newReferenceProvider creates a reference provider from a FASTA file path.
// Returns an error if the file does not exist.
func newReferenceProvider(fastaPath string) (*referenceProvider, error) {
	if _, err := os.Stat(fastaPath); err != nil {
		return nil, fmt.Errorf("reference FASTA not found: %s", fastaPath)
	}

	rp := &referenceProvider{fastaPath: fastaPath}

	// Try indexed mode (requires .fai).
	if r, err := seqio.NewIndexedFastaReader(fastaPath); err == nil {
		rp.indexed = r
	}
	// If no .fai, fall back to full load on first access.

	return rp, nil
}

// Close releases resources.
func (rp *referenceProvider) Close() error {
	if rp.indexed != nil {
		return rp.indexed.Close()
	}
	return nil
}

// getSequenceRange returns reference bases for [start, end) (0-based).
func (rp *referenceProvider) getSequenceRange(name string, start, end int) ([]byte, error) {
	if rp.indexed != nil {
		return rp.indexed.GetSequenceRange(name, start, end)
	}
	// Fallback: load full FASTA, return slice.
	seq, err := rp.getSequence(name)
	if err != nil {
		return nil, err
	}
	if end > len(seq) {
		end = len(seq)
	}
	if start < 0 {
		start = 0
	}
	if start >= end {
		return nil, nil
	}
	return seq[start:end], nil
}

// getSequence returns the full reference sequence for the given name.
func (rp *referenceProvider) getSequence(name string) ([]byte, error) {
	if rp.indexed != nil {
		return rp.indexed.GetSequence(name)
	}
	// Fallback: load full FASTA.
	if rp.seqs == nil {
		if err := rp.loadFullFasta(); err != nil {
			return nil, err
		}
	}
	seq, ok := rp.seqs[name]
	if !ok {
		return nil, fmt.Errorf("reference %q not found in %s", name, rp.fastaPath)
	}
	return seq, nil
}

// loadFullFasta reads the entire FASTA into memory (fallback when no .fai).
func (rp *referenceProvider) loadFullFasta() error {
	f, err := os.Open(rp.fastaPath)
	if err != nil {
		return fmt.Errorf("opening reference FASTA: %w", err)
	}
	defer f.Close()

	rp.seqs = make(map[string][]byte)
	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 0, 1024*1024), 10*1024*1024)

	var name string
	var seq bytes.Buffer

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			if name != "" {
				rp.seqs[name] = bytes.ToUpper(seq.Bytes())
			}
			fields := strings.Fields(line[1:])
			if len(fields) > 0 {
				name = fields[0]
			}
			seq.Reset()
		} else {
			seq.WriteString(strings.TrimSpace(line))
		}
	}
	if name != "" {
		rp.seqs[name] = bytes.ToUpper(seq.Bytes())
	}

	return scanner.Err()
}
