package cram

import (
	"fmt"
	"os"

	"github.com/compgen-io/cgkit/seqio"
)

// referenceProvider wraps a seqio.ReferenceReader for CRAM encoding/decoding.
type referenceProvider struct {
	ref seqio.ReferenceReader
}

// newReferenceProvider creates a reference provider from a FASTA file path.
// Uses seqio.OpenReference for auto-detection (indexed vs in-memory).
func newReferenceProvider(fastaPath string) (*referenceProvider, error) {
	if _, err := os.Stat(fastaPath); err != nil {
		return nil, fmt.Errorf("reference FASTA not found: %s", fastaPath)
	}

	ref, err := seqio.OpenReference(fastaPath)
	if err != nil {
		return nil, err
	}

	return &referenceProvider{ref: ref}, nil
}

// Close releases resources.
func (rp *referenceProvider) Close() error {
	return rp.ref.Close()
}

// getSequenceRange returns reference bases for [start, end) (0-based).
func (rp *referenceProvider) getSequenceRange(name string, start, end int) ([]byte, error) {
	return rp.ref.GetSequenceRange(name, start, end)
}

// getSequence returns the full reference sequence for the given name.
func (rp *referenceProvider) getSequence(name string) ([]byte, error) {
	return rp.ref.GetSequence(name)
}
