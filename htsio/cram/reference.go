package cram

import (
	"bufio"
	"bytes"
	"fmt"
	"os"
	"strings"
)

// referenceProvider loads reference sequences for CRAM decoding.
type referenceProvider struct {
	fastaPath string
	seqs      map[string][]byte // name → uppercase sequence
}

// newReferenceProvider creates a reference provider from a FASTA file path.
// The FASTA is loaded lazily on first access.
func newReferenceProvider(fastaPath string) *referenceProvider {
	return &referenceProvider{
		fastaPath: fastaPath,
	}
}

// getSequence returns the reference sequence for the given name, uppercased.
func (rp *referenceProvider) getSequence(name string) ([]byte, error) {
	if rp.seqs == nil {
		if err := rp.load(); err != nil {
			return nil, err
		}
	}
	seq, ok := rp.seqs[name]
	if !ok {
		return nil, fmt.Errorf("reference %q not found in %s", name, rp.fastaPath)
	}
	return seq, nil
}

func (rp *referenceProvider) load() error {
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
			// Parse name: first word after >
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
