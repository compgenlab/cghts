package cram

import (
	"bytes"
	"fmt"
	"io"
)

// compressionHeader holds the parsed compression header for a CRAM container.
type compressionHeader struct {
	// Preservation map
	readNamesPreserved bool
	apDelta            bool
	refRequired        bool
	substitutionMatrix [5][4]byte // for each ref base (A,C,G,T,N), 4 substitution bases
	tagDictionary      [][]tagKey // list of tag combos; each combo is a list of tag keys

	// Data series encoding map
	dataSeriesEncodings map[string]dataCodec

	// Tag encoding map
	tagEncodings map[int32]dataCodec
}

// tagKey represents a tag ID (2 chars) and its type (1 char), e.g. "RG" + 'Z'.
type tagKey struct {
	id   [2]byte
	typ  byte
}

// tagKeyToITF8 converts a tagKey to the 3-byte ITF8 key used in the tag encoding map.
func tagKeyToITF8(tk tagKey) int32 {
	return (int32(tk.id[0]) << 16) | (int32(tk.id[1]) << 8) | int32(tk.typ)
}

func readCompressionHeader(data []byte) (*compressionHeader, error) {
	r := bytes.NewReader(data)
	ch := &compressionHeader{
		dataSeriesEncodings: make(map[string]dataCodec),
		tagEncodings:        make(map[int32]dataCodec),
	}

	// Preservation map
	if err := ch.readPreservationMap(r); err != nil {
		return nil, fmt.Errorf("preservation map: %w", err)
	}

	// Data series encoding map
	if err := ch.readDataSeriesMap(r); err != nil {
		return nil, fmt.Errorf("data series map: %w", err)
	}

	// Tag encoding map
	if err := ch.readTagEncodingMap(r); err != nil {
		return nil, fmt.Errorf("tag encoding map: %w", err)
	}

	return ch, nil
}

func (ch *compressionHeader) readPreservationMap(r io.Reader) error {
	mapSize, err := readITF8(r)
	if err != nil {
		return err
	}
	// Read the map content from the known size
	mapData := make([]byte, mapSize)
	if _, err := io.ReadFull(r, mapData); err != nil {
		return err
	}
	mr := bytes.NewReader(mapData)

	numEntries, err := readITF8(mr)
	if err != nil {
		return err
	}

	// Defaults
	ch.readNamesPreserved = true
	ch.apDelta = true
	ch.refRequired = true

	for i := int32(0); i < numEntries; i++ {
		var key [2]byte
		if _, err := io.ReadFull(mr, key[:]); err != nil {
			return fmt.Errorf("reading preservation key: %w", err)
		}

		switch string(key[:]) {
		case "RN":
			var v [1]byte
			if _, err := io.ReadFull(mr, v[:]); err != nil {
				return err
			}
			ch.readNamesPreserved = v[0] != 0
		case "AP":
			var v [1]byte
			if _, err := io.ReadFull(mr, v[:]); err != nil {
				return err
			}
			ch.apDelta = v[0] != 0
		case "RR":
			var v [1]byte
			if _, err := io.ReadFull(mr, v[:]); err != nil {
				return err
			}
			ch.refRequired = v[0] != 0
		case "SM":
			var sm [5]byte
			if _, err := io.ReadFull(mr, sm[:]); err != nil {
				return err
			}
			ch.parseSubstitutionMatrix(sm)
		case "TD":
			// Tag dictionary: stored as a byte array
			tdLen, err := readITF8(mr)
			if err != nil {
				return err
			}
			tdData := make([]byte, tdLen)
			if _, err := io.ReadFull(mr, tdData); err != nil {
				return err
			}
			ch.parseTagDictionary(tdData)
		default:
			// Unknown key — skip one byte value
			var v [1]byte
			io.ReadFull(mr, v[:])
		}
	}

	return nil
}

// parseSubstitutionMatrix decodes the 5-byte substitution matrix.
// Each byte encodes 4 substitution bases (2 bits each) for ref base A, C, G, T, N.
func (ch *compressionHeader) parseSubstitutionMatrix(sm [5]byte) {
	// Each byte encodes 4 substitution base assignments for one reference base.
	// The 2-bit fields specify which substitution code (0-3) each alternate base
	// is assigned to. For ref A, others are C,G,T,N in that order.
	// htslib: matrix[refIdx][(byte>>6)&3] = first_other, etc.
	bases := [5]byte{'A', 'C', 'G', 'T', 'N'}
	for i := 0; i < 5; i++ {
		others := otherBases(bases[i])
		b := sm[i]
		ch.substitutionMatrix[i][(b>>6)&3] = others[0]
		ch.substitutionMatrix[i][(b>>4)&3] = others[1]
		ch.substitutionMatrix[i][(b>>2)&3] = others[2]
		ch.substitutionMatrix[i][(b>>0)&3] = others[3]
	}
}

// otherBases returns the 4 bases that are not refBase, in order A,C,G,T,N.
func otherBases(refBase byte) [4]byte {
	all := [5]byte{'A', 'C', 'G', 'T', 'N'}
	var result [4]byte
	j := 0
	for _, b := range all {
		if b != refBase {
			result[j] = b
			j++
		}
	}
	return result
}

// parseTagDictionary parses the TD preservation map value.
// Format: sequences of 3-byte tag keys (ID + type) separated by NUL bytes.
// Each NUL-terminated group is one tag combination.
func (ch *compressionHeader) parseTagDictionary(data []byte) {
	var current []tagKey
	pos := 0
	for pos < len(data) {
		if data[pos] == 0 {
			ch.tagDictionary = append(ch.tagDictionary, current)
			current = nil
			pos++
			continue
		}
		if pos+3 > len(data) {
			break
		}
		tk := tagKey{
			id:  [2]byte{data[pos], data[pos+1]},
			typ: data[pos+2],
		}
		current = append(current, tk)
		pos += 3
	}
	// If there's a trailing combo without NUL terminator
	if len(current) > 0 {
		ch.tagDictionary = append(ch.tagDictionary, current)
	}
}

func (ch *compressionHeader) readDataSeriesMap(r io.Reader) error {
	mapSize, err := readITF8(r)
	if err != nil {
		return err
	}
	mapData := make([]byte, mapSize)
	if _, err := io.ReadFull(r, mapData); err != nil {
		return err
	}
	mr := bytes.NewReader(mapData)

	numEntries, err := readITF8(mr)
	if err != nil {
		return err
	}

	for i := int32(0); i < numEntries; i++ {
		// Key: 2-byte data series name
		var key [2]byte
		if _, err := io.ReadFull(mr, key[:]); err != nil {
			return fmt.Errorf("reading data series key: %w", err)
		}
		name := string(key[:])

		enc, err := readEncoding(mr)
		if err != nil {
			return fmt.Errorf("reading encoding for %s: %w", name, err)
		}

		ch.dataSeriesEncodings[name] = enc
	}

	return nil
}

func (ch *compressionHeader) readTagEncodingMap(r io.Reader) error {
	mapSize, err := readITF8(r)
	if err != nil {
		return err
	}
	mapData := make([]byte, mapSize)
	if _, err := io.ReadFull(r, mapData); err != nil {
		return err
	}
	mr := bytes.NewReader(mapData)

	numEntries, err := readITF8(mr)
	if err != nil {
		return err
	}

	for i := int32(0); i < numEntries; i++ {
		key, err := readITF8(mr)
		if err != nil {
			return fmt.Errorf("reading tag key: %w", err)
		}

		enc, err := readEncoding(mr)
		if err != nil {
			return fmt.Errorf("reading tag encoding for %06x: %w", key, err)
		}

		ch.tagEncodings[key] = enc
	}

	return nil
}

// getIntCodec returns the integer codec for a data series, or an error if not found.
func (ch *compressionHeader) getIntCodec(name string) (intCodec, error) {
	c, ok := ch.dataSeriesEncodings[name]
	if !ok {
		return nil, fmt.Errorf("data series %q not found", name)
	}
	ic, ok := c.(intCodec)
	if !ok {
		return nil, fmt.Errorf("data series %q is not an int codec", name)
	}
	return ic, nil
}

// getByteCodec returns the byte codec for a data series, or an error if not found.
func (ch *compressionHeader) getByteCodec(name string) (byteCodec, error) {
	c, ok := ch.dataSeriesEncodings[name]
	if !ok {
		return nil, fmt.Errorf("data series %q not found", name)
	}
	bc, ok := c.(byteCodec)
	if !ok {
		return nil, fmt.Errorf("data series %q is not a byte codec", name)
	}
	return bc, nil
}

// getByteArrayCodec returns the byte array codec for a data series, or an error if not found.
func (ch *compressionHeader) getByteArrayCodec(name string) (byteArrayCodec, error) {
	c, ok := ch.dataSeriesEncodings[name]
	if !ok {
		return nil, fmt.Errorf("data series %q not found", name)
	}
	bac, ok := c.(byteArrayCodec)
	if !ok {
		return nil, fmt.Errorf("data series %q is not a byte array codec", name)
	}
	return bac, nil
}

// lookupSubstitution returns the substituted base given a reference base and substitution code.
func (ch *compressionHeader) lookupSubstitution(refBase byte, code byte) byte {
	var idx int
	switch refBase {
	case 'A':
		idx = 0
	case 'C':
		idx = 1
	case 'G':
		idx = 2
	case 'T':
		idx = 3
	case 'N':
		idx = 4
	default:
		idx = 4
	}
	if int(code) >= len(ch.substitutionMatrix[idx]) {
		return 'N'
	}
	return ch.substitutionMatrix[idx][code]
}
