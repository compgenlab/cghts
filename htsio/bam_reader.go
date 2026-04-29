package htsio

import (
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"

	"github.com/compgen-io/cgltk/htsio/bgzf"
)

// BAM CIGAR operation codes (4-bit).
const (
	bamCigarM    = 0 // alignment match (can be seq match or mismatch)
	bamCigarI    = 1 // insertion to the reference
	bamCigarD    = 2 // deletion from the reference
	bamCigarN    = 3 // skipped region from the reference
	bamCigarS    = 4 // soft clipping
	bamCigarH    = 5 // hard clipping
	bamCigarP    = 6 // padding
	bamCigarEq   = 7 // sequence match
	bamCigarX    = 8 // sequence mismatch
	bamCigarBack = 9 // reserved/back (unused in practice)
)

var cigarOpChar = [10]byte{'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B'}

// seqDecode maps 4-bit BAM sequence codes to ASCII bases.
var seqDecode = [16]byte{'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'}

// bamRefInfo holds a reference sequence name and length from the BAM header.
type bamRefInfo struct {
	name   string
	length int32
}

// BamReader reads BAM files natively using the bgzf package.
// It implements the SamReader interface.
type BamReader struct {
	r    *bgzf.Reader
	src  io.ReadCloser // underlying file/reader, closed on Close()
	refs []bamRefInfo  // reference sequences from the BAM header
	hdr  *SamHeader
	opts *SamReaderOpts
	err  error // sticky error
}

// NewBamReader creates a BAM reader from an io.ReadCloser.
// The reader must be positioned at the start of a BAM file.
func NewBamReader(rc io.ReadCloser, opts ...*SamReaderOpts) (*BamReader, error) {
	var o *SamReaderOpts
	if len(opts) > 0 {
		o = opts[0]
	} else {
		o = NewSamReaderOpts()
	}

	br := &BamReader{
		r:    bgzf.NewReader(rc),
		src:  rc,
		opts: o,
	}

	if err := br.readHeader(); err != nil {
		rc.Close()
		return nil, fmt.Errorf("bam: reading header: %w", err)
	}

	return br, nil
}

// Header returns the parsed SAM header.
func (b *BamReader) Header() (*SamHeader, error) {
	return b.hdr, nil
}

// Close releases resources.
func (b *BamReader) Close() error {
	if b.src != nil {
		return b.src.Close()
	}
	return nil
}

// readHeader reads the BAM magic, header text, and reference sequence dictionary.
func (b *BamReader) readHeader() error {
	// Magic: BAM\1
	var magic [4]byte
	if _, err := io.ReadFull(b.r, magic[:]); err != nil {
		return fmt.Errorf("reading magic: %w", err)
	}
	if magic != [4]byte{'B', 'A', 'M', 1} {
		return fmt.Errorf("invalid BAM magic: %x", magic)
	}

	// Header text length + text.
	var headerLen int32
	if err := binary.Read(b.r, binary.LittleEndian, &headerLen); err != nil {
		return fmt.Errorf("reading header length: %w", err)
	}
	if headerLen < 0 {
		return fmt.Errorf("negative header length: %d", headerLen)
	}

	headerText := make([]byte, headerLen)
	if _, err := io.ReadFull(b.r, headerText); err != nil {
		return fmt.Errorf("reading header text: %w", err)
	}

	b.hdr = NewSamHeader()
	for _, line := range strings.Split(string(headerText), "\n") {
		line = strings.TrimRight(line, "\r")
		if line != "" {
			b.hdr.AddLine(line)
		}
	}

	// Number of reference sequences.
	var nRef int32
	if err := binary.Read(b.r, binary.LittleEndian, &nRef); err != nil {
		return fmt.Errorf("reading ref count: %w", err)
	}
	if nRef < 0 {
		return fmt.Errorf("negative ref count: %d", nRef)
	}

	b.refs = make([]bamRefInfo, nRef)
	for i := int32(0); i < nRef; i++ {
		var nameLen int32
		if err := binary.Read(b.r, binary.LittleEndian, &nameLen); err != nil {
			return fmt.Errorf("reading ref name length [%d]: %w", i, err)
		}
		nameBuf := make([]byte, nameLen)
		if _, err := io.ReadFull(b.r, nameBuf); err != nil {
			return fmt.Errorf("reading ref name [%d]: %w", i, err)
		}
		// Name is NUL-terminated.
		name := string(nameBuf[:nameLen-1])

		var refLen int32
		if err := binary.Read(b.r, binary.LittleEndian, &refLen); err != nil {
			return fmt.Errorf("reading ref length [%d]: %w", i, err)
		}
		b.refs[i] = bamRefInfo{name: name, length: refLen}
	}

	return nil
}

// Next returns the next SamRecord. Returns nil, io.EOF when done.
func (b *BamReader) Next() (*SamRecord, error) {
	if b.err != nil {
		return nil, b.err
	}

	for {
		rec, err := b.readRecord()
		if err != nil {
			b.err = err
			return nil, err
		}
		if b.passesFilters(rec) {
			return rec, nil
		}
	}
}

// passesFilters checks flag, mapq, and tag filters.
func (b *BamReader) passesFilters(rec *SamRecord) bool {
	if b.opts.flagReq != 0 && rec.Flag&b.opts.flagReq != b.opts.flagReq {
		return false
	}
	if b.opts.flagFilter != 0 && rec.Flag&b.opts.flagFilter != 0 {
		return false
	}
	if b.opts.minMapQ != 0 && rec.MapQ < b.opts.minMapQ {
		return false
	}
	for _, f := range b.opts.tagFilters {
		if !f.matchesRecord(rec) {
			return false
		}
	}
	return true
}

// readRecord is a convenience wrapper for the BamReader.
func (b *BamReader) readRecord() (*SamRecord, error) {
	return readBamRecord(b.r, b.refs)
}

// readBamRecord reads a single BAM alignment record from r and converts it
// to a SamRecord, using refs to resolve reference IDs to names.
func readBamRecord(r io.Reader, refs []bamRefInfo) (*SamRecord, error) {
	// Block size (excludes the 4-byte block_size field itself).
	var blockSize int32
	if err := binary.Read(r, binary.LittleEndian, &blockSize); err != nil {
		if err == io.EOF || err == io.ErrUnexpectedEOF {
			return nil, io.EOF
		}
		return nil, fmt.Errorf("reading block size: %w", err)
	}
	if blockSize < 32 {
		return nil, fmt.Errorf("bam: record block too small: %d", blockSize)
	}

	// Read the entire record into a buffer for efficient parsing.
	buf := make([]byte, blockSize)
	if _, err := io.ReadFull(r, buf); err != nil {
		if err == io.ErrUnexpectedEOF {
			return nil, fmt.Errorf("bam: truncated record")
		}
		return nil, fmt.Errorf("reading record: %w", err)
	}

	// Fixed-length fields (32 bytes).
	refID := int32(binary.LittleEndian.Uint32(buf[0:4]))
	pos := int32(binary.LittleEndian.Uint32(buf[4:8]))
	nameLen := buf[8] // l_read_name (includes NUL)
	mapq := buf[9]
	// bin := binary.LittleEndian.Uint16(buf[10:12]) // unused for now
	nCigarOps := binary.LittleEndian.Uint16(buf[12:14])
	flag := binary.LittleEndian.Uint16(buf[14:16])
	seqLen := int32(binary.LittleEndian.Uint32(buf[16:20]))
	nextRefID := int32(binary.LittleEndian.Uint32(buf[20:24]))
	nextPos := int32(binary.LittleEndian.Uint32(buf[24:28]))
	tlen := int32(binary.LittleEndian.Uint32(buf[28:32]))

	offset := 32

	// Read name (NUL-terminated).
	if offset+int(nameLen) > len(buf) {
		return nil, fmt.Errorf("bam: record truncated at read name")
	}
	readName := string(buf[offset : offset+int(nameLen)-1]) // exclude NUL
	offset += int(nameLen)

	// CIGAR.
	cigarBytes := int(nCigarOps) * 4
	if offset+cigarBytes > len(buf) {
		return nil, fmt.Errorf("bam: record truncated at cigar")
	}
	cigar := decodeCigar(buf[offset:offset+cigarBytes], nCigarOps)
	offset += cigarBytes

	// Sequence (4-bit encoded, 2 bases per byte).
	seqBytes := (int(seqLen) + 1) / 2
	if offset+seqBytes > len(buf) {
		return nil, fmt.Errorf("bam: record truncated at seq")
	}
	seq := decodeSeq(buf[offset:offset+seqBytes], int(seqLen))
	offset += seqBytes

	// Quality (Phred+33 in SAM, raw Phred in BAM).
	if offset+int(seqLen) > len(buf) {
		return nil, fmt.Errorf("bam: record truncated at qual")
	}
	qual := decodeQual(buf[offset:offset+int(seqLen)])
	offset += int(seqLen)

	// Auxiliary tags.
	tags := decodeTags(buf[offset:])

	// Resolve reference names.
	refName := "*"
	if refID >= 0 && int(refID) < len(refs) {
		refName = refs[refID].name
	}
	refNext := "*"
	if nextRefID >= 0 && int(nextRefID) < len(refs) {
		refNext = refs[nextRefID].name
	}
	// SAM uses "=" for mate on same reference.
	if nextRefID == refID && refID >= 0 {
		refNext = "="
	}

	rec := &SamRecord{
		ReadName:  readName,
		Flag:      int(flag),
		RefName:   refName,
		Pos:       int(pos) + 1, // BAM is 0-based, SAM is 1-based
		MapQ:      int(mapq),
		Cigar:     cigar,
		RefNext:   refNext,
		PosNext:   int(nextPos) + 1,
		InsertLen: int(tlen),
		Seq:       seq,
		Qual:      qual,
		Tags:      tags,
	}

	return rec, nil
}

// decodeCigar converts packed BAM CIGAR operations to a SAM CIGAR string.
func decodeCigar(data []byte, nOps uint16) string {
	if nOps == 0 {
		return "*"
	}
	var sb strings.Builder
	sb.Grow(int(nOps) * 3) // rough estimate
	for i := 0; i < int(nOps); i++ {
		packed := binary.LittleEndian.Uint32(data[i*4:])
		opLen := packed >> 4
		opCode := packed & 0xf
		if opCode > 9 {
			opCode = 0 // treat unknown as M
		}
		sb.WriteString(strconv.FormatUint(uint64(opLen), 10))
		sb.WriteByte(cigarOpChar[opCode])
	}
	return sb.String()
}

// decodeSeq converts 4-bit packed BAM sequence to a string.
func decodeSeq(data []byte, seqLen int) string {
	if seqLen == 0 {
		return "*"
	}
	buf := make([]byte, seqLen)
	for i := 0; i < seqLen; i++ {
		if i%2 == 0 {
			buf[i] = seqDecode[data[i/2]>>4]
		} else {
			buf[i] = seqDecode[data[i/2]&0xf]
		}
	}
	return string(buf)
}

// decodeQual converts BAM quality values (raw Phred) to SAM format (Phred+33).
// If all values are 0xFF, returns "*" (quality not stored).
func decodeQual(data []byte) string {
	allMissing := true
	for _, v := range data {
		if v != 0xFF {
			allMissing = false
			break
		}
	}
	if allMissing || len(data) == 0 {
		return "*"
	}
	buf := make([]byte, len(data))
	for i, v := range data {
		buf[i] = v + 33
	}
	return string(buf)
}

// decodeTags parses the auxiliary data section of a BAM record into a tag map.
func decodeTags(data []byte) map[string]SamTag {
	tags := make(map[string]SamTag)
	pos := 0
	for pos+3 <= len(data) {
		tag := string(data[pos : pos+2])
		valType := data[pos+2]
		pos += 3

		var samTag SamTag
		switch valType {
		case 'A': // printable character
			if pos >= len(data) {
				return tags
			}
			samTag = SamTag{Type: 'A', Value: string(data[pos : pos+1])}
			pos++

		case 'c': // int8
			if pos >= len(data) {
				return tags
			}
			samTag = SamTag{Type: 'i', Value: strconv.Itoa(int(int8(data[pos])))}
			pos++
		case 'C': // uint8
			if pos >= len(data) {
				return tags
			}
			samTag = SamTag{Type: 'i', Value: strconv.Itoa(int(data[pos]))}
			pos++
		case 's': // int16
			if pos+2 > len(data) {
				return tags
			}
			samTag = SamTag{Type: 'i', Value: strconv.Itoa(int(int16(binary.LittleEndian.Uint16(data[pos:]))))}
			pos += 2
		case 'S': // uint16
			if pos+2 > len(data) {
				return tags
			}
			samTag = SamTag{Type: 'i', Value: strconv.Itoa(int(binary.LittleEndian.Uint16(data[pos:])))}
			pos += 2
		case 'i': // int32
			if pos+4 > len(data) {
				return tags
			}
			samTag = SamTag{Type: 'i', Value: strconv.Itoa(int(int32(binary.LittleEndian.Uint32(data[pos:]))))}
			pos += 4
		case 'I': // uint32
			if pos+4 > len(data) {
				return tags
			}
			samTag = SamTag{Type: 'i', Value: strconv.FormatUint(uint64(binary.LittleEndian.Uint32(data[pos:])), 10)}
			pos += 4

		case 'f': // float32
			if pos+4 > len(data) {
				return tags
			}
			bits := binary.LittleEndian.Uint32(data[pos:])
			samTag = SamTag{Type: 'f', Value: strconv.FormatFloat(float64(math.Float32frombits(bits)), 'g', -1, 32)}
			pos += 4

		case 'Z': // NUL-terminated string
			end := pos
			for end < len(data) && data[end] != 0 {
				end++
			}
			samTag = SamTag{Type: 'Z', Value: string(data[pos:end])}
			pos = end + 1 // skip NUL

		case 'H': // hex string, NUL-terminated
			end := pos
			for end < len(data) && data[end] != 0 {
				end++
			}
			samTag = SamTag{Type: 'H', Value: string(data[pos:end])}
			pos = end + 1

		case 'B': // array
			if pos >= len(data) {
				return tags
			}
			samTag, pos = decodeArrayTag(data, pos)

		default:
			// Unknown type — skip rest of tags.
			return tags
		}

		tags[tag] = samTag
	}
	return tags
}

// decodeArrayTag decodes a B-type (array) auxiliary tag.
// Returns the decoded tag and the new position in data.
func decodeArrayTag(data []byte, pos int) (SamTag, int) {
	elemType := data[pos]
	pos++
	if pos+4 > len(data) {
		return SamTag{Type: 'B'}, pos
	}
	count := int(binary.LittleEndian.Uint32(data[pos:]))
	pos += 4

	var sb strings.Builder
	sb.WriteByte(elemType)

	for i := 0; i < count; i++ {
		sb.WriteByte(',')
		switch elemType {
		case 'c':
			if pos >= len(data) {
				return SamTag{Type: 'B', Value: sb.String()}, pos
			}
			sb.WriteString(strconv.Itoa(int(int8(data[pos]))))
			pos++
		case 'C':
			if pos >= len(data) {
				return SamTag{Type: 'B', Value: sb.String()}, pos
			}
			sb.WriteString(strconv.Itoa(int(data[pos])))
			pos++
		case 's':
			if pos+2 > len(data) {
				return SamTag{Type: 'B', Value: sb.String()}, pos
			}
			sb.WriteString(strconv.Itoa(int(int16(binary.LittleEndian.Uint16(data[pos:])))))
			pos += 2
		case 'S':
			if pos+2 > len(data) {
				return SamTag{Type: 'B', Value: sb.String()}, pos
			}
			sb.WriteString(strconv.Itoa(int(binary.LittleEndian.Uint16(data[pos:]))))
			pos += 2
		case 'i':
			if pos+4 > len(data) {
				return SamTag{Type: 'B', Value: sb.String()}, pos
			}
			sb.WriteString(strconv.Itoa(int(int32(binary.LittleEndian.Uint32(data[pos:])))))
			pos += 4
		case 'I':
			if pos+4 > len(data) {
				return SamTag{Type: 'B', Value: sb.String()}, pos
			}
			sb.WriteString(strconv.FormatUint(uint64(binary.LittleEndian.Uint32(data[pos:])), 10))
			pos += 4
		case 'f':
			if pos+4 > len(data) {
				return SamTag{Type: 'B', Value: sb.String()}, pos
			}
			bits := binary.LittleEndian.Uint32(data[pos:])
			sb.WriteString(strconv.FormatFloat(float64(math.Float32frombits(bits)), 'g', -1, 32))
			pos += 4
		}
	}

	return SamTag{Type: 'B', Value: sb.String()}, pos
}

