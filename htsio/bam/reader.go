package bam

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"io"
	"iter"
	"math"
	"os"
	"strconv"
	"strings"

	"github.com/compgen-io/cgkit/htsio"
	"github.com/compgen-io/cgkit/htsio/bgzf"
	"github.com/compgen-io/cgkit/htsio/tabix"
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

// Reader reads BAM files natively using the bgzf package.
// It implements the htsio.SamReader interface.
type Reader struct {
	r        *bgzf.Reader
	src      io.ReadCloser // underlying file/reader, closed on Close()
	filename string        // original filename (for finding .bai)
	refs     []bamRefInfo  // reference sequences from the BAM header
	refMap   map[string]int // ref name → index
	hdr      *htsio.SamHeader
	opts     *htsio.SamReaderOpts
	idx      *tabix.BinIndex     // BAI index, loaded lazily on first Query()
	ir       *bgzf.IndexedReader // created lazily for Query()
	err      error               // sticky error
}

func init() {
	htsio.RegisterReader(htsio.ReaderRegistration{
		Detect: func(magic []byte) bool {
			return bytes.HasPrefix(magic, []byte("BAM\x01"))
		},
		NewFromFile: func(filename string, opts *htsio.SamReaderOpts) (htsio.SamReader, error) {
			f, err := os.Open(filename)
			if err != nil {
				return nil, err
			}
			return NewReader(f, filename, opts)
		},
		NewFromStream: func(r io.ReadCloser, opts *htsio.SamReaderOpts) (htsio.SamReader, error) {
			return NewReader(r, "", opts)
		},
	})
}

// NewReader creates a BAM reader from an io.ReadCloser.
// The reader must be positioned at the start of a BAM file (possibly with
// peeked bytes prepended via io.MultiReader).
func NewReader(rc io.ReadCloser, filename string, opts *htsio.SamReaderOpts) (*Reader, error) {
	if opts == nil {
		opts = htsio.NewSamReaderOpts()
	}

	br := &Reader{
		r:        bgzf.NewReader(rc),
		src:      rc,
		filename: filename,
		opts:     opts,
	}

	if err := br.readHeader(); err != nil {
		rc.Close()
		return nil, fmt.Errorf("bam: reading header: %w", err)
	}

	// Build ref name → index map.
	br.refMap = make(map[string]int, len(br.refs))
	for i, r := range br.refs {
		br.refMap[r.name] = i
	}

	return br, nil
}

// Header returns the parsed SAM header.
func (b *Reader) Header() (*htsio.SamHeader, error) {
	return b.hdr, nil
}

// Close releases resources.
func (b *Reader) Close() error {
	if b.src != nil {
		return b.src.Close()
	}
	return nil
}

// Query returns an iterator over records overlapping the 0-based half-open
// region [start, end) on the given reference. Requires a .bai index file
// at filename.bai.
func (b *Reader) Query(ref string, start, end int) (iter.Seq2[*htsio.SamRecord, error], error) {
	if b.filename == "" {
		return nil, fmt.Errorf("bam: Query requires a file-backed reader")
	}

	refID, ok := b.refMap[ref]
	if !ok {
		return nil, fmt.Errorf("bam: unknown reference %q", ref)
	}

	// Lazily load the BAI index.
	if b.idx == nil {
		baiPath := b.filename + ".bai"
		idx, err := tabix.LoadBAI(baiPath)
		if err != nil {
			return nil, fmt.Errorf("bam: loading BAI index: %w", err)
		}
		b.idx = idx
	}

	// Lazily create the indexed reader (shares the underlying file).
	if b.ir == nil {
		f, err := os.Open(b.filename)
		if err != nil {
			return nil, err
		}
		b.ir = bgzf.NewIndexedReader(f)
	}

	chunks := b.idx.Query(refID, start, end)
	if len(chunks) == 0 {
		return func(yield func(*htsio.SamRecord, error) bool) {}, nil
	}

	return b.iterChunks(chunks, refID, start, end), nil
}

// iterChunks returns an iterator that reads BAM records from the given
// chunks, filtering to those overlapping [start, end) on refID.
func (b *Reader) iterChunks(chunks []tabix.Chunk, refID, start, end int) iter.Seq2[*htsio.SamRecord, error] {
	return func(yield func(*htsio.SamRecord, error) bool) {
		for ci, chunk := range chunks {
			if err := b.ir.SeekToVirtualOffset(chunk.Begin); err != nil {
				yield(nil, fmt.Errorf("bam: seeking to chunk: %w", err))
				return
			}

			for {
				// Check if we've passed the chunk end.
				if ci < len(chunks) && b.ir.VirtualTell() >= chunk.End {
					break
				}

				rec, err := readBamRecord(b.ir, b.refs)
				if err != nil {
					if err == io.EOF {
						break
					}
					if !yield(nil, err) {
						return
					}
					return
				}

				// Check reference.
				recRefID := -1
				if rec.RefName != "*" {
					if id, ok := b.refMap[rec.RefName]; ok {
						recRefID = id
					}
				}
				if recRefID != refID {
					return // past our reference
				}

				// Record position (0-based).
				recStart := rec.Pos - 1
				recEnd := recStart + htsio.CigarRefLen(rec.Cigar)

				if recEnd <= start {
					continue // before query region
				}
				if recStart >= end {
					return // past query region
				}

				if !b.opts.PassesFilters(rec) {
					continue
				}

				if !yield(rec, nil) {
					return
				}
			}
		}
	}
}

// readHeader reads the BAM magic, header text, and reference sequence dictionary.
func (b *Reader) readHeader() error {
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

	b.hdr = htsio.NewSamHeader()
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

// Records returns an iterator over all records in the BAM file.
func (b *Reader) Records() iter.Seq2[*htsio.SamRecord, error] {
	return func(yield func(*htsio.SamRecord, error) bool) {
		if b.err != nil {
			yield(nil, b.err)
			return
		}
		for {
			rec, err := readBamRecord(b.r, b.refs)
			if err != nil {
				if err != io.EOF {
					yield(nil, err)
				}
				b.err = err
				return
			}
			if b.opts.PassesFilters(rec) {
				if !yield(rec, nil) {
					return
				}
			}
		}
	}
}

// readBamRecord reads a single BAM alignment record from r and converts it
// to a SamRecord, using refs to resolve reference IDs to names.
func readBamRecord(r io.Reader, refs []bamRefInfo) (*htsio.SamRecord, error) {
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
	tags, tagOrder := decodeTags(buf[offset:])

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

	rec := &htsio.SamRecord{
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
		TagOrder:  tagOrder,
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

// decodeTags parses the auxiliary data section of a BAM record into a tag map
// and an ordered list of tag keys.
func decodeTags(data []byte) (map[string]htsio.SamTag, []string) {
	tags := make(map[string]htsio.SamTag)
	var order []string
	pos := 0
	for pos+3 <= len(data) {
		tag := string(data[pos : pos+2])
		valType := data[pos+2]
		pos += 3

		var samTag htsio.SamTag
		switch valType {
		case 'A': // printable character
			if pos >= len(data) {
				return tags, order
			}
			samTag = htsio.SamTag{Type: 'A', Value: string(data[pos : pos+1])}
			pos++

		case 'c': // int8
			if pos >= len(data) {
				return tags, order
			}
			samTag = htsio.SamTag{Type: 'i', Value: strconv.Itoa(int(int8(data[pos])))}
			pos++
		case 'C': // uint8
			if pos >= len(data) {
				return tags, order
			}
			samTag = htsio.SamTag{Type: 'i', Value: strconv.Itoa(int(data[pos]))}
			pos++
		case 's': // int16
			if pos+2 > len(data) {
				return tags, order
			}
			samTag = htsio.SamTag{Type: 'i', Value: strconv.Itoa(int(int16(binary.LittleEndian.Uint16(data[pos:]))))}
			pos += 2
		case 'S': // uint16
			if pos+2 > len(data) {
				return tags, order
			}
			samTag = htsio.SamTag{Type: 'i', Value: strconv.Itoa(int(binary.LittleEndian.Uint16(data[pos:])))}
			pos += 2
		case 'i': // int32
			if pos+4 > len(data) {
				return tags, order
			}
			samTag = htsio.SamTag{Type: 'i', Value: strconv.Itoa(int(int32(binary.LittleEndian.Uint32(data[pos:]))))}
			pos += 4
		case 'I': // uint32
			if pos+4 > len(data) {
				return tags, order
			}
			samTag = htsio.SamTag{Type: 'i', Value: strconv.FormatUint(uint64(binary.LittleEndian.Uint32(data[pos:])), 10)}
			pos += 4

		case 'f': // float32
			if pos+4 > len(data) {
				return tags, order
			}
			bits := binary.LittleEndian.Uint32(data[pos:])
			samTag = htsio.SamTag{Type: 'f', Value: strconv.FormatFloat(float64(math.Float32frombits(bits)), 'g', -1, 32)}
			pos += 4

		case 'Z': // NUL-terminated string
			end := pos
			for end < len(data) && data[end] != 0 {
				end++
			}
			samTag = htsio.SamTag{Type: 'Z', Value: string(data[pos:end])}
			pos = end + 1 // skip NUL

		case 'H': // hex string, NUL-terminated
			end := pos
			for end < len(data) && data[end] != 0 {
				end++
			}
			samTag = htsio.SamTag{Type: 'H', Value: string(data[pos:end])}
			pos = end + 1

		case 'B': // array
			if pos >= len(data) {
				return tags, order
			}
			samTag, pos = decodeArrayTag(data, pos)

		default:
			// Unknown type — skip rest of tags.
			return tags, order
		}

		tags[tag] = samTag
		order = append(order, tag)
	}
	return tags, order
}

// decodeArrayTag decodes a B-type (array) auxiliary tag.
// Returns the decoded tag and the new position in data.
func decodeArrayTag(data []byte, pos int) (htsio.SamTag, int) {
	elemType := data[pos]
	pos++
	if pos+4 > len(data) {
		return htsio.SamTag{Type: 'B'}, pos
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
				return htsio.SamTag{Type: 'B', Value: sb.String()}, pos
			}
			sb.WriteString(strconv.Itoa(int(int8(data[pos]))))
			pos++
		case 'C':
			if pos >= len(data) {
				return htsio.SamTag{Type: 'B', Value: sb.String()}, pos
			}
			sb.WriteString(strconv.Itoa(int(data[pos])))
			pos++
		case 's':
			if pos+2 > len(data) {
				return htsio.SamTag{Type: 'B', Value: sb.String()}, pos
			}
			sb.WriteString(strconv.Itoa(int(int16(binary.LittleEndian.Uint16(data[pos:])))))
			pos += 2
		case 'S':
			if pos+2 > len(data) {
				return htsio.SamTag{Type: 'B', Value: sb.String()}, pos
			}
			sb.WriteString(strconv.Itoa(int(binary.LittleEndian.Uint16(data[pos:]))))
			pos += 2
		case 'i':
			if pos+4 > len(data) {
				return htsio.SamTag{Type: 'B', Value: sb.String()}, pos
			}
			sb.WriteString(strconv.Itoa(int(int32(binary.LittleEndian.Uint32(data[pos:])))))
			pos += 4
		case 'I':
			if pos+4 > len(data) {
				return htsio.SamTag{Type: 'B', Value: sb.String()}, pos
			}
			sb.WriteString(strconv.FormatUint(uint64(binary.LittleEndian.Uint32(data[pos:])), 10))
			pos += 4
		case 'f':
			if pos+4 > len(data) {
				return htsio.SamTag{Type: 'B', Value: sb.String()}, pos
			}
			bits := binary.LittleEndian.Uint32(data[pos:])
			sb.WriteString(strconv.FormatFloat(float64(math.Float32frombits(bits)), 'g', -1, 32))
			pos += 4
		}
	}

	return htsio.SamTag{Type: 'B', Value: sb.String()}, pos
}
