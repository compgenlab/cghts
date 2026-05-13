package bam

import (
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"os"
	"strconv"
	"strings"
	"sync"

	"github.com/compgen-io/cgkit/htsio"
	"github.com/compgen-io/cgkit/htsio/bgzf"
	"github.com/compgen-io/cgkit/htsio/tabix"
)

// Writer writes BAM files natively using the bgzf package.
// It implements the htsio.SamWriter interface and is safe for concurrent use.
type Writer struct {
	w       *bgzf.Writer
	f       *os.File   // non-nil if we opened the file
	header  *htsio.SamHeader
	refs    []bamRefInfo
	refIdx  map[string]int32 // ref name → index
	started bool
	closed  bool
	writeCh chan *htsio.SamRecord
	writeWg sync.WaitGroup
	err     error
	mu      sync.Mutex
}

// NewWriter creates a native BAM writer for the given output file.
// If filename is "-", writes to stdout.
func NewWriter(filename string, header *htsio.SamHeader) (*Writer, error) {
	return NewWriterWithThreads(filename, header, 1)
}

// NewWriterWithThreads creates a native BAM writer with parallel BGZF
// compression using the given number of threads. If threads <= 1, compression
// is single-threaded.
func NewWriterWithThreads(filename string, header *htsio.SamHeader, threads int) (*Writer, error) {
	if filename == "-" {
		return NewWriterFromWriterWithThreads(os.Stdout, header, threads), nil
	}
	f, err := os.Create(filename)
	if err != nil {
		return nil, err
	}
	var bw *Writer
	if threads > 1 {
		bw = newWriter(bgzf.NewParallelWriter(f, threads), header)
	} else {
		bw = newWriter(bgzf.NewWriter(f), header)
	}
	bw.f = f
	return bw, nil
}

// NewWriterFromWriter creates a native BAM writer that writes to w.
func NewWriterFromWriter(w io.Writer, header *htsio.SamHeader) *Writer {
	return newWriter(bgzf.NewWriter(w), header)
}

// NewWriterFromWriterWithThreads creates a native BAM writer that writes to w
// with parallel BGZF compression.
func NewWriterFromWriterWithThreads(w io.Writer, header *htsio.SamHeader, threads int) *Writer {
	if threads > 1 {
		return newWriter(bgzf.NewParallelWriter(w, threads), header)
	}
	return newWriter(bgzf.NewWriter(w), header)
}

func newWriter(w *bgzf.Writer, header *htsio.SamHeader) *Writer {
	bw := &Writer{
		w:      w,
		header: header,
	}
	if header != nil {
		hrefs := header.References()
		bw.refs = make([]bamRefInfo, len(hrefs))
		bw.refIdx = make(map[string]int32, len(hrefs))
		for i, hr := range hrefs {
			bw.refs[i] = bamRefInfo{name: hr.Name, length: int32(hr.Length)}
			bw.refIdx[hr.Name] = int32(i)
		}
	} else {
		bw.refIdx = make(map[string]int32)
	}
	return bw
}

// start writes the BAM header and starts the async writer goroutine.
func (bw *Writer) start() error {
	bw.mu.Lock()
	defer bw.mu.Unlock()

	if bw.started {
		return nil
	}

	if err := bw.writeHeader(); err != nil {
		return err
	}

	bw.writeCh = make(chan *htsio.SamRecord, 1024)
	bw.writeWg.Add(1)
	go func() {
		defer bw.writeWg.Done()
		for rec := range bw.writeCh {
			if err := bw.encodeRecord(rec); err != nil {
				bw.err = fmt.Errorf("bam write: %w", err)
				for range bw.writeCh {
				}
				return
			}
		}
	}()

	bw.started = true
	return nil
}

// writeHeader writes the BAM magic, header text, and reference dictionary.
func (bw *Writer) writeHeader() error {
	// Magic
	if _, err := bw.w.Write([]byte("BAM\x01")); err != nil {
		return err
	}

	// Header text
	var headerText []byte
	if bw.header != nil {
		headerText = []byte(bw.header.Text())
	}
	if err := binary.Write(bw.w, binary.LittleEndian, int32(len(headerText))); err != nil {
		return err
	}
	if _, err := bw.w.Write(headerText); err != nil {
		return err
	}

	// Reference sequences
	if err := binary.Write(bw.w, binary.LittleEndian, int32(len(bw.refs))); err != nil {
		return err
	}
	for _, ref := range bw.refs {
		name := []byte(ref.name + "\x00")
		if err := binary.Write(bw.w, binary.LittleEndian, int32(len(name))); err != nil {
			return err
		}
		if _, err := bw.w.Write(name); err != nil {
			return err
		}
		if err := binary.Write(bw.w, binary.LittleEndian, ref.length); err != nil {
			return err
		}
	}

	return nil
}

// Write sends a SamRecord to the async writer goroutine.
// It is safe for concurrent use by multiple goroutines.
func (bw *Writer) Write(rec *htsio.SamRecord) error {
	if err := bw.start(); err != nil {
		return err
	}
	bw.mu.Lock()
	if bw.err != nil {
		bw.mu.Unlock()
		return bw.err
	}
	bw.mu.Unlock()
	bw.writeCh <- rec
	return nil
}

// Close drains the write buffer, flushes the BGZF stream, and closes the file.
func (bw *Writer) Close() error {
	if bw.closed {
		return nil
	}
	bw.closed = true

	if !bw.started {
		// Ensure header is written even if no records were added.
		if err := bw.start(); err != nil {
			return err
		}
	}

	close(bw.writeCh)
	bw.writeWg.Wait()

	if err := bw.w.Close(); err != nil && bw.err == nil {
		bw.err = err
	}
	if bw.f != nil {
		if err := bw.f.Close(); err != nil && bw.err == nil {
			bw.err = err
		}
	}
	return bw.err
}

// encodeRecord encodes a SamRecord as a BAM binary record and writes it.
func (bw *Writer) encodeRecord(rec *htsio.SamRecord) error {
	// Resolve reference IDs.
	refID := int32(-1)
	if rec.RefName != "*" {
		if idx, ok := bw.refIdx[rec.RefName]; ok {
			refID = idx
		}
	}

	nextRefID := int32(-1)
	if rec.RefNext == "=" {
		nextRefID = refID
	} else if rec.RefNext != "*" {
		if idx, ok := bw.refIdx[rec.RefNext]; ok {
			nextRefID = idx
		}
	}

	// Parse CIGAR
	cigarOps := encodeCigar(rec.Cigar)

	// Encode sequence
	seqBytes := encodeSeqBytes(rec.Seq)
	seqLen := len(rec.Seq)
	if rec.Seq == "*" {
		seqLen = 0
	}

	// Encode quality
	qualBytes := encodeQualBytes(rec.Qual, seqLen)

	// Encode aux tags
	auxBytes := encodeAuxTags(rec.Tags, rec.TagOrder)

	// Read name (NUL-terminated)
	nameBytes := append([]byte(rec.ReadName), 0)

	// Compute BAM bin
	pos := int32(rec.Pos - 1) // SAM 1-based → BAM 0-based
	refLen := htsio.CigarRefLen(rec.Cigar)
	bin := tabix.Reg2Bin(int(pos), int(pos)+refLen)

	// Block size = fixed (32) + name + cigar + seq + qual + aux
	blockSize := int32(32 + len(nameBytes) + len(cigarOps)*4 + len(seqBytes) + seqLen + len(auxBytes))

	// Write block size
	if err := binary.Write(bw.w, binary.LittleEndian, blockSize); err != nil {
		return err
	}

	// Fixed fields
	if err := binary.Write(bw.w, binary.LittleEndian, refID); err != nil {
		return err
	}
	if err := binary.Write(bw.w, binary.LittleEndian, pos); err != nil {
		return err
	}
	bw.w.Write([]byte{byte(len(nameBytes))}) // l_read_name
	bw.w.Write([]byte{byte(rec.MapQ)})       // MAPQ
	if err := binary.Write(bw.w, binary.LittleEndian, uint16(bin)); err != nil {
		return err
	}
	if err := binary.Write(bw.w, binary.LittleEndian, uint16(len(cigarOps))); err != nil {
		return err
	}
	if err := binary.Write(bw.w, binary.LittleEndian, uint16(rec.Flag)); err != nil {
		return err
	}
	if err := binary.Write(bw.w, binary.LittleEndian, int32(seqLen)); err != nil {
		return err
	}
	nextPos := int32(rec.PosNext - 1)
	if err := binary.Write(bw.w, binary.LittleEndian, nextRefID); err != nil {
		return err
	}
	if err := binary.Write(bw.w, binary.LittleEndian, nextPos); err != nil {
		return err
	}
	if err := binary.Write(bw.w, binary.LittleEndian, int32(rec.InsertLen)); err != nil {
		return err
	}

	// Variable-length data
	bw.w.Write(nameBytes)
	for _, op := range cigarOps {
		binary.Write(bw.w, binary.LittleEndian, op)
	}
	bw.w.Write(seqBytes)
	bw.w.Write(qualBytes)
	bw.w.Write(auxBytes)

	return nil
}

// encodeCigar parses a SAM CIGAR string into packed BAM CIGAR ops.
func encodeCigar(cigar string) []uint32 {
	if cigar == "*" || cigar == "" {
		return nil
	}
	var ops []uint32
	num := 0
	for i := 0; i < len(cigar); i++ {
		c := cigar[i]
		if c >= '0' && c <= '9' {
			num = num*10 + int(c-'0')
		} else {
			var code uint32
			switch c {
			case 'M':
				code = 0
			case 'I':
				code = 1
			case 'D':
				code = 2
			case 'N':
				code = 3
			case 'S':
				code = 4
			case 'H':
				code = 5
			case 'P':
				code = 6
			case '=':
				code = 7
			case 'X':
				code = 8
			}
			ops = append(ops, uint32(num)<<4|code)
			num = 0
		}
	}
	return ops
}

// encodeSeqBytes converts an ASCII sequence to 4-bit packed BAM encoding.
func encodeSeqBytes(seq string) []byte {
	if seq == "*" || seq == "" {
		return nil
	}
	n := (len(seq) + 1) / 2
	out := make([]byte, n)
	for i := 0; i < len(seq); i++ {
		code := seqEncodeBase(seq[i])
		if i%2 == 0 {
			out[i/2] = code << 4
		} else {
			out[i/2] |= code
		}
	}
	return out
}

func seqEncodeBase(b byte) byte {
	switch b {
	case '=':
		return 0
	case 'A', 'a':
		return 1
	case 'C', 'c':
		return 2
	case 'M', 'm':
		return 3
	case 'G', 'g':
		return 4
	case 'R', 'r':
		return 5
	case 'S', 's':
		return 6
	case 'V', 'v':
		return 7
	case 'T', 't':
		return 8
	case 'W', 'w':
		return 9
	case 'Y', 'y':
		return 10
	case 'H', 'h':
		return 11
	case 'K', 'k':
		return 12
	case 'D', 'd':
		return 13
	case 'B', 'b':
		return 14
	case 'N', 'n':
		return 15
	default:
		return 15
	}
}

// encodeQualBytes converts a SAM quality string (Phred+33) to BAM raw Phred.
func encodeQualBytes(qual string, seqLen int) []byte {
	if qual == "*" || qual == "" {
		out := make([]byte, seqLen)
		for i := range out {
			out[i] = 0xFF
		}
		return out
	}
	out := make([]byte, len(qual))
	for i := 0; i < len(qual); i++ {
		out[i] = qual[i] - 33
	}
	return out
}

// encodeAuxTags encodes SAM optional tags into BAM binary format.
func encodeAuxTags(tags map[string]htsio.SamTag, tagOrder []string) []byte {
	if len(tags) == 0 {
		return nil
	}
	var buf []byte

	if len(tagOrder) > 0 {
		for _, tag := range tagOrder {
			st, ok := tags[tag]
			if !ok {
				continue
			}
			buf = encodeOneTag(buf, tag, st)
		}
	} else {
		for tag, st := range tags {
			buf = encodeOneTag(buf, tag, st)
		}
	}
	return buf
}

func encodeOneTag(buf []byte, tag string, st htsio.SamTag) []byte {
	buf = append(buf, tag[0], tag[1])
	switch st.Type {
	case 'A':
		buf = append(buf, 'A')
		if len(st.Value) > 0 {
			buf = append(buf, st.Value[0])
		} else {
			buf = append(buf, 0)
		}
	case 'i':
		v, _ := strconv.ParseInt(st.Value, 10, 64)
		if v >= 0 && v <= 255 {
			buf = append(buf, 'C', byte(v))
		} else if v >= -128 && v <= 127 {
			buf = append(buf, 'c', byte(int8(v)))
		} else if v >= 0 && v <= 65535 {
			buf = append(buf, 'S')
			buf = binary.LittleEndian.AppendUint16(buf, uint16(v))
		} else if v >= -32768 && v <= 32767 {
			buf = append(buf, 's')
			buf = binary.LittleEndian.AppendUint16(buf, uint16(int16(v)))
		} else if v >= 0 && v <= 4294967295 {
			buf = append(buf, 'I')
			buf = binary.LittleEndian.AppendUint32(buf, uint32(v))
		} else {
			buf = append(buf, 'i')
			buf = binary.LittleEndian.AppendUint32(buf, uint32(int32(v)))
		}
	case 'f':
		buf = append(buf, 'f')
		v, _ := strconv.ParseFloat(st.Value, 32)
		buf = binary.LittleEndian.AppendUint32(buf, math.Float32bits(float32(v)))
	case 'Z':
		buf = append(buf, 'Z')
		buf = append(buf, st.Value...)
		buf = append(buf, 0)
	case 'H':
		buf = append(buf, 'H')
		buf = append(buf, st.Value...)
		buf = append(buf, 0)
	case 'B':
		buf = append(buf, 'B')
		buf = encodeArrayTagValue(buf, st.Value)
	}
	return buf
}

// encodeArrayTagValue encodes a B-type array tag value (e.g. "C,1,2,3") into
// BAM binary and appends it to buf.
func encodeArrayTagValue(buf []byte, value string) []byte {
	parts := strings.Split(value, ",")
	if len(parts) < 2 {
		return buf
	}
	elemType := parts[0][0]
	buf = append(buf, elemType)
	count := len(parts) - 1
	buf = binary.LittleEndian.AppendUint32(buf, uint32(count))

	for _, p := range parts[1:] {
		switch elemType {
		case 'c':
			v, _ := strconv.ParseInt(p, 10, 8)
			buf = append(buf, byte(int8(v)))
		case 'C':
			v, _ := strconv.ParseUint(p, 10, 8)
			buf = append(buf, byte(v))
		case 's':
			v, _ := strconv.ParseInt(p, 10, 16)
			buf = binary.LittleEndian.AppendUint16(buf, uint16(int16(v)))
		case 'S':
			v, _ := strconv.ParseUint(p, 10, 16)
			buf = binary.LittleEndian.AppendUint16(buf, uint16(v))
		case 'i':
			v, _ := strconv.ParseInt(p, 10, 32)
			buf = binary.LittleEndian.AppendUint32(buf, uint32(int32(v)))
		case 'I':
			v, _ := strconv.ParseUint(p, 10, 32)
			buf = binary.LittleEndian.AppendUint32(buf, uint32(v))
		case 'f':
			v, _ := strconv.ParseFloat(p, 32)
			buf = binary.LittleEndian.AppendUint32(buf, math.Float32bits(float32(v)))
		}
	}
	return buf
}

