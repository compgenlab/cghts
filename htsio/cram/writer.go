package cram

import (
	"bytes"
	"compress/gzip"
	"crypto/md5"
	"encoding/binary"
	"fmt"
	"hash/crc32"
	"io"
	"math"
	"os"
	"strings"

	"github.com/compgen-io/cgltk/htsio"
)

// Version identifies a CRAM version for the writer.
type Version struct {
	major, minor byte
}

var (
	V2  = Version{2, 1}
	V3  = Version{3, 0}
	V31 = Version{3, 1}
)

// WriterOpts configures CRAM writing.
type WriterOpts struct {
	version         Version
	refPath         string
	level           int // compression level 0-9
	recordsPerSlice int
}

// NewWriterOpts returns default writer options (v3.1, level 6, 10000 records/slice).
func NewWriterOpts() *WriterOpts {
	return &WriterOpts{
		version:         V31,
		level:           6,
		recordsPerSlice: 10000,
	}
}

func (o *WriterOpts) SetVersion(v Version) *WriterOpts         { o.version = v; return o }
func (o *WriterOpts) Reference(path string) *WriterOpts        { o.refPath = path; return o }
func (o *WriterOpts) Level(n int) *WriterOpts                  { o.level = n; return o }
func (o *WriterOpts) RecordsPerSlice(n int) *WriterOpts        { o.recordsPerSlice = n; return o }

// Writer writes CRAM files. Implements htsio.SamWriter.
type Writer struct {
	w              io.Writer
	closer         io.Closer
	opts           *WriterOpts
	header         *htsio.SamHeader
	refs           []refInfo
	refMap         map[string]int32
	readGroupMap   map[string]int32
	refProv        *referenceProvider
	recordBuf      []*htsio.SamRecord
	recordCounter  int64
	headerWritten  bool
}

// NewWriter creates a CRAM writer. If filename is "-", writes to stdout.
func NewWriter(filename string, header *htsio.SamHeader, opts ...*WriterOpts) (*Writer, error) {
	var o *WriterOpts
	if len(opts) > 0 && opts[0] != nil {
		o = opts[0]
	} else {
		o = NewWriterOpts()
	}

	var w io.Writer
	var closer io.Closer
	if filename == "-" {
		w = os.Stdout
	} else {
		f, err := os.Create(filename)
		if err != nil {
			return nil, err
		}
		w = f
		closer = f
	}

	return newWriter(w, closer, header, o)
}

// NewWriterFromWriter creates a CRAM writer from an io.Writer.
func NewWriterFromWriter(w io.Writer, header *htsio.SamHeader, opts ...*WriterOpts) (*Writer, error) {
	var o *WriterOpts
	if len(opts) > 0 && opts[0] != nil {
		o = opts[0]
	} else {
		o = NewWriterOpts()
	}
	return newWriter(w, nil, header, o)
}

func newWriter(w io.Writer, closer io.Closer, header *htsio.SamHeader, opts *WriterOpts) (*Writer, error) {
	cw := &Writer{
		w:            w,
		closer:       closer,
		opts:         opts,
		header:       header,
		refMap:       make(map[string]int32),
		readGroupMap: make(map[string]int32),
	}

	// Build ref map from header.
	headerRefs := header.References()
	cw.refs = make([]refInfo, len(headerRefs))
	for i, hr := range headerRefs {
		cw.refs[i] = refInfo{name: hr.Name, length: hr.Length}
		cw.refMap[hr.Name] = int32(i)
	}

	// Build read group map.
	for i, rg := range header.ReadGroups() {
		cw.readGroupMap[rg] = int32(i)
	}

	// Set up reference provider.
	if opts.refPath != "" {
		cw.refProv = newReferenceProvider(opts.refPath)
	}

	return cw, nil
}

func (cw *Writer) version() Version { return cw.opts.version }
func (cw *Writer) majorVersion() byte { return cw.opts.version.major }

// Write buffers a record and flushes when the buffer is full.
func (cw *Writer) Write(rec *htsio.SamRecord) error {
	if !cw.headerWritten {
		if err := cw.writeFileHeader(); err != nil {
			return err
		}
		cw.headerWritten = true
	}

	cw.recordBuf = append(cw.recordBuf, rec)
	if len(cw.recordBuf) >= cw.opts.recordsPerSlice {
		return cw.flush()
	}
	return nil
}

// Close flushes remaining records, writes EOF, and closes the output.
func (cw *Writer) Close() error {
	if !cw.headerWritten {
		if err := cw.writeFileHeader(); err != nil {
			return err
		}
		cw.headerWritten = true
	}

	if len(cw.recordBuf) > 0 {
		if err := cw.flush(); err != nil {
			return err
		}
	}

	if err := cw.writeEOF(); err != nil {
		return err
	}

	if cw.closer != nil {
		return cw.closer.Close()
	}
	return nil
}

// writeFileHeader writes the file definition and header container.
func (cw *Writer) writeFileHeader() error {
	// File definition: "CRAM" + version + 20-byte file ID.
	var fd [26]byte
	fd[0], fd[1], fd[2], fd[3] = 'C', 'R', 'A', 'M'
	fd[4] = cw.opts.version.major
	fd[5] = cw.opts.version.minor
	// File ID: zeros (or could be random/hash)
	if _, err := cw.w.Write(fd[:]); err != nil {
		return fmt.Errorf("cram: writing file definition: %w", err)
	}

	// Header container.
	headerText := cw.header.Text()
	// Build header block data: int32 length prefix + header text.
	var headerBlock bytes.Buffer
	hLen := int32(len(headerText))
	binary.Write(&headerBlock, binary.LittleEndian, hLen)
	headerBlock.WriteString(headerText)

	// Write the header block as a single raw block in a container.
	blockData := headerBlock.Bytes()
	blockBytes, err := cw.encodeBlock(blockContentFileHeader, 0, blockMethodRaw, blockData)
	if err != nil {
		return err
	}

	return cw.writeContainer(-1, 0, 0, 0, 0, [][]byte{blockBytes}, nil)
}

// flush encodes buffered records into a container and writes it.
func (cw *Writer) flush() error {
	if len(cw.recordBuf) == 0 {
		return nil
	}

	records := cw.recordBuf
	cw.recordBuf = nil

	// Group records by reference. For simplicity, use multi-ref slice if mixed.
	refID := int32(-2) // multi-ref
	allSameRef := true
	firstRefID := int32(-1)
	if len(records) > 0 && records[0].RefName != "*" {
		if id, ok := cw.refMap[records[0].RefName]; ok {
			firstRefID = id
		}
	}
	for _, rec := range records {
		rid := int32(-1)
		if rec.RefName != "*" {
			if id, ok := cw.refMap[rec.RefName]; ok {
				rid = id
			}
		}
		if rid != firstRefID {
			allSameRef = false
			break
		}
	}
	if allSameRef {
		refID = firstRefID
	}

	// Convert SamRecords to cramRecords.
	cramRecords := make([]cramRecord, len(records))
	var minPos, maxEndPos int32
	minPos = 0x7FFFFFFF
	for i, rec := range records {
		cr := cw.samToCram(rec, refID == -2)
		cramRecords[i] = cr
		if cr.alignPos > 0 && cr.alignPos < minPos {
			minPos = cr.alignPos
		}
		endPos := cr.alignPos + cr.readLen
		if endPos > maxEndPos {
			maxEndPos = endPos
		}
	}
	if minPos == 0x7FFFFFFF {
		minPos = 0
	}
	alignmentSpan := int32(0)
	if maxEndPos > minPos {
		alignmentSpan = maxEndPos - minPos
	}

	// Build compression header and encode records.
	sliceData, err := cw.encodeSlice(cramRecords, refID, minPos, alignmentSpan)
	if err != nil {
		return fmt.Errorf("cram: encoding slice: %w", err)
	}

	numBases := int64(0)
	for _, r := range cramRecords {
		numBases += int64(r.readLen)
	}

	// Compute landmarks: byte offset from container content start to each slice header.
	// sliceData = [compHdrBlock, sliceHeaderBlock, extBlock1, extBlock2, ...]
	// The landmark for the single slice is the offset past the compression header block.
	var landmarks []int32
	if len(sliceData) > 1 {
		landmarks = []int32{int32(len(sliceData[0]))} // offset past compression header block
	}

	// Write container.
	err = cw.writeContainer(refID, minPos, alignmentSpan, int32(len(cramRecords)), numBases, sliceData, landmarks)
	if err != nil {
		return fmt.Errorf("cram: writing container: %w", err)
	}

	cw.recordCounter += int64(len(cramRecords))
	return nil
}

// samToCram converts a SamRecord to a cramRecord.
func (cw *Writer) samToCram(rec *htsio.SamRecord, multiRef bool) cramRecord {
	cr := cramRecord{
		bamFlags: int32(rec.Flag),
		readLen:  int32(len(rec.Seq)),
		alignPos: int32(rec.Pos),
		readName: rec.ReadName,
		mapQ:     int32(rec.MapQ),
	}

	if rec.Seq == "*" {
		cr.readLen = 0
	}

	// Reference ID.
	cr.refID = -1
	if rec.RefName != "*" {
		if id, ok := cw.refMap[rec.RefName]; ok {
			cr.refID = id
		}
	}

	// Read group.
	cr.readGroup = -1
	if rg, ok := rec.Tags["RG"]; ok {
		if id, ok := cw.readGroupMap[rg.Value]; ok {
			cr.readGroup = id
		}
	}

	// Mate info — always detached for simplicity.
	isPaired := rec.Flag&0x1 != 0
	if isPaired {
		cr.cramFlags |= 0x2 // detached
		if rec.Flag&0x8 != 0 {
			cr.mateFlags |= 0x1 // mate unmapped
		}
		if rec.Flag&0x20 != 0 {
			cr.mateFlags |= 0x2 // mate reverse
		}
		cr.mateRefID = -1
		mateRefName := rec.RefNext
		if mateRefName == "=" {
			cr.mateRefID = cr.refID
		} else if mateRefName != "*" {
			if id, ok := cw.refMap[mateRefName]; ok {
				cr.mateRefID = id
			}
		}
		cr.matePos = int32(rec.PosNext)
		cr.templateSize = int32(rec.InsertLen)
	}

	// Quality scores.
	isUnmapped := rec.Flag&0x4 != 0
	if rec.Qual != "*" {
		cr.cramFlags |= 0x1 // quality scores stored
		cr.qualScores = make([]byte, len(rec.Qual))
		for i, q := range rec.Qual {
			cr.qualScores[i] = byte(q) - 33
		}
	}

	if isUnmapped {
		// Store bases directly.
		if rec.Seq != "*" {
			cr.bases = []byte(rec.Seq)
		}
	} else {
		// Build features from CIGAR and sequence.
		cr.features = cw.buildFeatures(rec, &cr)
	}

	// Build tag values.
	cr.tagValues = make(map[int32][]byte)
	for _, tagName := range rec.TagOrder {
		if tagName == "RG" {
			continue // handled separately
		}
		tag := rec.Tags[tagName]
		tk := tagKey{id: [2]byte{tagName[0], tagName[1]}}
		tk.typ = samTypeToBAMType(tag.Type, tag.Value)
		key := tagKeyToITF8(tk)
		cr.tagValues[key] = samTagValueToBytes(tag)
	}
	// Also handle tags not in TagOrder but in Tags map.
	for tagName, tag := range rec.Tags {
		if tagName == "RG" {
			continue
		}
		found := false
		for _, n := range rec.TagOrder {
			if n == tagName {
				found = true
				break
			}
		}
		if found {
			continue
		}
		tk := tagKey{id: [2]byte{tagName[0], tagName[1]}}
		tk.typ = samTypeToBAMType(tag.Type, tag.Value)
		key := tagKeyToITF8(tk)
		cr.tagValues[key] = samTagValueToBytes(tag)
	}

	return cr
}

// buildFeatures parses the CIGAR string and builds CRAM features.
func (cw *Writer) buildFeatures(rec *htsio.SamRecord, cr *cramRecord) []cramFeature {
	cigar := rec.Cigar
	if cigar == "*" || cigar == "" {
		return nil
	}

	seq := rec.Seq

	// Get reference sequence for substitution detection.
	var refSeq []byte
	if cr.refID >= 0 && cw.refProv != nil && int(cr.refID) < len(cw.refs) {
		refSeq, _ = cw.refProv.getSequence(cw.refs[cr.refID].name)
	}

	var features []cramFeature
	readPos := 0  // 0-based position in read
	refPos := int(cr.alignPos) - 1 // 0-based ref position

	// Parse CIGAR.
	num := 0
	for i := 0; i < len(cigar); i++ {
		c := cigar[i]
		if c >= '0' && c <= '9' {
			num = num*10 + int(c-'0')
			continue
		}
		switch c {
		case 'M', '=', 'X':
			// Check each position for mismatches against reference.
			for j := 0; j < num; j++ {
				if readPos < len(seq) && refSeq != nil && refPos >= 0 && refPos < len(refSeq) {
					readBase := seq[readPos]
					refBase := refSeq[refPos]
					if readBase != refBase {
						features = append(features, cramFeature{
							code:    'X',
							pos:     int32(readPos + 1), // 1-based
							subCode: cw.substitutionCode(refBase, readBase),
						})
					}
				} else if readPos < len(seq) && refSeq == nil {
					// No reference: use read base feature
					features = append(features, cramFeature{
						code: 'B',
						pos:  int32(readPos + 1),
						base: seq[readPos],
					})
				}
				readPos++
				refPos++
			}
		case 'I':
			if num == 1 && readPos < len(seq) {
				features = append(features, cramFeature{
					code: 'i',
					pos:  int32(readPos + 1),
					base: seq[readPos],
				})
			} else if readPos+num <= len(seq) {
				bases := make([]byte, num)
				copy(bases, seq[readPos:readPos+num])
				features = append(features, cramFeature{
					code:  'I',
					pos:   int32(readPos + 1),
					bases: bases,
				})
			}
			readPos += num
		case 'D':
			features = append(features, cramFeature{
				code:   'D',
				pos:    int32(readPos + 1),
				length: int32(num),
			})
			refPos += num
		case 'N':
			features = append(features, cramFeature{
				code:   'N',
				pos:    int32(readPos + 1),
				length: int32(num),
			})
			refPos += num
		case 'S':
			if readPos+num <= len(seq) {
				bases := make([]byte, num)
				copy(bases, seq[readPos:readPos+num])
				features = append(features, cramFeature{
					code:  'S',
					pos:   int32(readPos + 1),
					bases: bases,
				})
			}
			readPos += num
		case 'H':
			features = append(features, cramFeature{
				code:   'H',
				pos:    int32(readPos + 1),
				length: int32(num),
			})
		case 'P':
			features = append(features, cramFeature{
				code:   'P',
				pos:    int32(readPos + 1),
				length: int32(num),
			})
		}
		num = 0
	}

	return features
}

// substitutionCode computes the 2-bit substitution code for a base replacement.
// This must match the substitution matrix that will be built.
func (cw *Writer) substitutionCode(refBase, readBase byte) byte {
	others := otherBases(refBase)
	for i, b := range others {
		if b == readBase {
			return byte(i)
		}
	}
	return 0
}

// encodeSlice encodes a set of cramRecords into slice blocks.
// Returns the compression header block bytes followed by slice block bytes.
func (cw *Writer) encodeSlice(records []cramRecord, refID, startPos, alignmentSpan int32) ([][]byte, error) {
	// Build substitution matrix from the actual mismatches.
	subMatrix := cw.buildSubstitutionMatrix(records)

	// Build tag dictionary.
	tagDict, tagComboMap := cw.buildTagDictionary(records)

	// Assign tagDictIdx to each record.
	for i := range records {
		records[i].tagDictIdx = tagComboMap[tagComboKey(records[i])]
	}

	// Build compression header.
	compHdr := cw.buildCompressionHeader(subMatrix, tagDict, records, refID)

	// Encode compression header block.
	compHdrData := cw.serializeCompressionHeader(compHdr, subMatrix, tagDict)
	compHdrBlock, err := cw.encodeBlock(blockContentCompressionHeader, 0, blockMethodRaw, compHdrData)
	if err != nil {
		return nil, err
	}

	// Encode the actual records into core + external blocks.
	sliceBlocks, err := cw.encodeRecords(records, compHdr, refID, startPos, alignmentSpan)
	if err != nil {
		return nil, err
	}

	// Combine: compression header block + slice blocks.
	result := make([][]byte, 0, 1+len(sliceBlocks))
	result = append(result, compHdrBlock)
	result = append(result, sliceBlocks...)
	return result, nil
}

// buildSubstitutionMatrix builds the 5-byte substitution matrix from record features.
func (cw *Writer) buildSubstitutionMatrix(records []cramRecord) [5][4]byte {
	// Use the default ordering: for each ref base, others in ACGTN order.
	var matrix [5][4]byte
	bases := [5]byte{'A', 'C', 'G', 'T', 'N'}
	for i := 0; i < 5; i++ {
		others := otherBases(bases[i])
		for j := 0; j < 4; j++ {
			matrix[i][j] = others[j]
		}
	}
	return matrix
}

// buildTagDictionary builds the tag dictionary from records.
// Returns the list of tag combos and a map from combo key to index.
func (cw *Writer) buildTagDictionary(records []cramRecord) ([][]tagKey, map[string]int32) {
	comboMap := make(map[string]int32)
	var combos [][]tagKey

	for i := range records {
		key := tagComboKey(records[i])
		if _, ok := comboMap[key]; !ok {
			comboMap[key] = int32(len(combos))
			combo := tagComboFromRecord(records[i])
			combos = append(combos, combo)
		}
	}
	return combos, comboMap
}

// tagComboKey returns a string key for the tag combination of a record.
func tagComboKey(rec cramRecord) string {
	var sb strings.Builder
	for key := range rec.tagValues {
		fmt.Fprintf(&sb, "%06x,", key)
	}
	return sb.String()
}

// tagComboFromRecord returns the tag keys for a record.
func tagComboFromRecord(rec cramRecord) []tagKey {
	var keys []tagKey
	for key := range rec.tagValues {
		tk := tagKey{
			id:  [2]byte{byte(key >> 16), byte(key >> 8)},
			typ: byte(key),
		}
		keys = append(keys, tk)
	}
	return keys
}

// Data series block IDs.
const (
	blockIDBF = 1
	blockIDCF = 2
	blockIDRI = 3
	blockIDRL = 4
	blockIDAP = 5
	blockIDRG = 6
	blockIDRN = 7
	blockIDMF = 8
	blockIDNS = 9
	blockIDNP = 10
	blockIDTS = 11
	blockIDNF = 12
	blockIDBA = 13
	blockIDQS = 14
	blockIDBS = 15
	blockIDIN = 16
	blockIDSC = 17
	blockIDDL = 18
	blockIDRS = 19
	blockIDHC = 20
	blockIDPD = 21
	blockIDMQ = 22
	blockIDTL = 23
	blockIDFN = 24
	blockIDFC = 25
	blockIDFP = 26
	blockIDBB = 27
	blockIDQQ = 28
	blockIDTagBase = 100 // tags start at 100+
)

type writerCompressionHeader struct {
	readNamesPreserved bool
	apDelta            bool
	refRequired        bool
	dataSeriesEncodings map[string]encodingDescriptor
	tagEncodings       map[int32]encodingDescriptor
}

type encodingDescriptor struct {
	codecID  int32
	params   []byte
}

func (cw *Writer) buildCompressionHeader(subMatrix [5][4]byte, tagDict [][]tagKey, records []cramRecord, refID int32) *writerCompressionHeader {
	ch := &writerCompressionHeader{
		readNamesPreserved: true,
		apDelta:            true,
		refRequired:        cw.refProv != nil,
		dataSeriesEncodings: make(map[string]encodingDescriptor),
		tagEncodings:       make(map[int32]encodingDescriptor),
	}

	// All data series use external encoding.
	ext := func(blockID int32) encodingDescriptor {
		var buf bytes.Buffer
		writeITF8(&buf, blockID)
		return encodingDescriptor{codecID: codecExternal, params: buf.Bytes()}
	}

	ch.dataSeriesEncodings["BF"] = ext(blockIDBF)
	ch.dataSeriesEncodings["CF"] = ext(blockIDCF)
	ch.dataSeriesEncodings["RL"] = ext(blockIDRL)
	ch.dataSeriesEncodings["AP"] = ext(blockIDAP)
	ch.dataSeriesEncodings["RG"] = ext(blockIDRG)
	ch.dataSeriesEncodings["MF"] = ext(blockIDMF)
	ch.dataSeriesEncodings["NS"] = ext(blockIDNS)
	ch.dataSeriesEncodings["NP"] = ext(blockIDNP)
	ch.dataSeriesEncodings["TS"] = ext(blockIDTS)
	ch.dataSeriesEncodings["NF"] = ext(blockIDNF)
	ch.dataSeriesEncodings["BA"] = ext(blockIDBA)
	ch.dataSeriesEncodings["QS"] = ext(blockIDQS)
	ch.dataSeriesEncodings["BS"] = ext(blockIDBS)
	ch.dataSeriesEncodings["DL"] = ext(blockIDDL)
	ch.dataSeriesEncodings["RS"] = ext(blockIDRS)
	ch.dataSeriesEncodings["HC"] = ext(blockIDHC)
	ch.dataSeriesEncodings["PD"] = ext(blockIDPD)
	ch.dataSeriesEncodings["MQ"] = ext(blockIDMQ)
	ch.dataSeriesEncodings["TL"] = ext(blockIDTL)
	ch.dataSeriesEncodings["FN"] = ext(blockIDFN)
	ch.dataSeriesEncodings["FC"] = ext(blockIDFC)
	ch.dataSeriesEncodings["FP"] = ext(blockIDFP)
	ch.dataSeriesEncodings["BB"] = ext(blockIDBB)
	ch.dataSeriesEncodings["QQ"] = ext(blockIDQQ)

	if refID == -2 {
		ch.dataSeriesEncodings["RI"] = ext(blockIDRI)
	}

	// RN: read name — byte array stop (NUL-terminated).
	{
		var buf bytes.Buffer
		buf.WriteByte(0) // stop byte = NUL
		writeITF8(&buf, blockIDRN)
		ch.dataSeriesEncodings["RN"] = encodingDescriptor{codecID: codecByteArrayStop, params: buf.Bytes()}
	}

	// IN: insertion bases — byte array stop (NUL-terminated).
	{
		var buf bytes.Buffer
		buf.WriteByte(0)
		writeITF8(&buf, blockIDIN)
		ch.dataSeriesEncodings["IN"] = encodingDescriptor{codecID: codecByteArrayStop, params: buf.Bytes()}
	}

	// SC: soft clip — byte array stop (NUL-terminated).
	{
		var buf bytes.Buffer
		buf.WriteByte(0)
		writeITF8(&buf, blockIDSC)
		ch.dataSeriesEncodings["SC"] = encodingDescriptor{codecID: codecByteArrayStop, params: buf.Bytes()}
	}

	// Tag encodings: byte array len with external length codec + external value codec.
	tagBlockID := int32(blockIDTagBase)
	for _, combo := range tagDict {
		for _, tk := range combo {
			key := tagKeyToITF8(tk)
			if _, ok := ch.tagEncodings[key]; ok {
				continue
			}
			// BYTE_ARRAY_LEN: len codec = EXTERNAL(tagBlockID), val codec = EXTERNAL(tagBlockID)
			var buf bytes.Buffer
			// Length codec: external
			writeITF8(&buf, codecExternal)
			var lenParams bytes.Buffer
			writeITF8(&lenParams, tagBlockID)
			writeITF8(&buf, int32(len(lenParams.Bytes())))
			buf.Write(lenParams.Bytes())
			// Value codec: external
			writeITF8(&buf, codecExternal)
			var valParams bytes.Buffer
			writeITF8(&valParams, tagBlockID)
			writeITF8(&buf, int32(len(valParams.Bytes())))
			buf.Write(valParams.Bytes())

			ch.tagEncodings[key] = encodingDescriptor{codecID: codecByteArrayLen, params: buf.Bytes()}
			tagBlockID++
		}
	}

	return ch
}

// serializeCompressionHeader serializes the compression header to bytes.
func (cw *Writer) serializeCompressionHeader(ch *writerCompressionHeader, subMatrix [5][4]byte, tagDict [][]tagKey) []byte {
	var buf bytes.Buffer

	// Preservation map.
	var pmBuf bytes.Buffer
	numPM := int32(5) // RN, AP, RR, SM, TD
	writeITF8(&pmBuf, numPM)

	// RN
	pmBuf.Write([]byte("RN"))
	if ch.readNamesPreserved {
		pmBuf.WriteByte(1)
	} else {
		pmBuf.WriteByte(0)
	}

	// AP
	pmBuf.Write([]byte("AP"))
	if ch.apDelta {
		pmBuf.WriteByte(1)
	} else {
		pmBuf.WriteByte(0)
	}

	// RR
	pmBuf.Write([]byte("RR"))
	if ch.refRequired {
		pmBuf.WriteByte(1)
	} else {
		pmBuf.WriteByte(0)
	}

	// SM: substitution matrix (5 bytes).
	pmBuf.Write([]byte("SM"))
	smBytes := encodeSubstitutionMatrix(subMatrix)
	pmBuf.Write(smBytes[:])

	// TD: tag dictionary.
	pmBuf.Write([]byte("TD"))
	tdData := encodeTagDictionary(tagDict)
	writeITF8(&pmBuf, int32(len(tdData)))
	pmBuf.Write(tdData)

	writeITF8(&buf, int32(pmBuf.Len()))
	buf.Write(pmBuf.Bytes())

	// Data series encoding map.
	var dsBuf bytes.Buffer
	writeITF8(&dsBuf, int32(len(ch.dataSeriesEncodings)))
	for name, enc := range ch.dataSeriesEncodings {
		dsBuf.Write([]byte(name[:2]))
		writeITF8(&dsBuf, enc.codecID)
		writeITF8(&dsBuf, int32(len(enc.params)))
		dsBuf.Write(enc.params)
	}
	writeITF8(&buf, int32(dsBuf.Len()))
	buf.Write(dsBuf.Bytes())

	// Tag encoding map.
	var teBuf bytes.Buffer
	writeITF8(&teBuf, int32(len(ch.tagEncodings)))
	for key, enc := range ch.tagEncodings {
		writeITF8(&teBuf, key)
		writeITF8(&teBuf, enc.codecID)
		writeITF8(&teBuf, int32(len(enc.params)))
		teBuf.Write(enc.params)
	}
	writeITF8(&buf, int32(teBuf.Len()))
	buf.Write(teBuf.Bytes())

	return buf.Bytes()
}

// encodeSubstitutionMatrix encodes the substitution matrix into 5 bytes.
func encodeSubstitutionMatrix(matrix [5][4]byte) [5]byte {
	var sm [5]byte
	bases := [5]byte{'A', 'C', 'G', 'T', 'N'}
	for i := 0; i < 5; i++ {
		others := otherBases(bases[i])
		var b byte
		for j, other := range others {
			// Find the position of this other base in the matrix row.
			for k := 0; k < 4; k++ {
				if matrix[i][k] == other {
					b |= byte(k) << uint(6-j*2)
					break
				}
			}
		}
		sm[i] = b
	}
	return sm
}

// encodeTagDictionary serializes the tag dictionary.
func encodeTagDictionary(tagDict [][]tagKey) []byte {
	var buf bytes.Buffer
	for _, combo := range tagDict {
		for _, tk := range combo {
			buf.Write(tk.id[:])
			buf.WriteByte(tk.typ)
		}
		buf.WriteByte(0) // NUL terminator
	}
	return buf.Bytes()
}

// encodeRecords encodes cramRecords into slice header + core + external blocks.
func (cw *Writer) encodeRecords(records []cramRecord, ch *writerCompressionHeader, refID, startPos, alignmentSpan int32) ([][]byte, error) {
	// External block buffers.
	externals := make(map[int32]*bytes.Buffer)
	getExt := func(id int32) *bytes.Buffer {
		if b, ok := externals[id]; ok {
			return b
		}
		b := &bytes.Buffer{}
		externals[id] = b
		return b
	}

	// Tag block ID map: tag key → block ID.
	tagBlockIDs := make(map[int32]int32)
	nextTagBlockID := int32(blockIDTagBase)
	for key := range ch.tagEncodings {
		tagBlockIDs[key] = nextTagBlockID
		nextTagBlockID++
	}

	prevAlignPos := startPos

	for _, rec := range records {
		// BF
		writeITF8(getExt(blockIDBF), rec.bamFlags)
		// CF
		writeITF8(getExt(blockIDCF), rec.cramFlags)
		// RI (multi-ref only)
		if refID == -2 {
			writeITF8(getExt(blockIDRI), rec.refID)
		}
		// RL
		writeITF8(getExt(blockIDRL), rec.readLen)
		// AP (delta encoded)
		if ch.apDelta {
			writeITF8(getExt(blockIDAP), rec.alignPos-prevAlignPos)
			prevAlignPos = rec.alignPos
		} else {
			writeITF8(getExt(blockIDAP), rec.alignPos)
		}
		// RG
		writeITF8(getExt(blockIDRG), rec.readGroup)

		// RN (read name, NUL-terminated)
		if ch.readNamesPreserved {
			rnBuf := getExt(blockIDRN)
			rnBuf.WriteString(rec.readName)
			rnBuf.WriteByte(0) // NUL stop
		}

		// Mate info.
		if rec.cramFlags&0x2 != 0 {
			// Detached.
			writeITF8(getExt(blockIDMF), rec.mateFlags)
			if !ch.readNamesPreserved {
				rnBuf := getExt(blockIDRN)
				rnBuf.WriteString(rec.readName)
				rnBuf.WriteByte(0)
			}
			writeITF8(getExt(blockIDNS), rec.mateRefID)
			writeITF8(getExt(blockIDNP), rec.matePos)
			writeITF8(getExt(blockIDTS), rec.templateSize)
		} else if rec.cramFlags&0x4 != 0 {
			writeITF8(getExt(blockIDNF), rec.recordsToNextFrag)
		}

		// TL (tag dictionary index)
		writeITF8(getExt(blockIDTL), rec.tagDictIdx)

		// Tags.
		for key, val := range rec.tagValues {
			blockID := tagBlockIDs[key]
			tb := getExt(blockID)
			writeITF8(tb, int32(len(val)))
			tb.Write(val)
		}

		// Mapped vs unmapped.
		isUnmapped := rec.bamFlags&0x4 != 0
		if !isUnmapped {
			// FN (number of features)
			writeITF8(getExt(blockIDFN), int32(len(rec.features)))

			prevFeatPos := int32(0)
			for _, feat := range rec.features {
				// FC (feature code)
				getExt(blockIDFC).WriteByte(feat.code)
				// FP (delta position)
				writeITF8(getExt(blockIDFP), feat.pos-prevFeatPos)
				prevFeatPos = feat.pos

				switch feat.code {
				case 'X': // Substitution
					getExt(blockIDBS).WriteByte(feat.subCode)
				case 'B': // Base + quality
					getExt(blockIDBA).WriteByte(feat.base)
					getExt(blockIDQS).WriteByte(feat.qual)
				case 'I': // Multi-base insertion
					inBuf := getExt(blockIDIN)
					inBuf.Write(feat.bases)
					inBuf.WriteByte(0) // NUL stop
				case 'i': // Single insertion
					getExt(blockIDBA).WriteByte(feat.base)
				case 'D': // Deletion
					writeITF8(getExt(blockIDDL), feat.length)
				case 'S': // Soft clip
					scBuf := getExt(blockIDSC)
					scBuf.Write(feat.bases)
					scBuf.WriteByte(0) // NUL stop
				case 'N': // Ref skip
					writeITF8(getExt(blockIDRS), feat.length)
				case 'H': // Hard clip
					writeITF8(getExt(blockIDHC), feat.length)
				case 'P': // Padding
					writeITF8(getExt(blockIDPD), feat.length)
				case 'b': // Base stretch
					bbBuf := getExt(blockIDBB)
					writeITF8(bbBuf, int32(len(feat.bases)))
					bbBuf.Write(feat.bases)
				case 'q': // Quality stretch
					qqBuf := getExt(blockIDQQ)
					writeITF8(qqBuf, int32(len(feat.quals)))
					qqBuf.Write(feat.quals)
				case 'Q': // Single quality
					getExt(blockIDQS).WriteByte(feat.qual)
				}
			}

			// MQ (mapping quality)
			writeITF8(getExt(blockIDMQ), rec.mapQ)

			// QS (quality scores)
			if rec.cramFlags&0x1 != 0 {
				qsBuf := getExt(blockIDQS)
				for _, q := range rec.qualScores {
					qsBuf.WriteByte(q)
				}
			}
		} else {
			// BA (unmapped bases)
			baBuf := getExt(blockIDBA)
			for _, b := range rec.bases {
				baBuf.WriteByte(b)
			}

			// QS (quality scores)
			if rec.cramFlags&0x1 != 0 {
				qsBuf := getExt(blockIDQS)
				for _, q := range rec.qualScores {
					qsBuf.WriteByte(q)
				}
			}
		}
	}

	// Compute reference MD5 for slice header.
	var refMD5 [16]byte
	if refID >= 0 && cw.refProv != nil && int(refID) < len(cw.refs) {
		if refSeq, err := cw.refProv.getSequence(cw.refs[refID].name); err == nil {
			// MD5 of the reference region covered by this slice.
			start := int(startPos) - 1
			end := start + int(alignmentSpan)
			if start < 0 {
				start = 0
			}
			if end > len(refSeq) {
				end = len(refSeq)
			}
			if start < end {
				refMD5 = md5.Sum(refSeq[start:end])
			}
		}
	}

	// Build external block content IDs.
	var contentIDs []int32
	for id := range externals {
		contentIDs = append(contentIDs, id)
	}

	// Build slice header.
	sh := &sliceHeader{
		refSeqID:        refID,
		alignmentStart:  startPos,
		alignmentSpan:   alignmentSpan,
		numRecords:      int32(len(records)),
		recordCounter:   cw.recordCounter,
		numBlocks:       int32(1 + len(contentIDs)), // core block + external blocks
		blockContentIDs: contentIDs,
		embeddedRefID:   -1,
		referenceMD5:    refMD5,
	}

	// Serialize slice header.
	shData := cw.serializeSliceHeader(sh)
	shBlock, err := cw.encodeBlock(blockContentSliceHeader, 0, blockMethodRaw, shData)
	if err != nil {
		return nil, err
	}

	// Build result blocks: slice header, core data block (empty), then external blocks.
	var blocks [][]byte
	blocks = append(blocks, shBlock)

	// Core data block (empty — we use only external encodings).
	coreBlock, err := cw.encodeBlock(blockContentCoreData, 0, blockMethodRaw, nil)
	if err != nil {
		return nil, err
	}
	blocks = append(blocks, coreBlock)

	for _, id := range contentIDs {
		data := externals[id].Bytes()
		blk, err := cw.compressAndEncodeBlock(blockContentExternalData, id, data)
		if err != nil {
			return nil, err
		}
		blocks = append(blocks, blk)
	}

	return blocks, nil
}

// serializeSliceHeader serializes a slice header to bytes.
func (cw *Writer) serializeSliceHeader(sh *sliceHeader) []byte {
	var buf bytes.Buffer
	writeITF8(&buf, sh.refSeqID)
	writeITF8(&buf, sh.alignmentStart)
	writeITF8(&buf, sh.alignmentSpan)
	writeITF8(&buf, sh.numRecords)
	writeLTF8(&buf, sh.recordCounter)
	writeITF8(&buf, sh.numBlocks)
	writeITF8Array(&buf, sh.blockContentIDs)
	writeITF8(&buf, sh.embeddedRefID)
	buf.Write(sh.referenceMD5[:])
	return buf.Bytes()
}

// compressAndEncodeBlock compresses data using multiple methods and picks the smallest.
func (cw *Writer) compressAndEncodeBlock(contentType byte, contentID int32, data []byte) ([]byte, error) {
	if len(data) == 0 {
		return cw.encodeBlock(contentType, contentID, blockMethodRaw, data)
	}

	// Try raw.
	best, err := cw.encodeBlock(contentType, contentID, blockMethodRaw, data)
	if err != nil {
		return nil, err
	}

	// Try gzip.
	if candidate, err := cw.encodeBlock(contentType, contentID, blockMethodGzip, data); err == nil {
		if len(candidate) < len(best) {
			best = candidate
		}
	}

	// Try rANS 4x8 (v3.0+).
	if cw.majorVersion() >= 3 {
		if candidate, err := cw.encodeBlock(contentType, contentID, blockMethodRans4x8, data); err == nil {
			if len(candidate) < len(best) {
				best = candidate
			}
		}
	}

	// Try rANS Nx16 (v3.1+).
	if cw.majorVersion() >= 3 && cw.opts.version.minor >= 1 {
		if candidate, err := cw.encodeBlock(contentType, contentID, blockMethodRans4x16, data); err == nil {
			if len(candidate) < len(best) {
				best = candidate
			}
		}
	}

	return best, nil
}

// encodeBlock creates the binary representation of a CRAM block.
func (cw *Writer) encodeBlock(contentType byte, contentID int32, method byte, data []byte) ([]byte, error) {
	var compData []byte
	switch method {
	case blockMethodGzip:
		var buf bytes.Buffer
		gz, err := gzip.NewWriterLevel(&buf, cw.opts.level)
		if err != nil {
			return nil, err
		}
		if _, err := gz.Write(data); err != nil {
			gz.Close()
			return nil, err
		}
		gz.Close()
		compData = buf.Bytes()
		if len(compData) >= len(data) {
			method = blockMethodRaw
			compData = data
		}
	case blockMethodRans4x8:
		// Try order-0 and order-1, pick the smaller.
		enc0 := encodeRans4x8(data, 0)
		enc1 := encodeRans4x8(data, 1)
		if len(enc1) < len(enc0) {
			compData = enc1
		} else {
			compData = enc0
		}
		if len(compData) >= len(data) {
			method = blockMethodRaw
			compData = data
		}
	case blockMethodRans4x16:
		compData = encodeRansNx16(data)
		if len(compData) >= len(data) {
			method = blockMethodRaw
			compData = data
		}
	default:
		compData = data
	}

	var buf bytes.Buffer

	// For v3+, we need CRC32 of the block header + data (not including the CRC itself).
	var h hash32
	var tw io.Writer
	if cw.majorVersion() >= 3 {
		h = crc32.NewIEEE()
		tw = io.MultiWriter(&buf, h)
	} else {
		tw = &buf
	}

	// Method byte.
	tw.Write([]byte{method})
	// Content type.
	tw.Write([]byte{contentType})
	// Content ID.
	writeITF8(tw, contentID)
	// Compressed size.
	writeITF8(tw, int32(len(compData)))
	// Raw size.
	writeITF8(tw, int32(len(data)))
	// Compressed data.
	tw.Write(compData)

	// CRC32 for v3+.
	if cw.majorVersion() >= 3 {
		crc := h.Sum32()
		var crcBuf [4]byte
		binary.LittleEndian.PutUint32(crcBuf[:], crc)
		buf.Write(crcBuf[:])
	}

	return buf.Bytes(), nil
}

// writeContainer writes a complete container (header + blocks).
// landmarks are byte offsets from the start of the container content to each slice header.
func (cw *Writer) writeContainer(refSeqID, startPos, alignmentSpan, numRecords int32, numBases int64, blocks [][]byte, landmarks []int32) error {
	// Compute total block content length.
	blockLen := int32(0)
	for _, b := range blocks {
		blockLen += int32(len(b))
	}

	// Build the complete header as a single byte slice.
	var hdrBuf bytes.Buffer

	// Length: int32 LE (= total size of all blocks).
	binary.Write(&hdrBuf, binary.LittleEndian, blockLen)
	// Ref seq ID.
	writeITF8(&hdrBuf, refSeqID)
	// Start pos.
	writeITF8(&hdrBuf, startPos)
	// Alignment span.
	writeITF8(&hdrBuf, alignmentSpan)
	// Num records.
	writeITF8(&hdrBuf, numRecords)
	// Record counter.
	if cw.majorVersion() >= 3 {
		writeLTF8(&hdrBuf, cw.recordCounter)
	} else {
		writeITF8(&hdrBuf, int32(cw.recordCounter))
	}
	// Bases.
	writeLTF8(&hdrBuf, numBases)
	// Num blocks.
	writeITF8(&hdrBuf, int32(len(blocks)))
	// Landmarks.
	writeITF8Array(&hdrBuf, landmarks)

	hdrBytes := hdrBuf.Bytes()

	// Write header.
	if _, err := cw.w.Write(hdrBytes); err != nil {
		return err
	}

	// CRC32 for v3+.
	if cw.majorVersion() >= 3 {
		crc := crc32.ChecksumIEEE(hdrBytes)
		var crcBuf [4]byte
		binary.LittleEndian.PutUint32(crcBuf[:], crc)
		if _, err := cw.w.Write(crcBuf[:]); err != nil {
			return err
		}
	}

	// Write blocks.
	for _, b := range blocks {
		if _, err := cw.w.Write(b); err != nil {
			return err
		}
	}

	return nil
}

// writeEOF writes the EOF container.
func (cw *Writer) writeEOF() error {
	if cw.majorVersion() >= 3 {
		// Use the exact htslib v3 EOF byte sequence.
		// This is the standard 38-byte EOF marker that samtools expects.
		eof := []byte{
			0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f,
			0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00, 0x00, 0x01,
			0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06,
			0x06, 0x01, 0x00, 0x01, 0x00, 0x01, 0x00, 0xee, 0x63,
			0x01, 0x4b,
		}
		_, err := cw.w.Write(eof)
		return err
	}

	// v2 EOF container: matches htslib format.
	// Block data is a minimal empty compression header: 3 empty maps.
	minCompHdr := []byte{0x01, 0x00, 0x01, 0x00, 0x01, 0x00}
	emptyBlock, err := cw.encodeBlock(blockContentCompressionHeader, 0, blockMethodRaw, minCompHdr)
	if err != nil {
		return err
	}
	return cw.writeContainer(-1, 0x454f46, 0, 0, 0, [][]byte{emptyBlock}, nil)
}

// samTypeToBAMType converts a SAM tag type character to the BAM/CRAM internal type.
func samTypeToBAMType(samType byte, value string) byte {
	switch samType {
	case 'i':
		// Determine the smallest integer type.
		v := 0
		fmt.Sscanf(value, "%d", &v)
		if v >= 0 && v <= 255 {
			return 'C'
		}
		if v >= -128 && v <= 127 {
			return 'c'
		}
		if v >= 0 && v <= 65535 {
			return 'S'
		}
		if v >= -32768 && v <= 32767 {
			return 's'
		}
		if v >= 0 {
			return 'I'
		}
		return 'i'
	default:
		return samType
	}
}

// samTagValueToBytes converts a SAM tag to its binary representation for CRAM.
func samTagValueToBytes(tag htsio.SamTag) []byte {
	switch tag.Type {
	case 'A':
		if len(tag.Value) > 0 {
			return []byte{tag.Value[0]}
		}
		return []byte{0}
	case 'i':
		v := 0
		fmt.Sscanf(tag.Value, "%d", &v)
		switch {
		case v >= 0 && v <= 255:
			return []byte{byte(v)}
		case v >= -128 && v <= 127:
			return []byte{byte(int8(v))}
		case v >= 0 && v <= 65535:
			var buf [2]byte
			binary.LittleEndian.PutUint16(buf[:], uint16(v))
			return buf[:]
		case v >= -32768 && v <= 32767:
			var buf [2]byte
			binary.LittleEndian.PutUint16(buf[:], uint16(int16(v)))
			return buf[:]
		default:
			var buf [4]byte
			binary.LittleEndian.PutUint32(buf[:], uint32(int32(v)))
			return buf[:]
		}
	case 'f':
		var f float64
		fmt.Sscanf(tag.Value, "%f", &f)
		var buf [4]byte
		binary.LittleEndian.PutUint32(buf[:], math.Float32bits(float32(f)))
		return buf[:]
	case 'Z', 'H':
		return []byte(tag.Value)
	case 'B':
		return []byte(tag.Value)
	default:
		return []byte(tag.Value)
	}
}
