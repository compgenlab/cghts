package cram

import (
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
)

var debugReconstruct = false

// cramRecord holds a decoded CRAM record before conversion to SamRecord fields.
type cramRecord struct {
	bamFlags  int32
	cramFlags int32
	refID     int32
	readLen   int32
	alignPos  int32 // 1-based
	readGroup int32
	readName  string
	mapQ      int32

	// Mate info (detached)
	mateFlags int32
	mateRefID int32
	matePos   int32
	templateSize int32

	// Mate info (downstream)
	recordsToNextFrag int32

	// Features
	features []cramFeature

	// Tags
	tagDictIdx int32
	tagValues  map[int32][]byte // key = ITF8 tag key, value = raw tag bytes

	// Quality and sequence
	qualScores []byte
	bases      []byte // for unmapped reads
}

// cramFeature represents a single read feature (substitution, insertion, etc.)
type cramFeature struct {
	code byte
	pos  int32 // 1-based position within the read

	// Fields used depend on code:
	base    byte   // B, X, i
	qual    byte   // B, Q
	bases   []byte // I, S, b
	quals   []byte // q
	subCode byte   // X (substitution code)
	length  int32  // D, N, P, H, RS
}

// decodeSliceRecords decodes all records in a slice.
func decodeSliceRecords(
	sh *sliceHeader,
	ch *compressionHeader,
	coreData []byte,
	externalBlocks map[int32][]byte,
	refs []refInfo,
	refSeqs map[int32][]byte, // refID → reference sequence bytes
) ([]cramRecord, error) {
	core := newBitReader(coreData)

	// Track positions in external blocks.
	extPos := make(map[int32]*int)
	for id := range externalBlocks {
		p := 0
		extPos[id] = &p
	}

	records := make([]cramRecord, 0, sh.numRecords)
	prevAlignPos := sh.alignmentStart

	for i := int32(0); i < sh.numRecords; i++ {
		rec, err := decodeRecord(core, externalBlocks, extPos, ch, sh, prevAlignPos, refs)
		if err != nil {
			return nil, fmt.Errorf("decoding record %d: %w", i, err)
		}
		prevAlignPos = rec.alignPos
		records = append(records, rec)
	}

	return records, nil
}

func decodeRecord(
	core *bitReader,
	external map[int32][]byte,
	extPos map[int32]*int,
	ch *compressionHeader,
	sh *sliceHeader,
	prevAlignPos int32,
	refs []refInfo,
) (cramRecord, error) {
	var rec cramRecord

	// BF: BAM bit flags
	bfCodec, err := ch.getIntCodec("BF")
	if err != nil {
		return rec, err
	}
	rec.bamFlags, err = bfCodec.decodeInt(core, external, extPos)
	if err != nil {
		return rec, fmt.Errorf("decoding BF: %w", err)
	}

	// CF: CRAM compression flags
	cfCodec, err := ch.getIntCodec("CF")
	if err != nil {
		return rec, err
	}
	rec.cramFlags, err = cfCodec.decodeInt(core, external, extPos)
	if err != nil {
		return rec, fmt.Errorf("decoding CF: %w", err)
	}

	// RI: Reference ID (only for multi-ref slices)
	if sh.refSeqID == -2 {
		riCodec, err := ch.getIntCodec("RI")
		if err != nil {
			return rec, err
		}
		rec.refID, err = riCodec.decodeInt(core, external, extPos)
		if err != nil {
			return rec, fmt.Errorf("decoding RI: %w", err)
		}
	} else {
		rec.refID = sh.refSeqID
	}

	// RL: Read length
	rlCodec, err := ch.getIntCodec("RL")
	if err != nil {
		return rec, err
	}
	rec.readLen, err = rlCodec.decodeInt(core, external, extPos)
	if err != nil {
		return rec, fmt.Errorf("decoding RL: %w", err)
	}

	// AP: Alignment position
	apCodec, err := ch.getIntCodec("AP")
	if err != nil {
		return rec, err
	}
	ap, err := apCodec.decodeInt(core, external, extPos)
	if err != nil {
		return rec, fmt.Errorf("decoding AP: %w", err)
	}
	if ch.apDelta {
		rec.alignPos = prevAlignPos + ap
	} else {
		rec.alignPos = ap
	}

	// RG: Read group
	rgCodec, err := ch.getIntCodec("RG")
	if err != nil {
		return rec, err
	}
	rec.readGroup, err = rgCodec.decodeInt(core, external, extPos)
	if err != nil {
		return rec, fmt.Errorf("decoding RG: %w", err)
	}

	// RN: Read name (if preserved)
	if ch.readNamesPreserved {
		rnCodec, err := ch.getByteArrayCodec("RN")
		if err != nil {
			return rec, err
		}
		nameBytes, err := rnCodec.decodeByteArray(core, external, extPos)
		if err != nil {
			return rec, fmt.Errorf("decoding RN: %w", err)
		}
		rec.readName = string(nameBytes)
	}

	// Mate info
	if rec.cramFlags&0x2 != 0 {
		// Detached: mate info stored inline
		mfCodec, err := ch.getIntCodec("MF")
		if err != nil {
			return rec, err
		}
		rec.mateFlags, err = mfCodec.decodeInt(core, external, extPos)
		if err != nil {
			return rec, fmt.Errorf("decoding MF: %w", err)
		}

		if !ch.readNamesPreserved {
			rnCodec, err := ch.getByteArrayCodec("RN")
			if err != nil {
				return rec, err
			}
			nameBytes, err := rnCodec.decodeByteArray(core, external, extPos)
			if err != nil {
				return rec, fmt.Errorf("decoding RN (detached): %w", err)
			}
			rec.readName = string(nameBytes)
		}

		nsCodec, err := ch.getIntCodec("NS")
		if err != nil {
			return rec, err
		}
		rec.mateRefID, err = nsCodec.decodeInt(core, external, extPos)
		if err != nil {
			return rec, fmt.Errorf("decoding NS: %w", err)
		}

		npCodec, err := ch.getIntCodec("NP")
		if err != nil {
			return rec, err
		}
		rec.matePos, err = npCodec.decodeInt(core, external, extPos)
		if err != nil {
			return rec, fmt.Errorf("decoding NP: %w", err)
		}

		tsCodec, err := ch.getIntCodec("TS")
		if err != nil {
			return rec, err
		}
		rec.templateSize, err = tsCodec.decodeInt(core, external, extPos)
		if err != nil {
			return rec, fmt.Errorf("decoding TS: %w", err)
		}
	} else if rec.cramFlags&0x4 != 0 {
		// Mate downstream
		nfCodec, err := ch.getIntCodec("NF")
		if err != nil {
			return rec, err
		}
		rec.recordsToNextFrag, err = nfCodec.decodeInt(core, external, extPos)
		if err != nil {
			return rec, fmt.Errorf("decoding NF: %w", err)
		}
	}

	// TL: Tag dictionary index
	tlCodec, err := ch.getIntCodec("TL")
	if err != nil {
		return rec, err
	}
	rec.tagDictIdx, err = tlCodec.decodeInt(core, external, extPos)
	if err != nil {
		return rec, fmt.Errorf("decoding TL: %w", err)
	}

	// Decode tags
	if int(rec.tagDictIdx) < len(ch.tagDictionary) {
		tagCombo := ch.tagDictionary[rec.tagDictIdx]
		rec.tagValues = make(map[int32][]byte, len(tagCombo))
		for _, tk := range tagCombo {
			key := tagKeyToITF8(tk)
			enc, ok := ch.tagEncodings[key]
			if !ok {
				return rec, fmt.Errorf("no tag encoding for %c%c:%c", tk.id[0], tk.id[1], tk.typ)
			}
			bac, ok := enc.(byteArrayCodec)
			if !ok {
				return rec, fmt.Errorf("tag encoding for %c%c is not byteArrayCodec", tk.id[0], tk.id[1])
			}
			val, err := bac.decodeByteArray(core, external, extPos)
			if err != nil {
				return rec, fmt.Errorf("decoding tag %c%c: %w", tk.id[0], tk.id[1], err)
			}
			rec.tagValues[key] = val
		}
	}

	// Mapped read features or unmapped bases
	isUnmapped := rec.bamFlags&0x4 != 0
	if !isUnmapped {
		// FN: Number of features
		fnCodec, err := ch.getIntCodec("FN")
		if err != nil {
			return rec, err
		}
		numFeatures, err := fnCodec.decodeInt(core, external, extPos)
		if err != nil {
			return rec, fmt.Errorf("decoding FN: %w", err)
		}

		rec.features = make([]cramFeature, 0, numFeatures)
		prevFeatPos := int32(0)
		for j := int32(0); j < numFeatures; j++ {
			feat, err := decodeFeature(core, external, extPos, ch, prevFeatPos)
			if err != nil {
				return rec, fmt.Errorf("decoding feature %d: %w", j, err)
			}
			prevFeatPos = feat.pos
			rec.features = append(rec.features, feat)
		}

		// MQ: Mapping quality
		mqCodec, err := ch.getIntCodec("MQ")
		if err != nil {
			return rec, err
		}
		rec.mapQ, err = mqCodec.decodeInt(core, external, extPos)
		if err != nil {
			return rec, fmt.Errorf("decoding MQ: %w", err)
		}

		// QS: Quality scores (if CF & 0x1)
		if rec.cramFlags&0x1 != 0 {
			qsCodec, err := ch.getByteCodec("QS")
			if err != nil {
				return rec, err
			}
			rec.qualScores = make([]byte, rec.readLen)
			for j := int32(0); j < rec.readLen; j++ {
				rec.qualScores[j], err = qsCodec.decodeByte(core, external, extPos)
				if err != nil {
					return rec, fmt.Errorf("decoding QS[%d]: %w", j, err)
				}
			}
		}
	} else {
		// Unmapped read: BA bases
		baCodec, err := ch.getByteCodec("BA")
		if err != nil {
			return rec, err
		}
		rec.bases = make([]byte, rec.readLen)
		for j := int32(0); j < rec.readLen; j++ {
			rec.bases[j], err = baCodec.decodeByte(core, external, extPos)
			if err != nil {
				return rec, fmt.Errorf("decoding BA[%d]: %w", j, err)
			}
		}

		// QS: Quality scores (if CF & 0x1)
		if rec.cramFlags&0x1 != 0 {
			qsCodec, err := ch.getByteCodec("QS")
			if err != nil {
				return rec, err
			}
			rec.qualScores = make([]byte, rec.readLen)
			for j := int32(0); j < rec.readLen; j++ {
				rec.qualScores[j], err = qsCodec.decodeByte(core, external, extPos)
				if err != nil {
					return rec, fmt.Errorf("decoding QS[%d]: %w", j, err)
				}
			}
		}
	}

	return rec, nil
}

func decodeFeature(
	core *bitReader,
	external map[int32][]byte,
	extPos map[int32]*int,
	ch *compressionHeader,
	prevPos int32,
) (cramFeature, error) {
	var feat cramFeature

	// FC: Feature code
	fcCodec, err := ch.getByteCodec("FC")
	if err != nil {
		return feat, err
	}
	feat.code, err = fcCodec.decodeByte(core, external, extPos)
	if err != nil {
		return feat, fmt.Errorf("decoding FC: %w", err)
	}

	// FP: Feature position (delta from previous)
	fpCodec, err := ch.getIntCodec("FP")
	if err != nil {
		return feat, err
	}
	fp, err := fpCodec.decodeInt(core, external, extPos)
	if err != nil {
		return feat, fmt.Errorf("decoding FP: %w", err)
	}
	feat.pos = prevPos + fp

	switch feat.code {
	case 'B': // Base + quality
		baCodec, err := ch.getByteCodec("BA")
		if err != nil {
			return feat, err
		}
		feat.base, err = baCodec.decodeByte(core, external, extPos)
		if err != nil {
			return feat, err
		}
		qsCodec, err := ch.getByteCodec("QS")
		if err != nil {
			return feat, err
		}
		feat.qual, err = qsCodec.decodeByte(core, external, extPos)
		if err != nil {
			return feat, err
		}

	case 'X': // Substitution
		bsCodec, err := ch.getByteCodec("BS")
		if err != nil {
			return feat, err
		}
		feat.subCode, err = bsCodec.decodeByte(core, external, extPos)
		if err != nil {
			return feat, err
		}

	case 'I': // Insertion (multi-base)
		inCodec, err := ch.getByteArrayCodec("IN")
		if err != nil {
			return feat, err
		}
		feat.bases, err = inCodec.decodeByteArray(core, external, extPos)
		if err != nil {
			return feat, err
		}

	case 'i': // Single inserted base
		baCodec, err := ch.getByteCodec("BA")
		if err != nil {
			return feat, err
		}
		feat.base, err = baCodec.decodeByte(core, external, extPos)
		if err != nil {
			return feat, err
		}

	case 'D': // Deletion
		dlCodec, err := ch.getIntCodec("DL")
		if err != nil {
			return feat, err
		}
		feat.length, err = dlCodec.decodeInt(core, external, extPos)
		if err != nil {
			return feat, err
		}

	case 'S': // Soft clip
		scCodec, err := ch.getByteArrayCodec("SC")
		if err != nil {
			return feat, err
		}
		feat.bases, err = scCodec.decodeByteArray(core, external, extPos)
		if err != nil {
			return feat, err
		}

	case 'N': // Reference skip
		rsCodec, err := ch.getIntCodec("RS")
		if err != nil {
			return feat, err
		}
		feat.length, err = rsCodec.decodeInt(core, external, extPos)
		if err != nil {
			return feat, err
		}

	case 'P': // Padding
		pdCodec, err := ch.getIntCodec("PD")
		if err != nil {
			return feat, err
		}
		feat.length, err = pdCodec.decodeInt(core, external, extPos)
		if err != nil {
			return feat, err
		}

	case 'H': // Hard clip
		hcCodec, err := ch.getIntCodec("HC")
		if err != nil {
			return feat, err
		}
		feat.length, err = hcCodec.decodeInt(core, external, extPos)
		if err != nil {
			return feat, err
		}

	case 'b': // Base stretch
		bbCodec, err := ch.getByteArrayCodec("BB")
		if err != nil {
			return feat, err
		}
		feat.bases, err = bbCodec.decodeByteArray(core, external, extPos)
		if err != nil {
			return feat, err
		}

	case 'q': // Quality stretch
		qqCodec, err := ch.getByteArrayCodec("QQ")
		if err != nil {
			return feat, err
		}
		feat.quals, err = qqCodec.decodeByteArray(core, external, extPos)
		if err != nil {
			return feat, err
		}

	case 'Q': // Single quality score
		qsCodec, err := ch.getByteCodec("QS")
		if err != nil {
			return feat, err
		}
		feat.qual, err = qsCodec.decodeByte(core, external, extPos)
		if err != nil {
			return feat, err
		}

	default:
		return feat, fmt.Errorf("unknown feature code: %c (0x%02x)", feat.code, feat.code)
	}

	return feat, nil
}

// refInfo holds reference name and length (mirrors bamRefInfo from htsio).
type refInfo struct {
	name   string
	length int
}

// reconstructSequence rebuilds the read sequence from features and the reference.
// refOffset is subtracted from reference positions (non-zero for embedded refs).
func reconstructSequence(rec *cramRecord, ch *compressionHeader, refSeq []byte, refOffset int) string {
	if rec.bamFlags&0x4 != 0 {
		// Unmapped
		if rec.readLen == 0 {
			return "*"
		}
		return string(rec.bases)
	}

	if rec.cramFlags&0x8 != 0 {
		// Sequence unknown
		return "*"
	}

	if debugReconstruct {
		fmt.Fprintf(os.Stderr, "DEBUG reconstruct: readLen=%d alignPos=%d numFeatures=%d\n",
			rec.readLen, rec.alignPos, len(rec.features))
		for i, f := range rec.features {
			fmt.Fprintf(os.Stderr, "  feat[%d]: code=%c pos=%d subCode=%d base=%c\n",
				i, f.code, f.pos, f.subCode, f.base)
		}
	}

	seq := make([]byte, rec.readLen)
	readPos := 0   // 0-based position in the read
	refPos := int(rec.alignPos) - 1 - refOffset // 0-based position in refSeq

	featureIdx := 0

	for readPos < int(rec.readLen) {
		if featureIdx < len(rec.features) {
			feat := rec.features[featureIdx]
			featReadPos := int(feat.pos) - 1 // convert to 0-based

			// Copy reference bases up to the feature position
			for readPos < featReadPos && readPos < int(rec.readLen) {
				if refPos >= 0 && refPos < len(refSeq) {
					seq[readPos] = refSeq[refPos]
				} else {
					seq[readPos] = 'N'
				}
				readPos++
				refPos++
			}

			featureIdx++

			switch feat.code {
			case 'X': // Substitution
				if readPos >= int(rec.readLen) {
					break
				}
				var refBase byte
				if refPos >= 0 && refPos < len(refSeq) {
					refBase = refSeq[refPos]
				} else {
					refBase = 'N'
				}
				seq[readPos] = ch.lookupSubstitution(refBase, feat.subCode)
				readPos++
				refPos++

			case 'B': // Base + quality (replaces ref base)
				if readPos >= int(rec.readLen) {
					break
				}
				seq[readPos] = feat.base
				readPos++
				refPos++

			case 'I': // Insertion
				for _, b := range feat.bases {
					if readPos < int(rec.readLen) {
						seq[readPos] = b
						readPos++
					}
				}

			case 'i': // Single insert
				if readPos < int(rec.readLen) {
					seq[readPos] = feat.base
					readPos++
				}

			case 'D': // Deletion (skip ref bases)
				refPos += int(feat.length)

			case 'N': // Reference skip
				refPos += int(feat.length)

			case 'S': // Soft clip
				for _, b := range feat.bases {
					if readPos < int(rec.readLen) {
						seq[readPos] = b
						readPos++
					}
				}

			case 'b': // Base stretch
				for _, b := range feat.bases {
					if readPos < int(rec.readLen) {
						seq[readPos] = b
						readPos++
						refPos++
					}
				}

			case 'P': // Padding — no change to read or ref
			case 'H': // Hard clip — no change to read or ref
			case 'Q', 'q': // Quality only — no change to seq
			}
		} else {
			// No more features — copy from reference
			if refPos >= 0 && refPos < len(refSeq) {
				seq[readPos] = refSeq[refPos]
			} else {
				seq[readPos] = 'N'
			}
			readPos++
			refPos++
		}
	}

	return string(seq)
}

// reconstructCigar builds a CIGAR string from the features.
func reconstructCigar(rec *cramRecord) string {
	if rec.bamFlags&0x4 != 0 {
		return "*"
	}

	type cigarOp struct {
		length int
		op     byte
	}

	var ops []cigarOp
	addOp := func(length int, op byte) {
		if length == 0 {
			return
		}
		if len(ops) > 0 && ops[len(ops)-1].op == op {
			ops[len(ops)-1].length += length
		} else {
			ops = append(ops, cigarOp{length, op})
		}
	}

	readPos := 0
	prevFeatReadPos := 0

	for _, feat := range rec.features {
		featReadPos := int(feat.pos) - 1 // 0-based

		// Match bases between features
		matchLen := featReadPos - prevFeatReadPos
		if matchLen < 0 {
			matchLen = 0
		}

		switch feat.code {
		case 'X': // Substitution — still consumes ref+read (M op)
			addOp(matchLen, 'M')
			addOp(1, 'M')
			readPos = featReadPos + 1
			prevFeatReadPos = readPos

		case 'B': // Base replacement — consumes ref+read (M op)
			addOp(matchLen, 'M')
			addOp(1, 'M')
			readPos = featReadPos + 1
			prevFeatReadPos = readPos

		case 'I': // Insertion
			addOp(matchLen, 'M')
			addOp(len(feat.bases), 'I')
			readPos = featReadPos + len(feat.bases)
			prevFeatReadPos = readPos

		case 'i': // Single insert
			addOp(matchLen, 'M')
			addOp(1, 'I')
			readPos = featReadPos + 1
			prevFeatReadPos = readPos

		case 'D': // Deletion
			addOp(matchLen, 'M')
			addOp(int(feat.length), 'D')
			readPos = featReadPos
			prevFeatReadPos = readPos

		case 'N': // Reference skip (intron)
			addOp(matchLen, 'M')
			addOp(int(feat.length), 'N')
			readPos = featReadPos
			prevFeatReadPos = readPos

		case 'S': // Soft clip
			addOp(matchLen, 'M')
			addOp(len(feat.bases), 'S')
			readPos = featReadPos + len(feat.bases)
			prevFeatReadPos = readPos

		case 'P': // Padding
			addOp(matchLen, 'M')
			addOp(int(feat.length), 'P')
			readPos = featReadPos
			prevFeatReadPos = readPos

		case 'H': // Hard clip
			addOp(matchLen, 'M')
			addOp(int(feat.length), 'H')
			readPos = featReadPos
			prevFeatReadPos = readPos

		case 'b': // Base stretch — matches with explicit bases
			addOp(matchLen, 'M')
			addOp(len(feat.bases), 'M')
			readPos = featReadPos + len(feat.bases)
			prevFeatReadPos = readPos

		case 'Q', 'q': // Quality modifications — don't affect CIGAR
			// readPos stays same, no CIGAR op
		}
	}

	// Trailing match
	remaining := int(rec.readLen) - prevFeatReadPos
	if remaining > 0 {
		addOp(remaining, 'M')
	}

	if len(ops) == 0 {
		return fmt.Sprintf("%dM", rec.readLen)
	}

	var sb strings.Builder
	for _, op := range ops {
		sb.WriteString(strconv.Itoa(op.length))
		sb.WriteByte(op.op)
	}
	return sb.String()
}

// reconstructQual builds the quality string from the record.
func reconstructQual(rec *cramRecord) string {
	if rec.qualScores == nil {
		return "*"
	}
	// Quality scores in CRAM are raw Phred values; SAM uses Phred+33.
	buf := make([]byte, len(rec.qualScores))
	allFF := true
	for i, q := range rec.qualScores {
		if q != 0xFF {
			allFF = false
		}
		buf[i] = q + 33
	}
	if allFF {
		return "*"
	}
	return string(buf)
}

// cramTagToSAM converts a CRAM tag value (raw bytes) to SAM tag string format.
// The type byte determines the interpretation of the raw bytes.
func cramTagToSAM(tk tagKey, raw []byte) (samType byte, samValue string) {
	switch tk.typ {
	case 'A': // character
		if len(raw) >= 1 {
			return 'A', string(raw[0:1])
		}
		return 'A', ""

	case 'c': // int8
		if len(raw) >= 1 {
			return 'i', strconv.Itoa(int(int8(raw[0])))
		}
	case 'C': // uint8
		if len(raw) >= 1 {
			return 'i', strconv.Itoa(int(raw[0]))
		}
	case 's': // int16
		if len(raw) >= 2 {
			v := int16(raw[0]) | int16(raw[1])<<8
			return 'i', strconv.Itoa(int(v))
		}
	case 'S': // uint16
		if len(raw) >= 2 {
			v := uint16(raw[0]) | uint16(raw[1])<<8
			return 'i', strconv.Itoa(int(v))
		}
	case 'i': // int32
		if len(raw) >= 4 {
			v := int32(raw[0]) | int32(raw[1])<<8 | int32(raw[2])<<16 | int32(raw[3])<<24
			return 'i', strconv.Itoa(int(v))
		}
	case 'I': // uint32
		if len(raw) >= 4 {
			v := uint32(raw[0]) | uint32(raw[1])<<8 | uint32(raw[2])<<16 | uint32(raw[3])<<24
			return 'i', strconv.FormatUint(uint64(v), 10)
		}
	case 'f': // float32
		if len(raw) >= 4 {
			bits := uint32(raw[0]) | uint32(raw[1])<<8 | uint32(raw[2])<<16 | uint32(raw[3])<<24
			return 'f', strconv.FormatFloat(float64(math.Float32frombits(bits)), 'g', -1, 32)
		}
	case 'Z': // string
		return 'Z', string(raw)
	case 'H': // hex
		return 'H', string(raw)
	case 'B': // array — raw bytes include subtype + count + values
		return 'B', string(raw)
	}
	return 'Z', string(raw) // fallback
}
