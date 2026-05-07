package cram

import (
	"bytes"
	"fmt"
	"io"
	"iter"
	"os"
	"strings"

	"github.com/compgen-io/cgkit/htsio"
)

func init() {
	htsio.RegisterReader(htsio.ReaderRegistration{
		Detect: func(magic []byte) bool {
			return bytes.HasPrefix(magic, []byte("CRAM"))
		},
		NewFromFile: func(filename string, opts *htsio.SamReaderOpts) (htsio.SamReader, error) {
			return NewReader(filename, "")
		},
		NewFromStream: func(r io.ReadCloser, opts *htsio.SamReaderOpts) (htsio.SamReader, error) {
			return NewReaderFromStream(r, "", "")
		},
	})
}

// Reader reads CRAM files and implements htsio.SamReader.
type Reader struct {
	r        io.Reader
	src      io.ReadCloser
	filename string
	fileDef  *fileDefinition
	hdr      *htsio.SamHeader
	refs     []refInfo
	refMap   map[string]int // ref name → index
	refProv  *referenceProvider
	idx      *craiIndex // lazily loaded CRAI index for Query
	queryFh  *os.File   // separate file handle for Query seeks
}

// NewReader creates a CRAM reader from a file path.
// refPath is the path to the reference FASTA. If empty, the reader
// will attempt to find the reference from the UR field in @SQ header lines.
func NewReader(filename string, refPath string) (*Reader, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	return NewReaderFromStream(f, filename, refPath)
}

// NewReaderFromStream creates a CRAM reader from an io.ReadCloser.
// The stream must be positioned at the start of a CRAM file (possibly
// with peeked bytes prepended). filename is used for index lookups.
func NewReaderFromStream(rc io.ReadCloser, filename string, refPath string) (*Reader, error) {
	cr := &Reader{
		r:        rc,
		src:      rc,
		filename: filename,
	}

	// Read file definition.
	var err error
	cr.fileDef, err = readFileDefinition(cr.r)
	if err != nil {
		rc.Close()
		return nil, fmt.Errorf("cram: %w", err)
	}

	// Read header container.
	if err := cr.readHeaderContainer(); err != nil {
		rc.Close()
		return nil, fmt.Errorf("cram: reading header: %w", err)
	}

	// Set up reference provider.
	if refPath == "" {
		refPath = cr.findReferenceFromHeader()
	}
	if refPath != "" {
		cr.refProv = newReferenceProvider(refPath)
	}

	return cr, nil
}

// version returns the major version of the CRAM file.
func (cr *Reader) version() byte {
	return cr.fileDef.Major
}

// Header returns the parsed SAM header.
func (cr *Reader) Header() (*htsio.SamHeader, error) {
	return cr.hdr, nil
}

// Close releases resources.
func (cr *Reader) Close() error {
	if cr.queryFh != nil {
		cr.queryFh.Close()
	}
	if cr.src != nil {
		return cr.src.Close()
	}
	return nil
}

// Query returns an iterator over records overlapping the 0-based half-open
// region [start, end) on the given reference. Requires a .crai index file.
func (cr *Reader) Query(ref string, start, end int) (iter.Seq2[*htsio.SamRecord, error], error) {
	if cr.filename == "" {
		return nil, fmt.Errorf("cram: Query requires a file-backed reader")
	}

	seqID, ok := cr.refMap[ref]
	if !ok {
		return nil, fmt.Errorf("cram: unknown reference %q", ref)
	}

	// Lazily load the CRAI index.
	if cr.idx == nil {
		craiPath := cr.filename + ".crai"
		idx, err := loadCRAI(craiPath)
		if err != nil {
			return nil, fmt.Errorf("cram: loading CRAI index: %w", err)
		}
		cr.idx = idx
	}

	// Find overlapping slices.
	entries := cr.idx.query(seqID, start, end)
	if len(entries) == 0 {
		return func(yield func(*htsio.SamRecord, error) bool) {}, nil
	}

	// Lazily open a separate file handle for queries.
	if cr.queryFh == nil {
		f, err := os.Open(cr.filename)
		if err != nil {
			return nil, err
		}
		cr.queryFh = f
	}

	return cr.iterCraiEntries(entries, seqID, start, end), nil
}

// iterCraiEntries returns an iterator that reads records from the CRAI entries,
// filtering to those overlapping [start, end) on seqID.
func (cr *Reader) iterCraiEntries(entries []craiEntry, seqID, start, end int) iter.Seq2[*htsio.SamRecord, error] {
	return func(yield func(*htsio.SamRecord, error) bool) {
		for _, entry := range entries {
			// Seek to the container.
			if _, err := cr.queryFh.Seek(entry.containerOffset, io.SeekStart); err != nil {
				yield(nil, fmt.Errorf("cram: seeking to container: %w", err))
				return
			}

			// Read container header.
			ch, err := readContainerHeader(cr.queryFh, cr.version())
			if err != nil {
				yield(nil, fmt.Errorf("cram: reading container header: %w", err))
				return
			}
			if ch.isEOF() {
				continue
			}

			// Read compression header (first block).
			compHdrBlock, err := readBlock(cr.queryFh, cr.version())
			if err != nil {
				yield(nil, fmt.Errorf("cram: reading compression header: %w", err))
				return
			}
			if compHdrBlock.contentType != blockContentCompressionHeader {
				yield(nil, fmt.Errorf("cram: expected compression header, got type %d", compHdrBlock.contentType))
				return
			}
			compHdr, err := readCompressionHeader(compHdrBlock.data)
			if err != nil {
				yield(nil, fmt.Errorf("cram: parsing compression header: %w", err))
				return
			}

			// Read remaining blocks (slices).
			remainingBlocks := ch.NumBlocks - 1
			for remainingBlocks > 0 {
				sliceHdrBlock, err := readBlock(cr.queryFh, cr.version())
				if err != nil {
					yield(nil, fmt.Errorf("cram: reading slice header: %w", err))
					return
				}
				remainingBlocks--

				if sliceHdrBlock.contentType != blockContentSliceHeader {
					yield(nil, fmt.Errorf("cram: expected slice header, got type %d", sliceHdrBlock.contentType))
					return
				}

				sh, err := readSliceHeader(sliceHdrBlock.data)
				if err != nil {
					yield(nil, fmt.Errorf("cram: parsing slice header: %w", err))
					return
				}

				coreData := []byte{}
				externalBlocks := make(map[int32][]byte)
				for i := int32(0); i < sh.numBlocks; i++ {
					blk, err := readBlock(cr.queryFh, cr.version())
					if err != nil {
						yield(nil, fmt.Errorf("cram: reading slice block %d: %w", i, err))
						return
					}
					remainingBlocks--
					switch blk.contentType {
					case blockContentCoreData:
						coreData = blk.data
					case blockContentExternalData:
						externalBlocks[blk.contentID] = blk.data
					}
				}

				// Load reference.
				var refSeq []byte
				if sh.refSeqID >= 0 && cr.refProv != nil {
					if int(sh.refSeqID) < len(cr.refs) {
						if seq, err := cr.refProv.getSequence(cr.refs[sh.refSeqID].name); err == nil {
							refSeq = seq
						}
					}
				}

				records, err := decodeSliceRecords(sh, compHdr, coreData, externalBlocks, cr.refs, nil)
				if err != nil {
					yield(nil, fmt.Errorf("cram: decoding slice records: %w", err))
					return
				}

				for i := range records {
					rec := &records[i]
					recRefSeq := refSeq
					if sh.refSeqID == -2 && cr.refProv != nil && rec.refID >= 0 {
						if int(rec.refID) < len(cr.refs) {
							if seq, err := cr.refProv.getSequence(cr.refs[rec.refID].name); err == nil {
								recRefSeq = seq
							}
						}
					}

					// Filter: record must be on the queried reference and overlap [start, end).
					if int(rec.refID) != seqID {
						continue
					}
					recStart := int(rec.alignPos) - 1 // convert 1-based to 0-based
					recEnd := recStart + int(rec.readLen)
					if recStart >= end || recEnd <= start {
						continue
					}

					samRec := cr.cramToSam(rec, compHdr, recRefSeq)
					if !yield(samRec, nil) {
						return
					}
				}
			}
		}
	}
}

// Records returns an iterator over all records in the CRAM file.
func (cr *Reader) Records() iter.Seq2[*htsio.SamRecord, error] {
	return func(yield func(*htsio.SamRecord, error) bool) {
		for {
			// Read next container header.
			ch, err := readContainerHeader(cr.r, cr.version())
			if err != nil {
				if err == io.EOF || err == io.ErrUnexpectedEOF {
					return
				}
				yield(nil, fmt.Errorf("cram: reading container: %w", err))
				return
			}

			if ch.isEOF() {
				return
			}

			// Read all blocks in the container.
			if !cr.processContainer(ch, yield) {
				return
			}
		}
	}
}

// processContainer reads and processes a single data container.
// Returns false if yield returned false (caller should stop).
func (cr *Reader) processContainer(ch *containerHeader, yield func(*htsio.SamRecord, error) bool) bool {
	// First block: compression header.
	compHdrBlock, err := readBlock(cr.r, cr.version())
	if err != nil {
		return yield(nil, fmt.Errorf("cram: reading compression header block: %w", err))
	}
	if compHdrBlock.contentType != blockContentCompressionHeader {
		return yield(nil, fmt.Errorf("cram: expected compression header block, got type %d", compHdrBlock.contentType))
	}

	compHdr, err := readCompressionHeader(compHdrBlock.data)
	if err != nil {
		return yield(nil, fmt.Errorf("cram: parsing compression header: %w", err))
	}

	// Remaining blocks are organized into slices.
	// Each slice starts with a slice header block, followed by a core data block,
	// then external data blocks.
	remainingBlocks := ch.NumBlocks - 1 // minus the compression header

	for remainingBlocks > 0 {
		// Slice header block.
		sliceHdrBlock, err := readBlock(cr.r, cr.version())
		if err != nil {
			return yield(nil, fmt.Errorf("cram: reading slice header block: %w", err))
		}
		remainingBlocks--

		if sliceHdrBlock.contentType != blockContentSliceHeader {
			return yield(nil, fmt.Errorf("cram: expected slice header block, got type %d", sliceHdrBlock.contentType))
		}

		sh, err := readSliceHeader(sliceHdrBlock.data)
		if err != nil {
			return yield(nil, fmt.Errorf("cram: parsing slice header: %w", err))
		}

		// Read the core block and external blocks for this slice.
		coreData := []byte{}
		externalBlocks := make(map[int32][]byte)

		for i := int32(0); i < sh.numBlocks; i++ {
			blk, err := readBlock(cr.r, cr.version())
			if err != nil {
				return yield(nil, fmt.Errorf("cram: reading slice block %d: %w", i, err))
			}
			remainingBlocks--

			switch blk.contentType {
			case blockContentCoreData:
				coreData = blk.data
			case blockContentExternalData:
				externalBlocks[blk.contentID] = blk.data
			}
		}

		// Load reference sequence for this slice.
		var refSeq []byte
		if sh.refSeqID >= 0 && cr.refProv != nil {
			refName := ""
			if int(sh.refSeqID) < len(cr.refs) {
				refName = cr.refs[sh.refSeqID].name
			}
			if refName != "" {
				refSeq, err = cr.refProv.getSequence(refName)
				if err != nil {
					// Non-fatal: sequence reconstruction will use 'N' for missing ref bases
					refSeq = nil
				}
			}
		}

		// Decode records.
		records, err := decodeSliceRecords(sh, compHdr, coreData, externalBlocks, cr.refs, nil)
		if err != nil {
			return yield(nil, fmt.Errorf("cram: decoding slice records: %w", err))
		}

		// Convert to SamRecords and yield.
		for i := range records {
			recRefSeq := refSeq
			// For multi-ref slices, load the correct reference per record.
			if sh.refSeqID == -2 && cr.refProv != nil && records[i].refID >= 0 {
				if int(records[i].refID) < len(cr.refs) {
					rn := cr.refs[records[i].refID].name
					if seq, err := cr.refProv.getSequence(rn); err == nil {
						recRefSeq = seq
					}
				}
			}
			samRec := cr.cramToSam(&records[i], compHdr, recRefSeq)
			if !yield(samRec, nil) {
				return false
			}
		}
	}

	return true
}

// cramToSam converts a CRAM record to a SamRecord.
func (cr *Reader) cramToSam(rec *cramRecord, ch *compressionHeader, refSeq []byte) *htsio.SamRecord {
	// Reference name
	refName := "*"
	if rec.refID >= 0 && int(rec.refID) < len(cr.refs) {
		refName = cr.refs[rec.refID].name
	}

	// Mate reference name
	mateRefName := "*"
	if rec.cramFlags&0x2 != 0 { // detached
		if rec.mateRefID >= 0 && int(rec.mateRefID) < len(cr.refs) {
			mateRefName = cr.refs[rec.mateRefID].name
		}
		if rec.mateRefID == rec.refID && rec.refID >= 0 {
			mateRefName = "="
		}
	}

	// Sequence
	seq := reconstructSequence(rec, ch, refSeq)

	// CIGAR
	cigar := reconstructCigar(rec)

	// Quality
	qual := reconstructQual(rec)

	// Read group tag
	tags := make(map[string]htsio.SamTag)

	// Add read group
	if rec.readGroup >= 0 {
		rgs := cr.hdr.ReadGroups()
		if int(rec.readGroup) < len(rgs) {
			tags["RG"] = htsio.SamTag{Type: 'Z', Value: rgs[rec.readGroup]}
		}
	}

	// Add other tags from CRAM tag values
	if int(rec.tagDictIdx) < len(ch.tagDictionary) {
		tagCombo := ch.tagDictionary[rec.tagDictIdx]
		for _, tk := range tagCombo {
			key := tagKeyToITF8(tk)
			raw, ok := rec.tagValues[key]
			if !ok {
				continue
			}
			tagName := string(tk.id[:])
			samType, samValue := cramTagToSAM(tk, raw)
			tags[tagName] = htsio.SamTag{Type: samType, Value: samValue}
		}
	}

	return &htsio.SamRecord{
		ReadName:  rec.readName,
		Flag:      int(rec.bamFlags),
		RefName:   refName,
		Pos:       int(rec.alignPos),
		MapQ:      int(rec.mapQ),
		Cigar:     cigar,
		RefNext:   mateRefName,
		PosNext:   int(rec.matePos),
		InsertLen: int(rec.templateSize),
		Seq:       seq,
		Qual:      qual,
		Tags:      tags,
	}
}

// readHeaderContainer reads the first (header) container from the CRAM file.
func (cr *Reader) readHeaderContainer() error {
	ch, err := readContainerHeader(cr.r, cr.version())
	if err != nil {
		return fmt.Errorf("reading header container: %w", err)
	}

	// The header container should have one block containing the SAM header text.
	blk, err := readBlock(cr.r, cr.version())
	if err != nil {
		return fmt.Errorf("reading header block: %w", err)
	}

	cr.hdr = htsio.NewSamHeader()

	// The FILE_HEADER block contains an int32 length prefix followed by the SAM text.
	headerData := blk.data
	if len(headerData) >= 4 {
		headerLen := int(headerData[0]) | int(headerData[1])<<8 | int(headerData[2])<<16 | int(headerData[3])<<24
		if headerLen > 0 && headerLen+4 <= len(headerData) {
			headerData = headerData[4 : 4+headerLen]
		} else {
			headerData = headerData[4:]
		}
	}

	headerText := string(headerData)
	for _, line := range strings.Split(headerText, "\n") {
		line = strings.TrimRight(line, "\r")
		if line != "" {
			cr.hdr.AddLine(line)
		}
	}

	// Build refs from header @SQ lines.
	headerRefs := cr.hdr.References()
	cr.refs = make([]refInfo, len(headerRefs))
	cr.refMap = make(map[string]int, len(headerRefs))
	for i, hr := range headerRefs {
		cr.refs[i] = refInfo{name: hr.Name, length: hr.Length}
		cr.refMap[hr.Name] = i
	}

	// Skip remaining blocks in header container if any.
	for i := int32(1); i < ch.NumBlocks; i++ {
		if _, err := readBlock(cr.r, cr.version()); err != nil {
			return fmt.Errorf("skipping header block %d: %w", i, err)
		}
	}

	return nil
}

// findReferenceFromHeader looks for a UR field in the first @SQ line
// and returns it if it looks like a local file path.
func (cr *Reader) findReferenceFromHeader() string {
	for _, line := range cr.hdr.Lines {
		if !strings.HasPrefix(line, "@SQ\t") {
			continue
		}
		for _, field := range strings.Split(line, "\t")[1:] {
			if strings.HasPrefix(field, "UR:") {
				uri := field[3:]
				// Only use local file paths, not URLs.
				if strings.HasPrefix(uri, "/") || strings.HasPrefix(uri, "file://") {
					path := strings.TrimPrefix(uri, "file://")
					if _, err := os.Stat(path); err == nil {
						return path
					}
				}
			}
		}
		break // only check first @SQ
	}
	return ""
}
