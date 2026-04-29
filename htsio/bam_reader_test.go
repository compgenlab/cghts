package htsio

import (
	"bytes"
	"encoding/binary"
	"io"
	"testing"

	"github.com/compgen-io/cgltk/htsio/bgzf"
)

// writeBamTestFile writes a minimal BAM file to a buffer with one reference
// and the given alignment records (already encoded as raw BAM record bytes).
func writeBamTestFile(headerText string, refs []bamRefInfo, records [][]byte) ([]byte, error) {
	var buf bytes.Buffer
	w := bgzf.NewWriter(&buf)

	// Magic
	w.Write([]byte("BAM\x01"))

	// Header text
	hdr := []byte(headerText)
	binary.Write(w, binary.LittleEndian, int32(len(hdr)))
	w.Write(hdr)

	// Reference sequences
	binary.Write(w, binary.LittleEndian, int32(len(refs)))
	for _, ref := range refs {
		name := []byte(ref.name + "\x00")
		binary.Write(w, binary.LittleEndian, int32(len(name)))
		w.Write(name)
		binary.Write(w, binary.LittleEndian, ref.length)
	}

	// Records
	for _, rec := range records {
		binary.Write(w, binary.LittleEndian, int32(len(rec)))
		w.Write(rec)
	}

	if err := w.Close(); err != nil {
		return nil, err
	}
	return buf.Bytes(), nil
}

// encodeBamRecord builds a raw BAM alignment record from readable fields.
func encodeBamRecord(
	readName string,
	flag uint16,
	refID int32,
	pos int32, // 0-based
	mapq uint8,
	cigarOps []uint32, // packed CIGAR ops
	nextRefID int32,
	nextPos int32,
	tlen int32,
	seq string, // ASCII bases
	qual []byte, // raw Phred values
	auxTags []byte,
) []byte {
	var buf bytes.Buffer

	nameBytes := append([]byte(readName), 0) // NUL-terminated

	// Compute bin from pos and cigar ref length.
	bin := uint16(0) // simplified: use bin 0 for testing

	// Fixed fields
	binary.Write(&buf, binary.LittleEndian, refID)
	binary.Write(&buf, binary.LittleEndian, pos)
	buf.WriteByte(byte(len(nameBytes))) // l_read_name
	buf.WriteByte(mapq)
	binary.Write(&buf, binary.LittleEndian, bin)
	binary.Write(&buf, binary.LittleEndian, uint16(len(cigarOps)))
	binary.Write(&buf, binary.LittleEndian, flag)
	binary.Write(&buf, binary.LittleEndian, int32(len(seq)))
	binary.Write(&buf, binary.LittleEndian, nextRefID)
	binary.Write(&buf, binary.LittleEndian, nextPos)
	binary.Write(&buf, binary.LittleEndian, tlen)

	// Read name
	buf.Write(nameBytes)

	// CIGAR
	for _, op := range cigarOps {
		binary.Write(&buf, binary.LittleEndian, op)
	}

	// Sequence (4-bit packed)
	seqBytes := (len(seq) + 1) / 2
	encoded := make([]byte, seqBytes)
	for i := 0; i < len(seq); i++ {
		code := seqEncode(seq[i])
		if i%2 == 0 {
			encoded[i/2] = code << 4
		} else {
			encoded[i/2] |= code
		}
	}
	buf.Write(encoded)

	// Quality
	if qual == nil {
		qual = make([]byte, len(seq))
		for i := range qual {
			qual[i] = 0xFF
		}
	}
	buf.Write(qual)

	// Aux tags
	buf.Write(auxTags)

	return buf.Bytes()
}

func seqEncode(base byte) byte {
	switch base {
	case '=':
		return 0
	case 'A':
		return 1
	case 'C':
		return 2
	case 'G':
		return 4
	case 'T':
		return 8
	case 'N':
		return 15
	default:
		return 15
	}
}

// packCigar creates a packed CIGAR op (length<<4 | opCode).
func packCigar(length uint32, op byte) uint32 {
	var code uint32
	switch op {
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
	case '=':
		code = 7
	case 'X':
		code = 8
	}
	return length<<4 | code
}

func TestBamReaderBasic(t *testing.T) {
	refs := []bamRefInfo{
		{name: "chr1", length: 100000},
	}
	headerText := "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n"

	// 10M cigar, seq=ACGTACGTAC, qual all 30
	cigar := []uint32{packCigar(10, 'M')}
	qual := make([]byte, 10)
	for i := range qual {
		qual[i] = 30
	}
	rec := encodeBamRecord("read1", 0, 0, 99, 60, cigar, -1, -1, 0, "ACGTACGTAC", qual, nil)
	records := [][]byte{rec}

	data, err := writeBamTestFile(headerText, refs, records)
	if err != nil {
		t.Fatalf("writeBamTestFile: %v", err)
	}

	reader, err := NewBamReader(io.NopCloser(bytes.NewReader(data)))
	if err != nil {
		t.Fatalf("NewBamReader: %v", err)
	}
	defer reader.Close()

	hdr, err := reader.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if hdr == nil {
		t.Fatal("expected non-nil header")
	}

	samRec, err := reader.Next()
	if err != nil {
		t.Fatalf("Next: %v", err)
	}

	if samRec.ReadName != "read1" {
		t.Errorf("ReadName: got %q, want %q", samRec.ReadName, "read1")
	}
	if samRec.RefName != "chr1" {
		t.Errorf("RefName: got %q, want %q", samRec.RefName, "chr1")
	}
	if samRec.Pos != 100 { // BAM 0-based 99 → SAM 1-based 100
		t.Errorf("Pos: got %d, want 100", samRec.Pos)
	}
	if samRec.MapQ != 60 {
		t.Errorf("MapQ: got %d, want 60", samRec.MapQ)
	}
	if samRec.Cigar != "10M" {
		t.Errorf("Cigar: got %q, want %q", samRec.Cigar, "10M")
	}
	if samRec.Seq != "ACGTACGTAC" {
		t.Errorf("Seq: got %q, want %q", samRec.Seq, "ACGTACGTAC")
	}
	if samRec.Flag != 0 {
		t.Errorf("Flag: got %d, want 0", samRec.Flag)
	}

	// Verify quality (BAM Phred 30 → SAM Phred+33 = '?')
	expectedQual := "??????????"
	if samRec.Qual != expectedQual {
		t.Errorf("Qual: got %q, want %q", samRec.Qual, expectedQual)
	}

	// Should be EOF
	_, err = reader.Next()
	if err != io.EOF {
		t.Errorf("expected io.EOF, got %v", err)
	}
}

func TestBamReaderMultipleRecords(t *testing.T) {
	refs := []bamRefInfo{
		{name: "chr1", length: 100000},
		{name: "chr2", length: 50000},
	}
	headerText := "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr2\tLN:50000\n"

	cigar := []uint32{packCigar(5, 'M')}
	rec1 := encodeBamRecord("readA", 0, 0, 0, 40, cigar, -1, -1, 0, "ACGTG", nil, nil)
	rec2 := encodeBamRecord("readB", 16, 1, 999, 30, cigar, -1, -1, 0, "TGCAA", nil, nil)

	data, err := writeBamTestFile(headerText, refs, [][]byte{rec1, rec2})
	if err != nil {
		t.Fatalf("writeBamTestFile: %v", err)
	}

	reader, err := NewBamReader(io.NopCloser(bytes.NewReader(data)))
	if err != nil {
		t.Fatalf("NewBamReader: %v", err)
	}
	defer reader.Close()

	// First record
	r1, err := reader.Next()
	if err != nil {
		t.Fatalf("Next[0]: %v", err)
	}
	if r1.ReadName != "readA" || r1.RefName != "chr1" || r1.Pos != 1 {
		t.Errorf("rec1: name=%q ref=%q pos=%d", r1.ReadName, r1.RefName, r1.Pos)
	}

	// Second record
	r2, err := reader.Next()
	if err != nil {
		t.Fatalf("Next[1]: %v", err)
	}
	if r2.ReadName != "readB" || r2.RefName != "chr2" || r2.Pos != 1000 {
		t.Errorf("rec2: name=%q ref=%q pos=%d", r2.ReadName, r2.RefName, r2.Pos)
	}
	if !r2.IsReverse() {
		t.Error("rec2 should be reverse strand")
	}

	// EOF
	_, err = reader.Next()
	if err != io.EOF {
		t.Errorf("expected io.EOF, got %v", err)
	}
}

func TestBamReaderAuxTags(t *testing.T) {
	refs := []bamRefInfo{{name: "chr1", length: 100000}}
	headerText := "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n"

	cigar := []uint32{packCigar(4, 'M')}

	// Build aux tags: RG:Z:sample1 NM:C:3
	var aux bytes.Buffer
	aux.WriteString("RG")
	aux.WriteByte('Z')
	aux.WriteString("sample1\x00")
	aux.WriteString("NM")
	aux.WriteByte('C') // uint8
	aux.WriteByte(3)

	rec := encodeBamRecord("read1", 0, 0, 0, 60, cigar, -1, -1, 0, "ACGT", nil, aux.Bytes())

	data, err := writeBamTestFile(headerText, refs, [][]byte{rec})
	if err != nil {
		t.Fatalf("writeBamTestFile: %v", err)
	}

	reader, err := NewBamReader(io.NopCloser(bytes.NewReader(data)))
	if err != nil {
		t.Fatalf("NewBamReader: %v", err)
	}
	defer reader.Close()

	r, err := reader.Next()
	if err != nil {
		t.Fatalf("Next: %v", err)
	}

	rg, ok := r.Tags["RG"]
	if !ok {
		t.Fatal("missing RG tag")
	}
	if rg.Type != 'Z' || rg.Value != "sample1" {
		t.Errorf("RG tag: type=%c value=%q", rg.Type, rg.Value)
	}

	nm, ok := r.Tags["NM"]
	if !ok {
		t.Fatal("missing NM tag")
	}
	if nm.Type != 'i' || nm.Value != "3" {
		t.Errorf("NM tag: type=%c value=%q", nm.Type, nm.Value)
	}
}

func TestBamReaderComplexCigar(t *testing.T) {
	refs := []bamRefInfo{{name: "chr1", length: 100000}}
	headerText := "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n"

	// 5S10M2I3M1D4M3S
	cigar := []uint32{
		packCigar(5, 'S'),
		packCigar(10, 'M'),
		packCigar(2, 'I'),
		packCigar(3, 'M'),
		packCigar(1, 'D'),
		packCigar(4, 'M'),
		packCigar(3, 'S'),
	}
	seq := "ACGTAACGTAACGTAACGTAACGTAAC" // 27 bases
	rec := encodeBamRecord("read1", 0, 0, 99, 60, cigar, -1, -1, 0, seq, nil, nil)

	data, err := writeBamTestFile(headerText, refs, [][]byte{rec})
	if err != nil {
		t.Fatalf("writeBamTestFile: %v", err)
	}

	reader, err := NewBamReader(io.NopCloser(bytes.NewReader(data)))
	if err != nil {
		t.Fatalf("NewBamReader: %v", err)
	}
	defer reader.Close()

	r, err := reader.Next()
	if err != nil {
		t.Fatalf("Next: %v", err)
	}

	if r.Cigar != "5S10M2I3M1D4M3S" {
		t.Errorf("Cigar: got %q, want %q", r.Cigar, "5S10M2I3M1D4M3S")
	}
}

func TestBamReaderFlagFilter(t *testing.T) {
	refs := []bamRefInfo{{name: "chr1", length: 100000}}
	headerText := "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n"

	cigar := []uint32{packCigar(4, 'M')}
	rec1 := encodeBamRecord("primary", 0, 0, 0, 60, cigar, -1, -1, 0, "ACGT", nil, nil)
	rec2 := encodeBamRecord("secondary", 0x100, 0, 100, 0, cigar, -1, -1, 0, "ACGT", nil, nil)
	rec3 := encodeBamRecord("primary2", 0, 0, 200, 50, cigar, -1, -1, 0, "ACGT", nil, nil)

	data, err := writeBamTestFile(headerText, refs, [][]byte{rec1, rec2, rec3})
	if err != nil {
		t.Fatalf("writeBamTestFile: %v", err)
	}

	// Filter out secondary alignments (0x100)
	opts := NewSamReaderOpts().FlagFilter(0x100)
	reader, err := NewBamReader(io.NopCloser(bytes.NewReader(data)), opts)
	if err != nil {
		t.Fatalf("NewBamReader: %v", err)
	}
	defer reader.Close()

	var names []string
	for {
		r, err := reader.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Next: %v", err)
		}
		names = append(names, r.ReadName)
	}

	if len(names) != 2 || names[0] != "primary" || names[1] != "primary2" {
		t.Errorf("expected [primary, primary2], got %v", names)
	}
}

func TestBamReaderUnmapped(t *testing.T) {
	refs := []bamRefInfo{{name: "chr1", length: 100000}}
	headerText := "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n"

	// Unmapped read: refID=-1, pos=-1, flag=4, cigar empty
	rec := encodeBamRecord("unmapped", 4, -1, -1, 0, nil, -1, -1, 0, "ACGT", nil, nil)

	data, err := writeBamTestFile(headerText, refs, [][]byte{rec})
	if err != nil {
		t.Fatalf("writeBamTestFile: %v", err)
	}

	reader, err := NewBamReader(io.NopCloser(bytes.NewReader(data)))
	if err != nil {
		t.Fatalf("NewBamReader: %v", err)
	}
	defer reader.Close()

	r, err := reader.Next()
	if err != nil {
		t.Fatalf("Next: %v", err)
	}
	if r.RefName != "*" {
		t.Errorf("RefName: got %q, want *", r.RefName)
	}
	if r.Cigar != "*" {
		t.Errorf("Cigar: got %q, want *", r.Cigar)
	}
	if !r.IsUnmapped() {
		t.Error("expected unmapped flag")
	}
}
