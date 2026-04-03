package htsio

import (
	"testing"
)

func TestParseSamLine(t *testing.T) {
	line := "read1\t0\tchr1\t100\t60\t50M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\t*\tNM:i:0\tMD:Z:50"
	rec, err := parseSamLine(line)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	if rec.QName != "read1" {
		t.Errorf("QName = %q, want %q", rec.QName, "read1")
	}
	if rec.Flag != 0 {
		t.Errorf("Flag = %d, want 0", rec.Flag)
	}
	if rec.RName != "chr1" {
		t.Errorf("RName = %q, want %q", rec.RName, "chr1")
	}
	if rec.Pos != 100 {
		t.Errorf("Pos = %d, want 100", rec.Pos)
	}
	if rec.MapQ != 60 {
		t.Errorf("MapQ = %d, want 60", rec.MapQ)
	}
	if rec.Cigar != "50M" {
		t.Errorf("Cigar = %q, want %q", rec.Cigar, "50M")
	}
	if len(rec.Tags) != 2 {
		t.Errorf("len(Tags) = %d, want 2", len(rec.Tags))
	}
	if nm, ok := rec.Tags["NM"]; !ok {
		t.Error("missing NM tag")
	} else {
		if nm.Type != 'i' {
			t.Errorf("NM.Type = %c, want 'i'", nm.Type)
		}
		if v, ok := nm.IntValue(); !ok || v != 0 {
			t.Errorf("NM.IntValue() = %d, %v; want 0, true", v, ok)
		}
	}
	if md, ok := rec.Tags["MD"]; !ok {
		t.Error("missing MD tag")
	} else {
		if md.Type != 'Z' {
			t.Errorf("MD.Type = %c, want 'Z'", md.Type)
		}
		if md.Value != "50" {
			t.Errorf("MD.Value = %q, want %q", md.Value, "50")
		}
	}
}

func TestParseSamLineMinimalFields(t *testing.T) {
	line := "read2\t16\tchr2\t200\t30\t30M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTAC\t*"
	rec, err := parseSamLine(line)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	if rec.QName != "read2" {
		t.Errorf("QName = %q, want %q", rec.QName, "read2")
	}
	if rec.Flag != 16 {
		t.Errorf("Flag = %d, want 16", rec.Flag)
	}
	if !rec.IsReverse() {
		t.Error("expected IsReverse() = true")
	}
	if len(rec.Tags) != 0 {
		t.Errorf("len(Tags) = %d, want 0", len(rec.Tags))
	}
}

func TestParseSamLineTooFewFields(t *testing.T) {
	line := "read1\t0\tchr1"
	_, err := parseSamLine(line)
	if err == nil {
		t.Error("expected error for too few fields")
	}
}

func TestSamRecordFlags(t *testing.T) {
	rec := &SamRecord{Flag: 0x4 | 0x10 | 0x100 | 0x400 | 0x800}
	if !rec.IsUnmapped() {
		t.Error("expected IsUnmapped()")
	}
	if !rec.IsReverse() {
		t.Error("expected IsReverse()")
	}
	if !rec.IsSecondary() {
		t.Error("expected IsSecondary()")
	}
	if !rec.IsDuplicate() {
		t.Error("expected IsDuplicate()")
	}
	if !rec.IsSupplementary() {
		t.Error("expected IsSupplementary()")
	}

	rec2 := &SamRecord{Flag: 0}
	if rec2.IsUnmapped() {
		t.Error("expected !IsUnmapped()")
	}
}

func TestSamReaderBuilder(t *testing.T) {
	r := NewSamReader("test.bam").
		Region("chr1:100-200").
		FlagRequired(0x2).
		FlagFilter(0x4 | 0x100).
		MinMapQ(20)

	if r.filename != "test.bam" {
		t.Errorf("filename = %q, want %q", r.filename, "test.bam")
	}
	if r.region != "chr1:100-200" {
		t.Errorf("region = %q, want %q", r.region, "chr1:100-200")
	}
	if r.flagReq != 0x2 {
		t.Errorf("flagReq = %d, want %d", r.flagReq, 0x2)
	}
	if r.flagFilter != 0x104 {
		t.Errorf("flagFilter = %d, want %d", r.flagFilter, 0x104)
	}
	if r.minMapQ != 20 {
		t.Errorf("minMapQ = %d, want %d", r.minMapQ, 20)
	}

	// verify it satisfies the SamReader interface
	var _ SamReader = r
}

func TestSamHeader(t *testing.T) {
	h := NewSamHeader()
	h.AddLine("@HD\tVN:1.6\tSO:coordinate")
	h.AddLine("@SQ\tSN:chr1\tLN:248956422")
	h.AddLine("@SQ\tSN:chr2\tLN:242193529")
	h.AddLine("@RG\tID:sample1\tSM:sample1")

	refs := h.References()
	if len(refs) != 2 {
		t.Fatalf("len(refs) = %d, want 2", len(refs))
	}
	if refs[0].Name != "chr1" || refs[0].Length != 248956422 {
		t.Errorf("refs[0] = %+v, want chr1/248956422", refs[0])
	}
	if refs[1].Name != "chr2" || refs[1].Length != 242193529 {
		t.Errorf("refs[1] = %+v, want chr2/242193529", refs[1])
	}

	rgs := h.ReadGroups()
	if len(rgs) != 1 || rgs[0] != "sample1" {
		t.Errorf("ReadGroups() = %v, want [sample1]", rgs)
	}

	text := h.Text()
	if !contains(text, "@HD") || !contains(text, "@SQ") || !contains(text, "@RG") {
		t.Errorf("Text() missing header lines")
	}
}

func TestSamRecordString(t *testing.T) {
	rec := &SamRecord{
		QName: "read1",
		Flag:  0,
		RName: "chr1",
		Pos:   100,
		MapQ:  60,
		Cigar: "50M",
		RNext: "*",
		PNext: 0,
		TLen:  0,
		Seq:   "ACGT",
		Qual:  "IIII",
		Tags:  map[string]SamTag{},
	}

	s := rec.String()
	if s != "read1\t0\tchr1\t100\t60\t50M\t*\t0\t0\tACGT\tIIII" {
		t.Errorf("String() = %q", s)
	}
}

func TestSamWriterBuilder(t *testing.T) {
	h := NewSamHeader()
	w := NewSamWriter("out.bam", h).
		Format(FormatCRAM).
		Reference("ref.fa")

	if w.filename != "out.bam" {
		t.Errorf("filename = %q, want %q", w.filename, "out.bam")
	}
	if w.format != FormatCRAM {
		t.Errorf("format = %d, want %d", w.format, FormatCRAM)
	}
	if w.refFile != "ref.fa" {
		t.Errorf("refFile = %q, want %q", w.refFile, "ref.fa")
	}

	// verify it satisfies the SamWriter interface
	var _ SamWriter = w
}

func TestStdoutSamWriterInterface(t *testing.T) {
	var _ SamWriter = NewStdoutSamWriter(nil)
}

func TestCigarRefLen(t *testing.T) {
	tests := []struct {
		cigar string
		want  int
	}{
		{"50M", 50},
		{"*", 0},
		{"10M2I10M", 20},
		{"10M3D10M", 23},
		{"5S10M5S", 10},
		{"10M1I5M2D3M", 20},
		{"100M", 100},
		{"10H20M10H", 20},
		{"5S10M2I3M1D5M5S", 19},
		{"10=5X", 15},
		{"10M100N10M", 120},
	}
	for _, tt := range tests {
		got := CigarRefLen(tt.cigar)
		if got != tt.want {
			t.Errorf("CigarRefLen(%q) = %d, want %d", tt.cigar, got, tt.want)
		}
	}
}

func contains(s, substr string) bool {
	return len(s) >= len(substr) && (s == substr || len(s) > 0 && containsHelper(s, substr))
}

func containsHelper(s, substr string) bool {
	for i := 0; i <= len(s)-len(substr); i++ {
		if s[i:i+len(substr)] == substr {
			return true
		}
	}
	return false
}
