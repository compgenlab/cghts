package htsio

import (
	"testing"
)

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

	opts := NewSamReaderOpts().
		FlagRequired(0x2).
		FlagFilter(0x4 | 0x100).
		MinMapQ(20)

	if opts.FlagReqValue() != 0x2 {
		t.Errorf("flagReq = %d, want %d", opts.FlagReqValue(), 0x2)
	}
	if opts.FlagFilterValue() != 0x104 {
		t.Errorf("flagFilter = %d, want %d", opts.FlagFilterValue(), 0x104)
	}
	if opts.MinMapQValue() != 20 {
		t.Errorf("minMapQ = %d, want %d", opts.MinMapQValue(), 20)
	}

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
		ReadName:  "read1",
		Flag:      0,
		RefName:   "chr1",
		Pos:       100,
		MapQ:      60,
		Cigar:     "50M",
		RefNext:   "*",
		PosNext:   0,
		InsertLen: 0,
		Seq:       "ACGT",
		Qual:      "IIII",
		Tags:      map[string]SamTag{},
	}

	s := rec.String()
	if s != "read1\t0\tchr1\t100\t60\t50M\t*\t0\t0\tACGT\tIIII" {
		t.Errorf("String() = %q", s)
	}
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

func TestParseRegion(t *testing.T) {
	tests := []struct {
		input string
		ref   string
		start int
		end   int
	}{
		{"chr1", "chr1", 0, -1},
		{"chr1:1000-2000", "chr1", 999, 2000},
		{"chr1:1000", "chr1", 999, -1},
		{"chr2:1,000-2,000", "chr2", 999, 2000},
	}

	for _, tt := range tests {
		ref, start, end, err := ParseRegion(tt.input)
		if err != nil {
			t.Errorf("ParseRegion(%q): %v", tt.input, err)
			continue
		}
		if ref != tt.ref || start != tt.start || end != tt.end {
			t.Errorf("ParseRegion(%q) = (%q, %d, %d), want (%q, %d, %d)",
				tt.input, ref, start, end, tt.ref, tt.start, tt.end)
		}
	}
}

