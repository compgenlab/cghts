package vcf

import "testing"

// The VarType/SVConnection tokens are read off INFO/SVTYPE and INFO/CT, and
// written back out by cgkit's vcf-tobedpe. A parser and a printer that disagree
// would produce a file that does not round-trip through this package's own reader.

func TestParseVarType(t *testing.T) {
	for _, tc := range []struct {
		in   string
		want VarType
	}{
		{"DEL", VarDEL},
		{"INS", VarINS},
		{"INV", VarINV},
		{"DUP", VarDUP},
		{"CNV", VarCNV},
		{"BND", VarBND},
		// TRA is the pre-4.x spelling for a translocation, which VCF 4.x models
		// as a breakend pair. Callers still emit it, so it aliases rather than
		// falling through to UNK.
		{"TRA", VarBND},
		{"", VarUNK},
		{"del", VarUNK}, // SVTYPE tokens are upper-case; no folding
		{"DEL:ME:ALU", VarUNK},
		{"nonsense", VarUNK},
	} {
		if got := ParseVarType(tc.in); got != tc.want {
			t.Errorf("ParseVarType(%q) = %v, want %v", tc.in, got, tc.want)
		}
	}
}

// Every type prints as a token that parses back to itself -- except VarSNV, which
// is the zero value used as "no SVTYPE said otherwise" rather than a token any
// caller writes. AltPositions never emits "SNV" into SVTYPE, so the asymmetry is
// unreachable; pinned here so a future round-trip assumption fails loudly.
func TestVarTypeRoundTrip(t *testing.T) {
	for _, vt := range []VarType{VarBND, VarDEL, VarINS, VarINV, VarDUP, VarCNV, VarUNK} {
		if got := ParseVarType(vt.String()); got != vt {
			t.Errorf("ParseVarType(%q) = %v, want %v", vt.String(), got, vt)
		}
	}
	if VarSNV.String() != "SNV" {
		t.Errorf("VarSNV.String() = %q, want SNV", VarSNV.String())
	}
	if got := ParseVarType("SNV"); got != VarUNK {
		t.Errorf("ParseVarType(SNV) = %v, want VarUNK -- SNV is not an SVTYPE token", got)
	}
}

func TestParseSVConnection(t *testing.T) {
	for _, tc := range []struct {
		in   string
		want SVConnection
	}{
		{"5to5", Conn5to5},
		{"5to3", Conn5to3},
		{"3to3", Conn3to3},
		{"3to5", Conn3to5},
		{"NtoN", ConnNtoN},
		// Unlike SVTYPE, CT tokens are matched case-insensitively -- callers
		// disagree on the capitalization of NtoN.
		{"NTON", ConnNtoN},
		{"nton", ConnNtoN},
		{"5TO5", Conn5to5},
		{"", ConnNA},
		{"nonsense", ConnNA},
	} {
		if got := ParseSVConnection(tc.in); got != tc.want {
			t.Errorf("ParseSVConnection(%q) = %v, want %v", tc.in, got, tc.want)
		}
	}
}

func TestSVConnectionRoundTrip(t *testing.T) {
	for _, c := range []SVConnection{Conn5to5, Conn5to3, Conn3to3, Conn3to5, ConnNtoN, ConnNA} {
		if got := ParseSVConnection(c.String()); got != c {
			t.Errorf("ParseSVConnection(%q) = %v, want %v", c.String(), got, c)
		}
	}
	// ConnUNK is the "symbolic <INV> with no CT yet" placeholder AltPositions
	// sets internally, not a token. It prints as UNK and parses back to NA,
	// which is correct -- an unrecognized CT is an absence, not a distinct
	// orientation.
	if ConnUNK.String() != "UNK" {
		t.Errorf("ConnUNK.String() = %q, want UNK", ConnUNK.String())
	}
	if got := ParseSVConnection("UNK"); got != ConnNA {
		t.Errorf("ParseSVConnection(UNK) = %v, want ConnNA", got)
	}
	// NA prints as empty, so a BEDPE writer emits a blank column rather than
	// the literal "NA".
	if ConnNA.String() != "" {
		t.Errorf("ConnNA.String() = %q, want empty", ConnNA.String())
	}
}
