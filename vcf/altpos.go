package vcf

import (
	"strconv"
	"strings"
)

// VarType classifies a variant or structural-variant allele.
type VarType int

// Variant types, mirroring ngsutilsj VCFRecord.VCFVarType.
const (
	VarSNV VarType = iota
	VarBND
	VarDEL
	VarINS
	VarINV
	VarDUP
	VarCNV
	VarUNK
)

// String returns the canonical token for the variant type.
func (t VarType) String() string {
	switch t {
	case VarSNV:
		return "SNV"
	case VarBND:
		return "BND"
	case VarDEL:
		return "DEL"
	case VarINS:
		return "INS"
	case VarINV:
		return "INV"
	case VarDUP:
		return "DUP"
	case VarCNV:
		return "CNV"
	default:
		return "UNK"
	}
}

// ParseVarType maps an SVTYPE token to a VarType. TRA maps to BND.
func ParseVarType(s string) VarType {
	switch s {
	case "DEL":
		return VarDEL
	case "BND", "TRA":
		return VarBND
	case "INS":
		return VarINS
	case "INV":
		return VarINV
	case "DUP":
		return VarDUP
	case "CNV":
		return VarCNV
	default:
		return VarUNK
	}
}

// SVConnection describes the breakend orientation of a structural variant.
type SVConnection int

// Breakend connection types, mirroring ngsutilsj VCFRecord.VCFSVConnection.
const (
	ConnNA SVConnection = iota
	Conn5to5
	Conn5to3
	Conn3to3
	Conn3to5
	ConnNtoN
	ConnUNK
)

// String returns the canonical token for the connection (empty for NA).
func (c SVConnection) String() string {
	switch c {
	case Conn5to5:
		return "5to5"
	case Conn5to3:
		return "5to3"
	case Conn3to3:
		return "3to3"
	case Conn3to5:
		return "3to5"
	case ConnNtoN:
		return "NToN"
	case ConnNA:
		return ""
	default:
		return "UNK"
	}
}

// ParseSVConnection maps a CT-style token to an SVConnection.
func ParseSVConnection(s string) SVConnection {
	switch strings.ToLower(s) {
	case "5to5":
		return Conn5to5
	case "5to3":
		return Conn5to3
	case "3to3":
		return Conn3to3
	case "3to5":
		return Conn3to5
	case "nton":
		return ConnNtoN
	default:
		return ConnNA
	}
}

// AltPos is a resolved alternate breakpoint for an ALT allele.
type AltPos struct {
	Chrom    string
	Pos      int // 1-based
	Type     VarType
	ConnType SVConnection
	Alt      string
	Extra    string
}

// AltPositions resolves the alternate position(s) for the record's ALT
// allele(s), interpreting structural-variant notations (BND brackets, symbolic
// <DEL>/<DUP>/<INV>/<CNV>/<INS>, and in-place indels). It ports
// VCFRecord.getAltPos.
//
// altChromKey/altPosKey/typeKey/connKey are optional INFO field names; pass ""
// to use the defaults (altPosKey defaults to END, typeKey to SVTYPE). When
// altChromKey is non-empty, a single AltPos is returned from the named INFO
// fields.
func (r *VcfRecord) AltPositions(altChromKey, altPosKey, typeKey, connKey string) []AltPos {
	if altPosKey == "" {
		altPosKey = "END"
	}
	if typeKey == "" {
		typeKey = "SVTYPE"
	}

	var ret []AltPos

	if altChromKey != "" {
		info := r.Info()
		ac, okC := info.Get(altChromKey)
		ap, okP := info.Get(altPosKey)
		if !okC || !okP {
			return ret
		}
		altPos, err := ap.Int()
		if err != nil {
			return ret
		}
		altType := VarSNV
		altConn := ConnNA
		if v, ok := info.Get(typeKey); ok {
			altType = ParseVarType(v.String())
		}
		if connKey != "" {
			if v, ok := info.Get(connKey); ok {
				altConn = ParseSVConnection(v.String())
			}
		}
		ret = append(ret, AltPos{Chrom: ac.String(), Pos: altPos, Type: altType, ConnType: altConn, Alt: r.AltOrig()})
		return ret
	}

	for _, alt := range r.Alt() {
		altType := VarSNV
		altConn := ConnNA
		altChrom := r.Chrom
		altPos := r.Pos
		extra := ""

		switch {
		case strings.HasPrefix(alt, "[") || strings.HasPrefix(alt, "]") ||
			strings.HasSuffix(alt, "[") || strings.HasSuffix(alt, "]"):
			// BND
			altType = VarBND
			var sub string
			switch {
			case strings.HasPrefix(alt, "["):
				altConn = Conn5to5
				idx := strings.Index(alt[1:], "[") + 1
				sub = alt[1:idx]
				extra = alt[idx+1:]
			case strings.HasPrefix(alt, "]"):
				altConn = Conn5to3
				idx := strings.Index(alt[1:], "]") + 1
				sub = alt[1:idx]
				extra = alt[idx+1:]
			case strings.HasSuffix(alt, "["):
				altConn = Conn3to5
				idx := strings.Index(alt, "[")
				sub = alt[idx+1 : len(alt)-1]
				extra = alt[:strings.Index(alt[1:], "[")+1]
			case strings.HasSuffix(alt, "]"):
				altConn = Conn3to3
				idx := strings.Index(alt, "]")
				sub = alt[idx+1 : len(alt)-1]
				extra = alt[:strings.Index(alt[1:], "]")+1]
			}
			spl := strings.Split(sub, ":")
			if len(spl) == 2 {
				altChrom = spl[0]
				if p, err := strconv.Atoi(spl[1]); err == nil {
					altPos = p
				}
			}

		case strings.HasPrefix(alt, "<") && !strings.HasPrefix(alt, "<INS"):
			// DEL, DUP, INV, CNV
			altChrom = r.Chrom
			switch {
			case strings.HasPrefix(alt, "<CNV"):
				altType = VarCNV
				altConn = ConnNA
			case strings.HasPrefix(alt, "<INV"):
				altType = VarINV
				altConn = ConnUNK
			case strings.HasPrefix(alt, "<DEL"):
				altType = VarDEL
				altConn = Conn3to5
			case strings.HasPrefix(alt, "<DUP"):
				altType = VarDUP
				altConn = Conn5to3
			}
			// These must specify an END value.
			if v, ok := r.Info().Get(altPosKey); ok {
				if p, err := v.Int(); err == nil {
					altPos = p
				}
			}

		default:
			// INS or SNV
			switch {
			case strings.HasPrefix(alt, "<INS") || len(alt) > 1:
				altChrom = r.Chrom
				altPos = r.Pos
				altType = VarINS
				altConn = ConnNtoN
			case len(r.Ref) > 1:
				// in-place deletion; VCF pos is 1-based
				altChrom = r.Chrom
				altPos = r.Pos + len(r.Ref)
				altType = VarDEL
				altConn = Conn3to5
			default:
				altChrom = r.Chrom
				altPos = r.Pos
				altType = VarSNV
				altConn = ConnNA
			}
		}

		info := r.Info()
		if v, ok := info.Get(altPosKey); ok {
			if p, err := v.Int(); err == nil {
				altPos = p
			}
		}
		if v, ok := info.Get(typeKey); ok {
			altType = ParseVarType(v.String())
		}
		if connKey != "" {
			if v, ok := info.Get(connKey); ok {
				altConn = ParseSVConnection(v.String())
			}
		} else if altType == VarINV && altConn == ConnUNK {
			if v, ok := info.Get("CT"); ok {
				altConn = ParseSVConnection(v.String())
			}
		}

		ret = append(ret, AltPos{Chrom: altChrom, Pos: altPos, Type: altType, ConnType: altConn, Alt: alt, Extra: extra})
	}

	return ret
}
