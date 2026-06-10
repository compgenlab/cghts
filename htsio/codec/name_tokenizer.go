package codec

import (
	"encoding/binary"
	"fmt"
	"strconv"
)

// Name tokenizer (tok3) decoder for CRAM v3.1.
// Decompresses read names that have been tokenized and entropy-coded.
// Based on the htscodecs tokenise_name3.c reference implementation.

// Token types — must match the enum order in htscodecs tokenise_name3.c.
const (
	tokType    = 0  // token type descriptor
	tokAlpha   = 1  // alphabetic string
	tokChar    = 2  // single character
	tokDigits0 = 3  // zero-padded integer
	tokDZLen   = 4  // zero-pad length for DIGITS0/DDELTA0
	tokDup     = 5  // duplicate entire previous name
	tokDiff    = 6  // different from previous
	tokDigits  = 7  // variable-width integer
	tokDDelta  = 8  // delta of variable-width integer
	tokDDelta0 = 9  // delta of zero-padded integer
	tokMatch   = 10 // copy token from previous name
	tokNop     = 11 // no-op
	tokEnd     = 12 // end of name
)

const maxTokens = 128

// lastToken stores the value of a token from the previous name.
type lastToken struct {
	tokenType byte
	intVal    uint32
	strOffset int // offset into name string
	strLen    int // length of string token
}

// nameContext holds state for decoding a sequence of names.
type nameContext struct {
	desc    []descriptor // (token << 4) | type → data
	maxTok  int
	counter int

	// Per-name state from previous names.
	lastNames  []string
	lastTokens [][]lastToken
	lastNTok   []int
}

type descriptor struct {
	data []byte
	pos  int
}

func DecodeNameTokenizer(data []byte) ([]byte, error) {
	if len(data) < 9 {
		return nil, fmt.Errorf("tok3: data too short")
	}

	ulen := binary.LittleEndian.Uint32(data[0:4])
	nreads := int(binary.LittleEndian.Uint32(data[4:8]))
	useArith := data[8]

	if useArith != 0 {
		return nil, fmt.Errorf("tok3: arithmetic coding not supported (only rANS)")
	}

	// nreads is attacker-controlled and is used to size several per-read slices
	// below. Each read must contribute at least one byte to the token streams
	// (every read ends with an explicit end token), so nreads can never exceed
	// the encoded length; reject anything larger before allocating to avoid a
	// huge make() on a crafted header.
	if nreads < 0 || nreads > len(data) {
		return nil, fmt.Errorf("tok3: implausible read count %d for %d bytes of input", nreads, len(data))
	}

	ctx := &nameContext{
		desc:       make([]descriptor, maxTokens*16),
		lastNames:  make([]string, nreads),
		lastTokens: make([][]lastToken, nreads),
		lastNTok:   make([]int, nreads),
	}

	// Unpack descriptor blocks.
	o := 9
	tnum := -1

	for o < len(data) {
		ttype := data[o]
		o++

		if ttype&64 != 0 {
			// Copy descriptor from another token.
			if o+2 > len(data) {
				return nil, fmt.Errorf("tok3: truncated copy descriptor")
			}
			srcIdx := int(data[o])<<4 + int(data[o+1])
			o += 2

			if ttype&128 != 0 {
				tnum++
				if tnum >= maxTokens {
					return nil, fmt.Errorf("tok3: too many tokens")
				}
				ctx.maxTok = tnum + 1
				// Initialize type descriptor with defaults.
				if ttype&15 != 0 {
					ctx.initTypeDesc(tnum, ttype&15, nreads)
				}
			}
			if tnum < 0 {
				return nil, fmt.Errorf("tok3: copy before first token")
			}

			dstIdx := (tnum << 4) | int(ttype&15)
			if srcIdx >= dstIdx || srcIdx >= len(ctx.desc) {
				return nil, fmt.Errorf("tok3: invalid copy source %d", srcIdx)
			}
			// Copy the descriptor data.
			src := ctx.desc[srcIdx]
			ctx.desc[dstIdx] = descriptor{
				data: make([]byte, len(src.data)),
				pos:  0,
			}
			copy(ctx.desc[dstIdx].data, src.data)
			continue
		}

		if ttype&128 != 0 {
			tnum++
			if tnum >= maxTokens {
				return nil, fmt.Errorf("tok3: too many tokens")
			}
			ctx.maxTok = tnum + 1
			// Initialize type descriptor with MATCH defaults.
			if ttype&15 != 0 {
				ctx.initTypeDesc(tnum, ttype&15, nreads)
			}
		}

		if tnum < 0 {
			return nil, fmt.Errorf("tok3: data before first token")
		}

		idx := (tnum << 4) | int(ttype&15)
		if idx >= len(ctx.desc) {
			return nil, fmt.Errorf("tok3: descriptor index %d out of range", idx)
		}

		// Read compressed data for this descriptor.
		// Format: varint(clen) + clen bytes of rANS Nx16 data
		if o >= len(data) {
			return nil, fmt.Errorf("tok3: truncated descriptor data")
		}
		n1, clen := varGetU32(data[o:])
		if n1 == 0 {
			return nil, fmt.Errorf("tok3: truncated clen")
		}
		o += n1

		if o+int(clen) > len(data) {
			return nil, fmt.Errorf("tok3: descriptor data overflows (clen=%d, remaining=%d)", clen, len(data)-o)
		}

		// Decompress using rANS Nx16.
		compData := data[o : o+int(clen)]
		decompressed, err := DecodeRansNx16(compData)
		if err != nil {
			return nil, fmt.Errorf("tok3: decompressing descriptor %d: %w", idx, err)
		}
		ctx.desc[idx] = descriptor{data: decompressed, pos: 0}
		o += int(clen)
	}

	// Reconstruct names.
	out := make([]byte, 0, ulen)
	for i := 0; i < nreads; i++ {
		name, err := ctx.decodeName(i)
		if err != nil {
			return nil, fmt.Errorf("tok3: decoding name %d: %w", i, err)
		}
		out = append(out, name...)
		out = append(out, 0) // null terminator
	}

	return out, nil
}

func (ctx *nameContext) initTypeDesc(tnum int, ttype byte, nreads int) {
	idx := tnum << 4
	buf := make([]byte, nreads)
	buf[0] = ttype
	for i := 1; i < nreads; i++ {
		buf[i] = tokMatch
	}
	ctx.desc[idx] = descriptor{data: buf, pos: 0}
}

func (ctx *nameContext) readByte(idx int) (byte, error) {
	d := &ctx.desc[idx]
	if d.pos >= len(d.data) {
		return 0, fmt.Errorf("descriptor %d exhausted", idx)
	}
	b := d.data[d.pos]
	d.pos++
	return b, nil
}

func (ctx *nameContext) readUint32(idx int) (uint32, error) {
	d := &ctx.desc[idx]
	if d.pos+4 > len(d.data) {
		return 0, fmt.Errorf("descriptor %d: not enough data for uint32", idx)
	}
	v := binary.LittleEndian.Uint32(d.data[d.pos:])
	d.pos += 4
	return v, nil
}

func (ctx *nameContext) readString(idx int) (string, error) {
	d := &ctx.desc[idx]
	start := d.pos
	for d.pos < len(d.data) {
		if d.data[d.pos] == 0 {
			s := string(d.data[start:d.pos])
			d.pos++ // skip null
			return s, nil
		}
		d.pos++
	}
	// No null terminator — take the rest.
	s := string(d.data[start:d.pos])
	return s, nil
}

func (ctx *nameContext) decodeName(cnum int) (string, error) {
	// Read token 0: type
	t0, err := ctx.readByte(0 << 4)
	if err != nil {
		return "", fmt.Errorf("reading token 0 type: %w", err)
	}

	// Read distance to reference name.
	dist := uint32(0)
	if t0 != tokEnd {
		dist, err = ctx.readUint32((0 << 4) | int(t0))
		if err != nil {
			return "", fmt.Errorf("reading token 0 dist: %w", err)
		}
	}
	pnum := cnum - int(dist)
	if pnum < 0 {
		pnum = 0
	}

	// Handle DUP: entire name is a copy.
	if t0 == tokDup {
		if pnum == cnum {
			return "", fmt.Errorf("DUP self-reference")
		}
		name := ctx.lastNames[pnum]
		ctx.lastNames[cnum] = name
		ctx.lastNTok[cnum] = ctx.lastNTok[pnum]
		ctx.lastTokens[cnum] = make([]lastToken, len(ctx.lastTokens[pnum]))
		copy(ctx.lastTokens[cnum], ctx.lastTokens[pnum])
		return name, nil
	}

	// Decode token by token.
	var name []byte
	tokens := make([]lastToken, maxTokens)

	for ntok := 1; ntok < maxTokens && ntok < ctx.maxTok; ntok++ {
		tok, err := ctx.readByte(ntok << 4)
		if err != nil {
			return "", fmt.Errorf("reading token %d type: %w", ntok, err)
		}

		switch tok {
		case tokChar:
			c, err := ctx.readByte((ntok << 4) | tokChar)
			if err != nil {
				return "", err
			}
			name = append(name, c)
			tokens[ntok] = lastToken{tokenType: tokChar, intVal: uint32(c)}

		case tokAlpha:
			s, err := ctx.readString((ntok << 4) | tokAlpha)
			if err != nil {
				return "", err
			}
			off := len(name)
			name = append(name, s...)
			tokens[ntok] = lastToken{tokenType: tokAlpha, strOffset: off, strLen: len(s)}

		case tokDigits0:
			vl, err := ctx.readByte((ntok << 4) | tokDZLen)
			if err != nil {
				return "", err
			}
			v, err := ctx.readUint32((ntok << 4) | tokDigits0)
			if err != nil {
				return "", err
			}
			s := appendUint32Fixed(v, int(vl))
			name = append(name, s...)
			tokens[ntok] = lastToken{tokenType: tokDigits0, intVal: v, strLen: int(vl)}

		case tokDigits:
			v, err := ctx.readUint32((ntok << 4) | tokDigits)
			if err != nil {
				return "", err
			}
			s := strconv.FormatUint(uint64(v), 10)
			name = append(name, s...)
			tokens[ntok] = lastToken{tokenType: tokDigits, intVal: v}

		case tokDDelta0:
			if pnum < 0 || pnum >= cnum || ntok >= ctx.lastNTok[pnum] {
				return "", fmt.Errorf("DDELTA0 invalid reference")
			}
			d, err := ctx.readByte((ntok << 4) | tokDDelta0)
			if err != nil {
				return "", err
			}
			v := ctx.lastTokens[pnum][ntok].intVal + uint32(d)
			padLen := ctx.lastTokens[pnum][ntok].strLen
			s := appendUint32Fixed(v, padLen)
			name = append(name, s...)
			tokens[ntok] = lastToken{tokenType: tokDigits0, intVal: v, strLen: padLen}

		case tokDDelta:
			if pnum < 0 || pnum >= cnum || ntok >= ctx.lastNTok[pnum] {
				return "", fmt.Errorf("DDELTA invalid reference")
			}
			d, err := ctx.readByte((ntok << 4) | tokDDelta)
			if err != nil {
				return "", err
			}
			v := ctx.lastTokens[pnum][ntok].intVal + uint32(d)
			s := strconv.FormatUint(uint64(v), 10)
			name = append(name, s...)
			tokens[ntok] = lastToken{tokenType: tokDigits, intVal: v}

		case tokMatch:
			if pnum < 0 || pnum >= cnum || ntok >= ctx.lastNTok[pnum] {
				return "", fmt.Errorf("MATCH invalid reference (pnum=%d, cnum=%d, ntok=%d, lastNTok=%d)", pnum, cnum, ntok, ctx.lastNTok[pnum])
			}
			prev := ctx.lastTokens[pnum][ntok]
			switch prev.tokenType {
			case tokChar:
				name = append(name, byte(prev.intVal))
				tokens[ntok] = lastToken{tokenType: tokChar, intVal: prev.intVal}
			case tokAlpha:
				prevName := ctx.lastNames[pnum]
				if prev.strOffset+prev.strLen <= len(prevName) {
					s := prevName[prev.strOffset : prev.strOffset+prev.strLen]
					off := len(name)
					name = append(name, s...)
					tokens[ntok] = lastToken{tokenType: tokAlpha, strOffset: off, strLen: prev.strLen}
				}
			case tokDigits:
				s := strconv.FormatUint(uint64(prev.intVal), 10)
				name = append(name, s...)
				tokens[ntok] = lastToken{tokenType: tokDigits, intVal: prev.intVal}
			case tokDigits0:
				s := appendUint32Fixed(prev.intVal, prev.strLen)
				name = append(name, s...)
				tokens[ntok] = lastToken{tokenType: tokDigits0, intVal: prev.intVal, strLen: prev.strLen}
			default:
				return "", fmt.Errorf("MATCH unsupported prev type %d", prev.tokenType)
			}

		case tokNop:
			tokens[ntok] = lastToken{tokenType: tokNop}

		case tokEnd:
			tokens[ntok] = lastToken{tokenType: tokEnd}
			ctx.lastNames[cnum] = string(name)
			ctx.lastNTok[cnum] = ntok
			ctx.lastTokens[cnum] = make([]lastToken, ntok+1)
			copy(ctx.lastTokens[cnum], tokens[:ntok+1])
			return string(name), nil

		default:
			return "", fmt.Errorf("unknown token type %d at token %d", tok, ntok)
		}
	}

	// Reached max tokens without END.
	ctx.lastNames[cnum] = string(name)
	ctx.lastNTok[cnum] = ctx.maxTok
	ctx.lastTokens[cnum] = make([]lastToken, ctx.maxTok)
	copy(ctx.lastTokens[cnum], tokens[:ctx.maxTok])
	return string(name), nil
}

func appendUint32Fixed(v uint32, width int) string {
	s := strconv.FormatUint(uint64(v), 10)
	for len(s) < width {
		s = "0" + s
	}
	return s
}
