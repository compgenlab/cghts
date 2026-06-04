package codec

import (
	"bytes"
	"encoding/binary"
	"strconv"
	"strings"
)

// EncodeNameTokenizer compresses read names using the tok3 name tokenizer codec.
// names should be the read names (without null terminators).
// Returns the compressed data in the format expected by DecodeNameTokenizer.
func EncodeNameTokenizer(names []string) []byte {
	if len(names) == 0 {
		return nil
	}

	nreads := len(names)

	// Phase 1: Tokenize all names.
	allTokens := make([][]nameToken, nreads)
	for i, name := range names {
		allTokens[i] = tokenizeName(name)
	}

	// Determine max token count.
	maxTok := 0
	for _, toks := range allTokens {
		if len(toks) > maxTok {
			maxTok = len(toks)
		}
	}
	if maxTok > maxTokens-1 {
		maxTok = maxTokens - 1
	}

	// Phase 2: Build descriptor streams by comparing with previous name.
	descs := make(map[int]*bytes.Buffer)
	getDesc := func(idx int) *bytes.Buffer {
		if b, ok := descs[idx]; ok {
			return b
		}
		b := &bytes.Buffer{}
		descs[idx] = b
		return b
	}
	writeByte := func(idx int, v byte) { getDesc(idx).WriteByte(v) }
	writeUint32 := func(idx int, v uint32) {
		var b [4]byte
		binary.LittleEndian.PutUint32(b[:], v)
		getDesc(idx).Write(b[:])
	}
	writeString := func(idx int, s string) {
		getDesc(idx).WriteString(s)
		getDesc(idx).WriteByte(0) // null terminator
	}

	// Track previous name's tokens for MATCH/DDELTA comparison.
	var prevTokens []nameToken

	for i := 0; i < nreads; i++ {
		tokens := allTokens[i]

		// Token 0: type/distance header.
		if i > 0 && names[i] == names[i-1] {
			// DUP: entire name is identical to previous.
			writeByte(0<<4, tokDup)
			writeUint32((0<<4)|tokDup, 1) // distance = 1
			// prevTokens stays the same.
			continue
		}

		// DIFF from previous (distance = 1 for simplicity).
		writeByte(0<<4, tokDiff)
		writeUint32((0<<4)|tokDiff, 1)

		// Tokens 1..N.
		for ntok := 0; ntok < len(tokens); ntok++ {
			tok := tokens[ntok]
			descTok := ntok + 1 // token 0 is type/dist, real tokens start at 1

			if i > 0 && ntok < len(prevTokens) && tokensMatch(tok, prevTokens[ntok]) {
				// MATCH: identical to previous name's token.
				writeByte(descTok<<4, tokMatch)
			} else if i > 0 && ntok < len(prevTokens) && canDeltaToken(tok, prevTokens[ntok]) {
				// DDELTA: numeric token with small delta.
				delta := byte(tok.intVal - prevTokens[ntok].intVal)
				if tok.tokType == tokDigits0 {
					writeByte(descTok<<4, tokDDelta0)
					writeByte((descTok<<4)|tokDDelta0, delta)
				} else {
					writeByte(descTok<<4, tokDDelta)
					writeByte((descTok<<4)|tokDDelta, delta)
				}
			} else {
				// Literal token.
				switch tok.tokType {
				case tokChar:
					writeByte(descTok<<4, tokChar)
					writeByte((descTok<<4)|tokChar, tok.charVal)
				case tokAlpha:
					writeByte(descTok<<4, tokAlpha)
					writeString((descTok<<4)|tokAlpha, tok.strVal)
				case tokDigits:
					writeByte(descTok<<4, tokDigits)
					writeUint32((descTok<<4)|tokDigits, tok.intVal)
				case tokDigits0:
					writeByte(descTok<<4, tokDigits0)
					writeByte((descTok<<4)|tokDZLen, byte(tok.padLen))
					writeUint32((descTok<<4)|tokDigits0, tok.intVal)
				}
			}
		}

		// END token.
		endTok := len(tokens) + 1
		writeByte(endTok<<4, tokEnd)

		prevTokens = tokens
	}

	// Track which descriptor indices are actually used.
	maxDescTok := maxTok + 2 // +1 for token offset, +1 for END
	if maxDescTok > maxTokens {
		maxDescTok = maxTokens
	}

	// Phase 3: Compress each descriptor and pack into output.
	return packNameTokenizer(descs, maxDescTok, nreads, names)
}

// packNameTokenizer serializes the descriptor streams into the tok3 format.
func packNameTokenizer(descs map[int]*bytes.Buffer, maxDescTok, nreads int, names []string) []byte {
	// Compute uncompressed output size (names + null terminators).
	ulen := 0
	for _, name := range names {
		ulen += len(name) + 1
	}

	// Header.
	var out bytes.Buffer
	var hdr [9]byte
	binary.LittleEndian.PutUint32(hdr[0:4], uint32(ulen))
	binary.LittleEndian.PutUint32(hdr[4:8], uint32(nreads))
	hdr[8] = 0 // useArith = 0 (rANS Nx16)
	out.Write(hdr[:])

	// Write descriptor blocks.
	// For each token position, write NEW_TOKEN flag on first descriptor,
	// then each used type's compressed data.
	for tnum := 0; tnum < maxDescTok; tnum++ {
		firstForToken := true

		// Determine which types are used for this token.
		for ttype := 0; ttype < 16; ttype++ {
			idx := (tnum << 4) | ttype
			buf, ok := descs[idx]
			if !ok || buf.Len() == 0 {
				continue
			}

			data := buf.Bytes()

			// Compress with rANS Nx16.
			compressed := EncodeRansNx16(data)

			// Write ttype byte.
			ttypeByte := byte(ttype)
			if firstForToken {
				ttypeByte |= 0x80 // NEW_TOKEN flag
				firstForToken = false
			}
			out.WriteByte(ttypeByte)

			// Write compressed length + data.
			out.Write(varPutU32Slice(nil, uint32(len(compressed))))
			out.Write(compressed)
		}
	}

	return out.Bytes()
}

// nameToken is an intermediate token from tokenizing a read name.
type nameToken struct {
	tokType byte
	charVal byte   // for tokChar
	strVal  string // for tokAlpha
	intVal  uint32 // for tokDigits, tokDigits0
	padLen  int    // for tokDigits0 (width of zero-padded field)
}

// tokenizeName splits a read name into typed tokens.
func tokenizeName(name string) []nameToken {
	var tokens []nameToken
	i := 0
	for i < len(name) {
		c := name[i]

		if isAlpha(c) {
			// Run of alphabetic characters.
			j := i + 1
			for j < len(name) && isAlpha(name[j]) {
				j++
			}
			s := name[i:j]
			if len(s) == 1 {
				tokens = append(tokens, nameToken{tokType: tokChar, charVal: s[0]})
			} else {
				tokens = append(tokens, nameToken{tokType: tokAlpha, strVal: s})
			}
			i = j
		} else if isDigit(c) {
			// Run of digits.
			j := i + 1
			for j < len(name) && isDigit(name[j]) {
				j++
			}
			s := name[i:j]
			v, _ := strconv.ParseUint(s, 10, 32)
			if len(s) > 1 && s[0] == '0' {
				// Zero-padded.
				tokens = append(tokens, nameToken{
					tokType: tokDigits0,
					intVal:  uint32(v),
					padLen:  len(s),
				})
			} else {
				tokens = append(tokens, nameToken{
					tokType: tokDigits,
					intVal:  uint32(v),
				})
			}
			i = j
		} else {
			// Single character (separator, etc.).
			tokens = append(tokens, nameToken{tokType: tokChar, charVal: c})
			i++
		}
	}
	return tokens
}

func isAlpha(c byte) bool {
	return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')
}

func isDigit(c byte) bool {
	return c >= '0' && c <= '9'
}

// tokensMatch returns true if two tokens are identical and can use tokMatch.
func tokensMatch(a, b nameToken) bool {
	if a.tokType != b.tokType {
		return false
	}
	switch a.tokType {
	case tokChar:
		return a.charVal == b.charVal
	case tokAlpha:
		return a.strVal == b.strVal
	case tokDigits:
		return a.intVal == b.intVal
	case tokDigits0:
		return a.intVal == b.intVal && a.padLen == b.padLen
	}
	return false
}

// canDeltaToken returns true if a numeric delta can encode the difference.
func canDeltaToken(a, b nameToken) bool {
	if a.tokType == tokDigits && b.tokType == tokDigits {
		delta := a.intVal - b.intVal
		return delta < 256 // fits in a byte
	}
	if a.tokType == tokDigits0 && b.tokType == tokDigits0 && a.padLen == b.padLen {
		delta := a.intVal - b.intVal
		return delta < 256
	}
	return false
}

// reformatNames joins names with null terminators (for comparison with decoder output).
func reformatNames(names []string) []byte {
	return []byte(strings.Join(names, "\x00") + "\x00")
}
