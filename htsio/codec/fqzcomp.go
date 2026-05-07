package codec

import (
	"fmt"
)

// FQZcomp quality score decoder (CRAM v3.1 method 7).
// Based on the htscodecs fqzcomp_qual.c reference implementation.

const (
	fqzVers    = 5
	fqzCtxBits = 16
	fqzCtxSize = 1 << fqzCtxBits
	fqzQMax    = 256
)

// Global flags.
const (
	gflagMultiParam = 1
	gflagHaveStab   = 2
	gflagDoRev      = 4
)

// Per-param flags.
const (
	pflagDoDedup  = 2
	pflagDoLen    = 4
	pflagDoSel    = 8
	pflagHaveQmap = 16
	pflagHavePtab = 32
	pflagHaveDtab = 64
	pflagHaveQtab = 128
)

type fqzParam struct {
	context  uint16
	pflags   uint
	doSel    bool
	doDedup  bool
	storeQmap bool
	fixedLen bool
	useQtab  bool
	useDtab  bool
	usePtab  bool

	qbits, qloc uint
	ploc, dloc  uint
	sloc        uint
	qshift      uint
	qmask       uint

	maxSym int

	qmap [256]uint
	qtab [256]uint
	ptab [1024]uint
	dtab [256]uint
}

type fqzGParams struct {
	vers    int
	gflags  uint
	nparam  int
	maxSel  int
	stab    [256]uint
	maxSym  int
	p       []fqzParam
}

type fqzState struct {
	qctx     uint
	p        uint // position (bytes remaining in current read)
	delta    uint
	prevq    uint
	s        uint // selector
	firstLen bool
	lastLen  uint
	rec      int
	ctx      uint
}

type fqzModels struct {
	qual    []*simpleModel // fqzCtxSize models for quality
	len     [4]*simpleModel
	revcomp *simpleModel
	sel     *simpleModel
	dup     *simpleModel
}

func DecodeFqzcomp(data []byte) ([]byte, error) {
	if len(data) < 4 {
		return nil, fmt.Errorf("fqzcomp: data too short")
	}

	// Read uncompressed size.
	n, outLen := varGetU32(data)
	if n == 0 {
		return nil, fmt.Errorf("fqzcomp: truncated size")
	}
	pos := n

	// Read parameters.
	gp, consumed, err := fqzReadParameters(data[pos:])
	if err != nil {
		return nil, fmt.Errorf("fqzcomp: %w", err)
	}
	pos += consumed

	// Pre-shift tables (optimization from reference).
	for i := 0; i < gp.nparam; i++ {
		pm := &gp.p[i]
		for j := 0; j < 1024; j++ {
			pm.ptab[j] <<= pm.ploc
		}
		for j := 0; j < 256; j++ {
			pm.dtab[j] <<= pm.dloc
		}
	}

	// Create models.
	models := fqzCreateModels(&gp)

	// Initialize range coder.
	rc := newRangeDecoder(data[pos:])

	// Decode.
	out := make([]byte, outLen)
	state := fqzState{
		firstLen: true,
		ctx:      uint(gp.p[0].context),
	}

	revA := make([]byte, 1000)
	lenA := make([]int, 1000)

	var last uint
	pm := &gp.p[0]

	for i := uint32(0); i < outLen; {
		// Grow arrays if needed.
		if state.rec >= len(revA) {
			newRevA := make([]byte, len(revA)*2)
			copy(newRevA, revA)
			revA = newRevA
			newLenA := make([]int, len(lenA)*2)
			copy(newLenA, lenA)
			lenA = newLenA
		}

		if state.p == 0 {
			// Start of a new read.
			skip, err := fqzDecompressNewRead(&state, &gp, pm, models, rc,
				out, outLen, &i, revA, lenA)
			if err != nil {
				return nil, fmt.Errorf("fqzcomp: %w", err)
			}
			if skip {
				continue
			}
			// Update pm based on selector.
			x := state.s
			if gp.gflags&gflagHaveStab != 0 {
				if x > 255 {
					x = 255
				}
				x = gp.stab[x]
			}
			if int(x) >= gp.nparam {
				return nil, fmt.Errorf("fqzcomp: selector %d out of range", x)
			}
			pm = &gp.p[x]
			last = state.ctx
		}

		// Decode quality symbols for current read.
		for state.p != 0 && i < outLen {
			q := uint(models.qual[last].decodeSymbol(rc))
			last = fqzUpdateCtx(pm, &state, q)
			out[i] = byte(pm.qmap[q])
			i++
		}
	}

	// Record final entry.
	rec := state.rec
	if rec >= len(revA) {
		newRevA := make([]byte, rec+1)
		copy(newRevA, revA)
		revA = newRevA
		newLenA := make([]int, rec+1)
		copy(newLenA, lenA)
		lenA = newLenA
	}

	// Reverse reads if needed.
	if gp.gflags&gflagDoRev != 0 {
		var idx uint32
		for r := 0; idx < outLen && r < len(lenA) && r <= state.rec; r++ {
			rlen := lenA[r]
			if revA[r] != 0 {
				// Reverse this read's quality.
				for a, b := int(idx), int(idx)+rlen-1; a < b; a, b = a+1, b-1 {
					out[a], out[b] = out[b], out[a]
				}
			}
			idx += uint32(rlen)
		}
	}

	if !rc.finish() {
		return nil, fmt.Errorf("fqzcomp: range coder error")
	}

	return out, nil
}

func fqzCreateModels(gp *fqzGParams) *fqzModels {
	m := &fqzModels{}
	m.qual = make([]*simpleModel, fqzCtxSize)
	for i := range m.qual {
		m.qual[i] = newSimpleModel(fqzQMax, gp.maxSym+1)
	}
	for i := 0; i < 4; i++ {
		m.len[i] = newSimpleModel(256, 256)
	}
	m.revcomp = newSimpleModel(2, 2)
	if gp.maxSel > 0 {
		m.sel = newSimpleModel(256, gp.maxSel+1)
	} else {
		m.sel = newSimpleModel(256, 1)
	}
	m.dup = newSimpleModel(2, 2)
	return m
}

func fqzUpdateCtx(pm *fqzParam, state *fqzState, q uint) uint {
	state.qctx = (state.qctx << pm.qshift) + pm.qtab[q]
	last := (state.qctx & pm.qmask) << pm.qloc

	p := state.p
	if p > 1023 {
		p = 1023
	}
	last += pm.ptab[p]

	d := state.delta
	if d > 255 {
		d = 255
	}
	last += pm.dtab[d]

	last += state.s << pm.sloc

	if state.prevq != q {
		state.delta++
	}
	state.prevq = q
	state.p--

	return last & (fqzCtxSize - 1)
}

func fqzDecompressNewRead(state *fqzState, gp *fqzGParams, pm *fqzParam,
	models *fqzModels, rc *rangeCoder,
	uncomp []byte, outSize uint32, iPtr *uint32,
	revA []byte, lenA []int) (bool, error) {

	i := *iPtr

	if pm.doSel || (gp.gflags&gflagMultiParam != 0) {
		state.s = uint(models.sel.decodeSymbol(rc))
	} else {
		state.s = 0
	}

	x := state.s
	if gp.gflags&gflagHaveStab != 0 {
		if x > 255 {
			x = 255
		}
		x = gp.stab[x]
	}
	if int(x) >= gp.nparam {
		return false, fmt.Errorf("selector %d out of range (nparam=%d)", x, gp.nparam)
	}
	pm = &gp.p[x]

	readLen := state.lastLen
	if !pm.fixedLen || state.firstLen {
		readLen = uint(models.len[0].decodeSymbol(rc))
		readLen |= uint(models.len[1].decodeSymbol(rc)) << 8
		readLen |= uint(models.len[2].decodeSymbol(rc)) << 16
		readLen |= uint(models.len[3].decodeSymbol(rc)) << 24
		state.firstLen = false
		state.lastLen = readLen
	}
	if readLen == 0 || uint32(readLen) > outSize-i {
		return false, fmt.Errorf("invalid read length %d", readLen)
	}

	if gp.gflags&gflagDoRev != 0 {
		rev := byte(models.revcomp.decodeSymbol(rc))
		if state.rec < len(revA) {
			revA[state.rec] = rev
			lenA[state.rec] = int(readLen)
		}
	}

	if pm.doDedup {
		if models.dup.decodeSymbol(rc) != 0 {
			// Duplicate of previous read.
			if uint32(readLen) > i {
				return false, fmt.Errorf("dedup: read length %d > position %d", readLen, i)
			}
			copy(uncomp[i:i+uint32(readLen)], uncomp[i-uint32(readLen):i])
			*iPtr = i + uint32(readLen)
			state.p = 0
			state.rec++
			return true, nil // skip
		}
	}

	state.rec++
	state.p = readLen
	state.delta = 0
	state.prevq = 0
	state.qctx = 0
	state.ctx = uint(pm.context)

	*iPtr = i
	return false, nil
}

func fqzReadParameters(data []byte) (fqzGParams, int, error) {
	var gp fqzGParams
	if len(data) < 2 {
		return gp, 0, fmt.Errorf("params too short")
	}

	pos := 0

	// Version.
	gp.vers = int(data[pos])
	pos++
	if gp.vers != fqzVers {
		return gp, 0, fmt.Errorf("unsupported fqzcomp version %d (want %d)", gp.vers, fqzVers)
	}

	// Global flags.
	gp.gflags = uint(data[pos])
	pos++

	// Number of param blocks.
	if gp.gflags&gflagMultiParam != 0 {
		if pos >= len(data) {
			return gp, 0, fmt.Errorf("truncated nparam")
		}
		gp.nparam = int(data[pos])
		pos++
	} else {
		gp.nparam = 1
	}
	if gp.nparam <= 0 {
		return gp, 0, fmt.Errorf("invalid nparam %d", gp.nparam)
	}
	gp.maxSel = 0
	if gp.nparam > 1 {
		gp.maxSel = gp.nparam
	}

	// Selector table.
	if gp.gflags&gflagHaveStab != 0 {
		if pos >= len(data) {
			return gp, 0, fmt.Errorf("truncated stab")
		}
		gp.maxSel = int(data[pos])
		pos++
		used := fqzReadArray(data[pos:], gp.stab[:], 256)
		if used < 0 {
			return gp, 0, fmt.Errorf("failed to read stab")
		}
		pos += used
	} else {
		for i := 0; i < gp.nparam; i++ {
			gp.stab[i] = uint(i)
		}
		for i := gp.nparam; i < 256; i++ {
			gp.stab[i] = uint(gp.nparam - 1)
		}
	}

	// Read per-param blocks.
	gp.p = make([]fqzParam, gp.nparam)
	gp.maxSym = 0
	for i := 0; i < gp.nparam; i++ {
		consumed, err := fqzReadParam1(&gp.p[i], data[pos:])
		if err != nil {
			return gp, 0, fmt.Errorf("param %d: %w", i, err)
		}
		if gp.p[i].doSel && gp.maxSel == 0 {
			return gp, 0, fmt.Errorf("param %d: do_sel without max_sel", i)
		}
		pos += consumed
		if gp.maxSym < gp.p[i].maxSym {
			gp.maxSym = gp.p[i].maxSym
		}
	}

	return gp, pos, nil
}

func fqzReadParam1(pm *fqzParam, data []byte) (int, error) {
	if len(data) < 7 {
		return 0, fmt.Errorf("param too short")
	}

	pos := 0

	// Starting context.
	pm.context = uint16(data[pos]) | uint16(data[pos+1])<<8
	pos += 2

	// Flags.
	pm.pflags = uint(data[pos])
	pos++
	pm.useQtab = pm.pflags&pflagHaveQtab != 0
	pm.useDtab = pm.pflags&pflagHaveDtab != 0
	pm.usePtab = pm.pflags&pflagHavePtab != 0
	pm.doSel = pm.pflags&pflagDoSel != 0
	pm.fixedLen = pm.pflags&pflagDoLen != 0
	pm.doDedup = pm.pflags&pflagDoDedup != 0
	pm.storeQmap = pm.pflags&pflagHaveQmap != 0

	// Max symbol.
	pm.maxSym = int(data[pos])
	pos++

	// Sub-context sizes and locations.
	pm.qbits = uint(data[pos] >> 4)
	pm.qmask = (1 << pm.qbits) - 1
	pm.qshift = uint(data[pos] & 15)
	pos++
	pm.qloc = uint(data[pos] >> 4)
	pm.sloc = uint(data[pos] & 15)
	pos++
	pm.ploc = uint(data[pos] >> 4)
	pm.dloc = uint(data[pos] & 15)
	pos++

	// Quality map.
	if pm.storeQmap {
		if pos+pm.maxSym > len(data) {
			return 0, fmt.Errorf("truncated qmap")
		}
		for i := 0; i < 256; i++ {
			pm.qmap[i] = ^uint(0) // sentinel
		}
		for i := 0; i < pm.maxSym; i++ {
			pm.qmap[i] = uint(data[pos])
			pos++
		}
	} else {
		for i := 0; i < 256; i++ {
			pm.qmap[i] = uint(i)
		}
	}

	// Quality context table.
	if pm.qbits != 0 {
		if pm.useQtab {
			used := fqzReadArray(data[pos:], pm.qtab[:], 256)
			if used < 0 {
				return 0, fmt.Errorf("failed to read qtab")
			}
			pos += used
		} else {
			for i := 0; i < 256; i++ {
				pm.qtab[i] = uint(i)
			}
		}
	}

	// Position table.
	if pm.usePtab {
		used := fqzReadArray(data[pos:], pm.ptab[:], 1024)
		if used < 0 {
			return 0, fmt.Errorf("failed to read ptab")
		}
		pos += used
	}

	// Delta table.
	if pm.useDtab {
		used := fqzReadArray(data[pos:], pm.dtab[:], 256)
		if used < 0 {
			return 0, fmt.Errorf("failed to read dtab")
		}
		pos += used
	}

	return pos, nil
}

// fqzReadArray reads a double-RLE encoded array.
// Returns number of bytes consumed, or -1 on error.
func fqzReadArray(in []byte, array []uint, size int) int {
	if size > 1024 {
		size = 1024
	}

	// Remove level one of run-length encoding.
	R := make([]int, 1024)
	var last, j, z int
	last = -1
	i := 0
	for z < size && i < len(in) {
		run := int(in[i])
		i++
		R[j] = run
		j++
		z += run
		if run == last {
			if i >= len(in) {
				return -1
			}
			copies := int(in[i])
			i++
			z += run * copies
			for copies > 0 && z <= size && j < 1024 {
				R[j] = run
				j++
				copies--
			}
		}
		if j >= 1024 {
			return -1
		}
		last = run
	}
	nb := i

	// Expand inner level of run-length encoding.
	rMax := j
	z = 0
	j = 0
	for sym := 0; j < size; sym++ {
		runLen := 0
		for {
			if z >= rMax {
				return -1
			}
			runPart := R[z]
			z++
			runLen += runPart
			if runPart != 255 || z >= rMax {
				break
			}
		}
		for runLen > 0 && j < size {
			array[j] = uint(sym)
			j++
			runLen--
		}
	}

	return nb
}
