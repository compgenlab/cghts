package codec

// FQZcomp quality score encoder (CRAM v3.1 method 7).
// Encodes quality scores using adaptive arithmetic coding with context modeling.
//
// This is the reverse of the decoder in fqzcomp.go. Both share the same
// data structures (fqzParam, fqzGParams, fqzState, fqzModels) and context
// update logic (fqzUpdateCtx).

// EncodeFqzcomp encodes quality scores using the fqzcomp codec.
// quals is the raw quality values (0-based, not ASCII+33).
// readLengths is the length of each read in order.
func EncodeFqzcomp(quals []byte, readLengths []int) []byte {
	if len(quals) == 0 {
		return nil
	}

	// Analyze quality distribution for parameter selection.
	maxSym := 0
	for _, q := range quals {
		if int(q) >= maxSym {
			maxSym = int(q) + 1
		}
	}

	// Build simple single-param configuration.
	gp := fqzBuildParams(maxSym, readLengths)

	// Serialize parameters BEFORE pre-shifting (file stores unshifted values).
	paramBytes := fqzWriteParameters(&gp)

	// Pre-shift tables (must match decoder behavior).
	for i := 0; i < gp.nparam; i++ {
		pm := &gp.p[i]
		for j := 0; j < 1024; j++ {
			pm.ptab[j] <<= pm.ploc
		}
		for j := 0; j < 256; j++ {
			pm.dtab[j] <<= pm.dloc
		}
	}

	// Create models (same as decoder).
	models := fqzCreateModels(&gp)

	// Initialize range encoder.
	re := newRangeEncoder()

	// Encode quality symbols.
	state := fqzState{
		firstLen: true,
		ctx:      uint(gp.p[0].context),
	}
	pm := &gp.p[0]
	last := state.ctx

	readIdx := 0
	for i := 0; i < len(quals); {
		if state.p == 0 {
			// Start of new read — encode read header.
			readLen := readLengths[readIdx]
			readIdx++

			// Encode selector (always 0 for single param).
			if pm.doSel || gp.gflags&gflagMultiParam != 0 {
				models.sel.encodeSymbol(re, uint16(state.s))
			}

			// Encode read length.
			if !pm.fixedLen || state.firstLen {
				models.len[0].encodeSymbol(re, uint16(readLen&0xFF))
				models.len[1].encodeSymbol(re, uint16((readLen>>8)&0xFF))
				models.len[2].encodeSymbol(re, uint16((readLen>>16)&0xFF))
				models.len[3].encodeSymbol(re, uint16((readLen>>24)&0xFF))
				state.firstLen = false
				state.lastLen = uint(readLen)
			}

			// No reverse, no dedup in simple mode.

			state.rec++
			state.p = uint(readLen)
			state.delta = 0
			state.prevq = 0
			state.qctx = 0
			state.ctx = uint(pm.context)
			last = state.ctx
		}

		// Encode quality symbols for current read.
		for state.p != 0 && i < len(quals) {
			q := uint(quals[i])
			models.qual[last].encodeSymbol(re, uint16(q))
			last = fqzUpdateCtx(pm, &state, q)
			i++
		}
	}

	// Finalize range coder.
	encoded := re.finish()

	// Build output: size + params + encoded data.
	var out []byte
	out = varPutU32Slice(out, uint32(len(quals)))
	out = append(out, paramBytes...)
	out = append(out, encoded...)
	return out
}

// fqzBuildParams creates a simple single-param configuration.
func fqzBuildParams(maxSym int, readLengths []int) fqzGParams {
	// Check if all reads have the same length.
	fixedLen := true
	if len(readLengths) > 1 {
		for _, l := range readLengths[1:] {
			if l != readLengths[0] {
				fixedLen = false
				break
			}
		}
	}

	pm := fqzParam{
		context:  0,
		pflags:   0,
		maxSym:   maxSym,
		qbits:    0,
		qshift:   0,
		qloc:     0,
		sloc:     0,
		ploc:     0,
		dloc:     0,
		qmask:    0,
		fixedLen: fixedLen,
	}

	if fixedLen {
		pm.pflags |= pflagDoLen
	}

	// Identity quality map (no remapping).
	for i := 0; i < 256; i++ {
		pm.qmap[i] = uint(i)
	}

	// Identity qtab.
	for i := 0; i < 256; i++ {
		pm.qtab[i] = uint(i)
	}

	// Zero dtab and ptab (not using position or delta context in simple mode).
	// dtab and ptab are all zero by default.

	return fqzGParams{
		vers:   fqzVers,
		gflags: 0, // no multi-param, no stab, no reverse
		nparam: 1,
		maxSym: maxSym,
		p:      []fqzParam{pm},
	}
}

// positionBucket maps a position (0-1023) to a small bucket value for context.
func positionBucket(pos int) int {
	if pos < 2 {
		return pos
	}
	// Log2-ish bucketing.
	b := 1
	for v := pos; v > 1; v >>= 1 {
		b++
	}
	return b
}

// fqzWriteParameters serializes the global and per-param configuration.
func fqzWriteParameters(gp *fqzGParams) []byte {
	var out []byte

	// Version.
	out = append(out, byte(gp.vers))

	// Global flags.
	out = append(out, byte(gp.gflags))

	// Number of params (only if multi-param).
	if gp.gflags&gflagMultiParam != 0 {
		out = append(out, byte(gp.nparam))
	}

	// Selector table (only if have stab).
	if gp.gflags&gflagHaveStab != 0 {
		out = append(out, byte(gp.maxSel))
		out = append(out, fqzWriteArray(gp.stab[:], 256)...)
	}

	// Per-param blocks.
	for i := 0; i < gp.nparam; i++ {
		out = append(out, fqzWriteParam1(&gp.p[i])...)
	}

	return out
}

// fqzWriteParam1 serializes one parameter block.
func fqzWriteParam1(pm *fqzParam) []byte {
	var out []byte

	// Starting context (LE uint16).
	out = append(out, byte(pm.context), byte(pm.context>>8))

	// Flags.
	out = append(out, byte(pm.pflags))

	// Max symbol.
	out = append(out, byte(pm.maxSym))

	// qbits/qshift.
	out = append(out, byte(pm.qbits<<4)|byte(pm.qshift&0xF))

	// qloc/sloc.
	out = append(out, byte(pm.qloc<<4)|byte(pm.sloc&0xF))

	// ploc/dloc.
	out = append(out, byte(pm.ploc<<4)|byte(pm.dloc&0xF))

	// Quality map.
	if pm.storeQmap {
		for i := 0; i < pm.maxSym; i++ {
			out = append(out, byte(pm.qmap[i]))
		}
	}

	// Quality context table.
	if pm.useQtab {
		out = append(out, fqzWriteArray(pm.qtab[:], 256)...)
	}

	// Position table.
	if pm.usePtab {
		out = append(out, fqzWriteArray(pm.ptab[:], 1024)...)
	}

	// Delta table.
	if pm.useDtab {
		out = append(out, fqzWriteArray(pm.dtab[:], 256)...)
	}

	return out
}

// fqzWriteArray encodes an array using double-RLE (reverse of fqzReadArray).
// The array maps positions → symbol values. The encoding is:
//   Inner layer: for each symbol 0,1,2,..., count its run length.
//     Runs > 255 are split into 255-byte continuation segments.
//   Outer layer: compress the resulting byte stream with consecutive-duplicate RLE.
func fqzWriteArray(array []uint, size int) []byte {
	// Step 1: Build inner RLE parts.
	// For each symbol value in order, count how many consecutive positions have that value.
	var parts []int
	pos := 0
	for sym := 0; pos < size; sym++ {
		// Count run of this symbol.
		runLen := 0
		for pos+runLen < size && array[pos+runLen] == uint(sym) {
			runLen++
		}
		pos += runLen

		// Encode run as 255-max continuation segments.
		// A run of 0 still emits one part (0).
		for {
			seg := runLen
			if seg > 255 {
				seg = 255
			}
			parts = append(parts, seg)
			runLen -= seg
			if seg < 255 || runLen == 0 {
				break
			}
		}

		if sym > 255 {
			break
		}
	}

	// Step 2: Outer RLE — compress consecutive identical values.
	// If value == previous value, write a copies byte (0-255 additional copies).
	var out []byte
	last := -1
	for i := 0; i < len(parts); i++ {
		v := parts[i]
		out = append(out, byte(v))
		if v == last {
			// Count additional consecutive identical values.
			copies := 0
			for i+1 < len(parts) && parts[i+1] == v && copies < 255 {
				copies++
				i++
			}
			out = append(out, byte(copies))
		}
		last = v
	}

	return out
}
