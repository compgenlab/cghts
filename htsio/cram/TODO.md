# Native CRAM Reader — TODO

## Must complete before merge

### Integration
- [x] Wire `cram.Reader` into `NewSamReader()` so `.cram` files use the native reader instead of `SamtoolsSamReader`
- [x] Verify that all existing callers of `NewSamReader()` with CRAM files work with the native reader

### Correctness
- [x] Validate CRC32 checksums on blocks and container headers
- [x] CRAM v2.1 support (no CRC32, itf8 record counter, different EOF)
- [x] Verify handling of unmapped/unplaced reads (refID = -1, no reference sequence)
- [x] Verify embedded reference sequences in CRAM containers (ref stored in the file itself, not external FASTA)

### CRAM v3.1 codec support
- [x] rANS Nx16 (method 5) — order-0 and order-1 working, with 4-way and 32-way (X32) interleaving
- [x] Name tokenizer (method 8) — fixed enum ordering, passes 200-read test
- [x] fqzcomp (method 7) — quality score compression decoder, verified with htscodecs test vectors
- [x] Adaptive arithmetic coder (method 6) — order-0, order-1, with RLE/PACK/STRIPE/CAT transforms, verified with htscodecs test vectors

### Reference providers (seqio.ReferenceReader)
- [x] Indexed FASTA (.fai) with 10MB chunk cache, LRU (1GB max)
- [x] Bgzip-compressed FASTA (.fa.gz + .gzi) via bgzf.IndexedReader
- [x] In-memory FASTA fallback (no .fai, for small refs / tests)
- [x] Auto-detection factory: seqio.OpenReference(path)
- [ ] Remote HTTP/HTTPS reference (range requests for .fai-indexed refs)
- [ ] GA4GH refget API (sequence retrieval by MD5/refget accession)
- [ ] htslib REF_CACHE / REF_PATH directory layout

### Not needed for merge
- [x] Tag round-trip testing — verified all tag types (A, i, f, Z) survive write→read roundtrip across v2.1/v3.0/v3.1
- [x] Native CRAM writer (v2.1, v3.0, v3.1 with gzip, rANS 4x8, and rANS Nx16 competitive compression)
- [ ] Performance benchmarks vs samtools

## Bugs fixed

### Writer: tag block ID assignment (FIXED)
Tag block IDs were re-assigned by iterating a map (non-deterministic order), causing
misalignment with the compression header. Fixed to extract block IDs from encoding
descriptor params.

### Reader: unmapped zero-length reads (FIXED)
`reconstructSequence` returned `""` instead of `"*"` for unmapped reads with no sequence.

## Debug notes (resolved)

### rANS Nx16 order-1 (FIXED)
Two bugs found:
1. **Output layout**: Order-1 uses NX-way split (each state decodes 1/NX of output), not round-robin like order-0. Fixed by using quarter/32-way split output indices.
2. **X32 flag ignored**: The X32 flag (0x04) means 32-way interleaving (32 states, 128 bytes for init). Was hardcoded to 4 states. Fixed by passing NX parameter based on X32 flag.
3. **Renormalization**: Order-1 was using 8-bit byte-by-byte renorm but Nx16 uses 16-bit LE uint16 renorm. Fixed to match order-0.

### Name tokenizer (FIXED)
Token type enum values were in wrong order. Corrected to match htscodecs.
