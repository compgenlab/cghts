# Native CRAM Reader — TODO

## Must complete before merge

### Integration
- [ ] Wire `cram.Reader` into `NewSamReader()` so `.cram` files use the native reader instead of `SamtoolsSamReader`
- [ ] Verify that all existing callers of `NewSamReader()` with CRAM files work with the native reader

### Correctness
- [x] Validate CRC32 checksums on blocks and container headers
- [x] CRAM v2.1 support (no CRC32, itf8 record counter, different EOF)
- [ ] Verify handling of unmapped/unplaced reads (refID = -1, no reference sequence)
- [ ] Verify embedded reference sequences in CRAM containers (ref stored in the file itself, not external FASTA)

### CRAM v3.1 codec support
- [x] rANS Nx16 (method 5) — order-0 and order-1 working, with 4-way and 32-way (X32) interleaving
- [x] Name tokenizer (method 8) — fixed enum ordering, passes 200-read test
- [ ] fqzcomp (method 7) — quality score compression, not yet implemented. Note: samtools uses rANS Nx16 (not fqzcomp) for quality scores with synthetic data; fqzcomp may only be used with real sequencing data patterns.
- [ ] Adaptive arithmetic coder (method 6) — rarely used, low priority

### Not needed for merge
- [ ] Tag round-trip testing (deferred)
- [ ] Native CRAM writer (currently CRAM writing goes through samtools)
- [ ] Performance benchmarks vs samtools

## Debug notes (resolved)

### rANS Nx16 order-1 (FIXED)
Two bugs found:
1. **Output layout**: Order-1 uses NX-way split (each state decodes 1/NX of output), not round-robin like order-0. Fixed by using quarter/32-way split output indices.
2. **X32 flag ignored**: The X32 flag (0x04) means 32-way interleaving (32 states, 128 bytes for init). Was hardcoded to 4 states. Fixed by passing NX parameter based on X32 flag.
3. **Renormalization**: Order-1 was using 8-bit byte-by-byte renorm but Nx16 uses 16-bit LE uint16 renorm. Fixed to match order-0.

### Name tokenizer (FIXED)
Token type enum values were in wrong order. Corrected to match htscodecs.
