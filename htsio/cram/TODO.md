# Native CRAM Reader — TODO

## Must complete before merge

### Integration
- [ ] Wire `cram.Reader` into `NewSamReader()` so `.cram` files use the native reader instead of `SamtoolsSamReader`
- [ ] Verify that all existing callers of `NewSamReader()` with CRAM files work with the native reader

### Correctness
- [x] Validate CRC32 checksums on blocks and container headers
- [ ] Tag round-trip testing — current tests only compare core fields (first 11 SAM columns), not tags
- [ ] Verify handling of unmapped/unplaced reads (refID = -1, no reference sequence)
- [ ] Verify embedded reference sequences in CRAM containers (ref stored in the file itself, not external FASTA)

### CRAM v3.1 codec support
- [x] rANS Nx16 (method 5) — order-0, order-1, with PACK, RLE, STRIPE, CAT, NOSIZE transforms
- [ ] Name tokenizer (method 8) — initial implementation exists but has a bug with DIGITS0/DZLen descriptor mapping; needs debugging against htscodecs reference (DZLen data appears to be stored at a different descriptor index than expected)
- [ ] fqzcomp (method 7) — quality score compression, not yet implemented
- [ ] Adaptive arithmetic coder (method 6) — rarely used, low priority

### Writing
- [ ] Native CRAM writer (currently CRAM writing goes through samtools)

## Nice to have (not blocking merge)
- [ ] CRAM v2 support
- [ ] Performance benchmarks vs samtools

## Debug notes

### Name tokenizer DZLen issue
When decoding DIGITS0 tokens, the C code reads the zero-pad width from descriptor `(ntok<<4)|12` (DZLen=12). In the test CRAM 3.1 file (200 reads), token 2 has 109 DIGITS0 entries but no descriptor at type 12. There IS a descriptor at type 4 (DIGITS) with 109 single-byte values, which may be the widths stored under a different type index. Need to check the htscodecs tok3 encoder to understand the actual descriptor numbering.

Observed descriptor layout for "readXXXX" names:
- Token 0: type=DDELTA(6), 800 bytes (200×4 uint32)
- Token 1: type=ALPHA(1)/MATCH(7), 5 bytes ("read")
- Token 2: type=DIGITS0(3)/END(9), desc[35]=436B, desc[36]=109B, desc[41]=91B
- Token 3: type=DZLEN(12) for all 200 reads (acts as END via default case)
