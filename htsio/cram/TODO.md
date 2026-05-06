# Native CRAM Reader — TODO

## Must complete before merge

### Integration
- [ ] Wire `cram.Reader` into `NewSamReader()` so `.cram` files use the native reader instead of `SamtoolsSamReader`
- [ ] Verify that all existing callers of `NewSamReader()` with CRAM files work with the native reader

### Correctness
- [ ] Validate CRC32 checksums on blocks (`block.go` skips this currently)
- [ ] Tag round-trip testing — current tests only compare core fields (first 11 SAM columns), not tags
- [ ] Verify handling of unmapped/unplaced reads (refID = -1, no reference sequence)
- [ ] Verify embedded reference sequences in CRAM containers (ref stored in the file itself, not external FASTA)

### Writing
- [ ] Native CRAM writer (currently CRAM writing goes through samtools)

## Nice to have (not blocking merge)
- [ ] CRAM v2 support
- [ ] Performance benchmarks vs samtools
