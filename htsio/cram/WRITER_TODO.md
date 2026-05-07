# Native CRAM Writer — TODO

## API

```go
type WriterOpts struct { ... }

func NewWriterOpts() *WriterOpts                          // defaults: v3.1, level 6
func (o *WriterOpts) Version(v Version) *WriterOpts       // V2, V3, V31
func (o *WriterOpts) Reference(path string) *WriterOpts
func (o *WriterOpts) Level(n int) *WriterOpts             // 0-9, controls which methods tried
func (o *WriterOpts) RecordsPerSlice(n int) *WriterOpts   // default 10000

func NewWriter(filename string, header *htsio.SamHeader, opts ...*WriterOpts) (*Writer, error)
// Writer implements htsio.SamWriter (Write + Close)
// Also support io.Writer and "-" for stdout
```

## Phase 1: Core encoder with gzip compression

### File structure
- [ ] Write file definition: magic "CRAM" + version bytes + 20-byte file ID
- [ ] Write header container: SAM header text as a single block in a container
- [ ] Write EOF container on Close (version-dependent format)
- [ ] Version-specific container headers: CRC32 for v3+, ITF8 record counter for v2, LTF8 for v3+

### Record encoding
- [ ] Buffer records, flush as containers when reaching RecordsPerSlice threshold
- [ ] Group records by reference sequence within a container
- [ ] Build substitution matrix from reference vs read base mismatches
- [ ] Compute alignment features from CIGAR + SEQ + QUAL:
  - Substitution (base mismatch)
  - Insertion (bases inserted)
  - Deletion (bases deleted from ref)
  - Soft clip
  - Hard clip
  - Read base (for unmapped/no-ref)
  - Quality score
  - Reference skip (N in CIGAR)

### Compression header
- [ ] Build encoding map (data series → codec + external block ID)
- [ ] Initial encoding strategy (Phase 1 — simple):
  - BF (bit flags): EXTERNAL
  - CF (CRAM flags): EXTERNAL
  - RI (ref ID): EXTERNAL
  - RL (read length): EXTERNAL
  - AP (alignment position): EXTERNAL (delta-encoded)
  - RG (read group): EXTERNAL
  - RN (read name): BYTE_ARRAY_STOP (NUL-terminated, external block)
  - MF (mate flags): EXTERNAL
  - NS (next segment ref): EXTERNAL
  - NP (next segment pos): EXTERNAL
  - TS (template size): EXTERNAL
  - NF (records to next fragment): EXTERNAL
  - BA (base substitution codes): EXTERNAL
  - QS (quality scores): EXTERNAL
  - IN (insertion bases): BYTE_ARRAY_STOP
  - SC (soft clip bases): BYTE_ARRAY_STOP
  - DL (deletion length): EXTERNAL
  - RS (reference skip length): EXTERNAL
  - HC (hard clip length): EXTERNAL
  - MQ (mapping quality): EXTERNAL
  - TL (tag line/combo ID): EXTERNAL
  - Tag values: BYTE_ARRAY_LEN per tag type, external blocks
- [ ] Write substitution matrix
- [ ] Write tag encoding map (one entry per tag:type combo seen)

### Slice/block encoding
- [ ] Core block: bit-level encoding of Huffman/Beta coded fields
- [ ] External blocks: one per data series, gzip compressed
- [ ] Block headers: method, content type, content ID, compressed/uncompressed sizes
- [ ] CRC32 on blocks and container headers (v3+)

### Tag encoding
- [ ] Group tags by key:type (e.g., "NM:i", "RG:Z")
- [ ] Assign tag line combos (TL field): each unique set of tag keys → integer ID
- [ ] Encode tag values into per-tag-type external blocks
- [ ] Handle all SAM tag types: A, i, f, Z, H, B arrays

### Reference handling
- [ ] Load reference sequences via referenceProvider (existing code)
- [ ] Compute MD5 for reference slices (needed for container header)
- [ ] Support embedded reference mode (ref stored in container, not external)
- [ ] Handle unmapped reads (refID = -1)

## Phase 2: rANS 4x8 encoder

- [ ] Implement rANS order-0 encoder
- [ ] Implement rANS order-1 encoder
- [ ] Frequency table construction from input data
- [ ] Competitive compression: try gzip vs rANS, keep smaller
- [ ] Covers v2 and v3.0 well

## Phase 3: v3.1 codec encoders + competitive compression

### rANS Nx16 encoder
- [ ] Order-0 with 4-way and 32-way (X32) interleaving
- [ ] Order-1 with NX-way split output
- [ ] PACK, RLE, STRIPE, CAT transforms
- [ ] 16-bit LE renormalization

### fqzcomp encoder (quality scores)
- [ ] Context model parameter selection
- [ ] Range coder encoder (reverse of decoder)
- [ ] Multiple parameter sets (selector table)
- [ ] Try fqzcomp vs other methods for quality blocks

### Name tokenizer encoder (read names)
- [ ] Tokenize read names into typed tokens
- [ ] Encode token streams
- [ ] Try name tokenizer vs gzip/bzip2 for name blocks

### Adaptive arithmetic encoder
- [ ] Simple model encoder (reverse of decoder)
- [ ] Order-0 and order-1 with RLE
- [ ] PACK/STRIPE/CAT transforms

### Competitive compression
- [ ] Per-block: try all enabled methods, keep smallest output
- [ ] Method selection based on compression level:
  - Level 0: no compression
  - Level 1: gzip_1 + rANS order-0 (fast)
  - Level 3-5: gzip + rANS + bzip2
  - Level 6+: all methods including lzma, full rANS Nx16 variants, fqzcomp, arith
- [ ] Quality scores: prefer fqzcomp at higher levels (v3.1)
- [ ] Read names: prefer name tokenizer (v3.1), exclude rANS
- [ ] Stats/metrics tracking to avoid re-trying methods that consistently lose

## Notes

### htslib reference behavior
- htslib tries ALL enabled compression methods per block and picks smallest
- Available methods depend on version + compression level
- Quality scores get fqzcomp (v3.1 only)
- Read names get name tokenizer (v3.1), NOT rANS
- Core block: gzip only, low level, only if >500 bytes

### Container structure (for reference)
```
Container header:
  - length, ref_seq_id, start_pos, alignment_span
  - num_records, record_counter, num_bases, num_blocks, num_landmarks
  - CRC32 (v3+)

Compression header block:
  - preservation map (substitution matrix, tag IDs, etc.)
  - data series encoding map
  - tag encoding map

Slice header block:
  - ref_seq_id, start_pos, alignment_span, num_records
  - record_counter, external block content IDs, ref_bases_md5

Core data block:
  - bit-level encoded fields (Huffman, Beta)

External blocks (one per data series):
  - compressed data for that series
```
