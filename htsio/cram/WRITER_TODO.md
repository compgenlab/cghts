# Native CRAM Writer — TODO

## API

```go
type WriterOpts struct { ... }

func NewWriterOpts() *WriterOpts                          // defaults: v3.1, level 6
func (o *WriterOpts) SetVersion(v Version) *WriterOpts    // V2, V3, V31
func (o *WriterOpts) Reference(path string) *WriterOpts
func (o *WriterOpts) Level(n int) *WriterOpts             // 0-9, controls which methods tried
func (o *WriterOpts) RecordsPerSlice(n int) *WriterOpts   // default 10000
func (o *WriterOpts) EmbedRef(v bool) *WriterOpts         // store ref in CRAM file

func NewWriter(filename string, header *htsio.SamHeader, opts ...*WriterOpts) (*Writer, error)
// Writer implements htsio.SamWriter (Write + Close)
// Also support io.Writer and "-" for stdout
```

## Phase 1: Core encoder with gzip compression — DONE

### File structure
- [x] Write file definition: magic "CRAM" + version bytes + 20-byte file ID
- [x] Write header container: SAM header text as a single block in a container
- [x] Write EOF container on Close (version-dependent format)
- [x] Version-specific container headers: CRC32 for v3+, ITF8 record counter for v2, LTF8 for v3+

### Record encoding
- [x] Buffer records, flush as containers when reaching RecordsPerSlice threshold
- [x] Group records by reference sequence within a container
- [x] Build substitution matrix from reference vs read base mismatches
- [x] Compute alignment features from CIGAR + SEQ + QUAL:
  - Substitution (base mismatch)
  - Insertion (bases inserted)
  - Deletion (bases deleted from ref)
  - Soft clip
  - Hard clip
  - Read base (for unmapped/no-ref)
  - Quality score
  - Reference skip (N in CIGAR)

### Compression header
- [x] Build encoding map (data series → codec + external block ID)
- [x] Initial encoding strategy (Phase 1 — simple):
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
- [x] Write substitution matrix
- [x] Write tag encoding map (one entry per tag:type combo seen)

### Slice/block encoding
- [x] Core block: bit-level encoding of Huffman/Beta coded fields
- [x] External blocks: one per data series, gzip compressed
- [x] Block headers: method, content type, content ID, compressed/uncompressed sizes
- [x] CRC32 on blocks and container headers (v3+)

### Tag encoding
- [x] Group tags by key:type (e.g., "NM:i", "RG:Z")
- [x] Assign tag line combos (TL field): each unique set of tag keys → integer ID
- [x] Encode tag values into per-tag-type external blocks
- [x] Handle all SAM tag types: A, i, f, Z, H, B arrays

### Reference handling
- [x] Load reference sequences via referenceProvider (existing code)
- [x] Compute MD5 for reference slices (needed for container header)
- [x] Support embedded reference mode (ref stored in container, not external)
- [x] Handle unmapped reads (refID = -1)

## Phase 2: rANS 4x8 encoder — DONE

- [x] Implement rANS order-0 encoder
- [x] Implement rANS order-1 encoder
- [x] Frequency table construction from input data
- [x] Competitive compression: try gzip vs rANS, keep smaller
- [x] Covers v2 and v3.0 well

## Phase 3: v3.1 codec encoders + competitive compression

### rANS Nx16 encoder
- [x] Order-0 with 4-way interleaving + 16-bit LE renormalization
- [x] Order-1 with NX-way split output + freqD frequency table format
- [x] PACK transform (symbol remapping + bit packing for ≤16 symbols)
- [ ] RLE transform (run-length encoding for repeated symbols)
- [ ] STRIPE transform (byte-interleaving into N sub-streams)
- [ ] CAT passthrough

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
- [x] Per-block: try all enabled methods, keep smallest output
- [x] Method selection based on version:
  - v2.1: raw + gzip
  - v3.0: raw + gzip + rANS 4x8
  - v3.1: raw + gzip + rANS 4x8 + rANS Nx16 (order-0, order-1, PACK)
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
