// Package codec implements the block-compression codecs used by the CRAM
// sequence-alignment format, plus the entropy primitives they are built from.
//
// The package is a clean-room Go port of the codecs in the htscodecs reference
// library (see ATTRIBUTION.md). It aims for bit-exact interoperability with
// samtools/htslib: any block samtools accepts must decode here, and blocks
// produced here are validated against samtools. The range coder follows Eugene
// Shelwien's public-domain carryless coder, the same one htscodecs uses.
//
// # Codecs
//
// Each codec corresponds to a CRAM block compression method (see the Method
// constants) and exposes plain []byte encode and/or decode entry points:
//
//   - rANS 4x8 (method 4): order-0 and order-1 range-asymmetric numeral
//     systems with four interleaved 8-bit-renormalized states. See
//     [EncodeRans4x8] and [DecodeRans4x8].
//   - rANS Nx16 (method 5): the CRAM v3.1 rANS variant with 16-bit
//     renormalization and optional PACK, RLE, STRIPE, and CAT transforms. See
//     [EncodeRansNx16] and [DecodeRansNx16].
//   - Adaptive arithmetic coding (method 6): order-0/order-1 modelling over a
//     range coder, with the same PACK/RLE/STRIPE/CAT transforms plus an EXT
//     (bzip2) escape. Decode only: see [DecodeArithDynamic].
//   - fqzcomp (method 7): a context-modelling codec specialised for FASTQ
//     quality scores. See [EncodeFqzcomp], [DecodeFqzcomp], and
//     [DecodeFqzcompSize].
//   - Name tokenizer, "tok3" (method 8): splits read names into typed tokens
//     and entropy-codes the per-token streams (with rANS Nx16). See
//     [EncodeNameTokenizer] and [DecodeNameTokenizer].
//
// Internally these share a carryless range coder and a small adaptive
// frequency model.
//
// # Streaming interface
//
// In addition to the []byte functions, [Decoder] and [Encoder] adapt several
// codecs to io.Reader and io.Writer. Each value handles one complete
// compressed block; the streaming wrappers do not perform incremental
// compression across blocks. A decoder reads and decodes the whole block on the
// first Read; an encoder buffers all writes and emits the compressed block on
// Close.
//
// Decoder usage:
//
//	decoder := codec.NewRans4x8Decoder(compressedReader)
//	decoded, err := io.ReadAll(decoder)
//
// Encoder usage:
//
//	encoder := codec.NewRans4x8Encoder(outputWriter, codec.Order0)
//	encoder.Write(data)
//	encoder.Close() // compresses and flushes to outputWriter
package codec
