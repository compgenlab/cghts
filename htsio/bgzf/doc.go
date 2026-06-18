// Package bgzf implements reading and writing of BGZF (Blocked GNU Zip Format)
// block compression as defined in the SAM/BAM specification.
//
// BGZF is a multi-member gzip format. Each member (block) holds at most 64 KiB
// (MaxUncompressedSize) of uncompressed data and carries an extra "BC" subfield
// (BSIZE) recording the total block size minus one. The total compressed block
// size is likewise bounded by 64 KiB (MaxBlockSize). Because blocks are
// independently compressed, BGZF supports random access into the uncompressed
// stream via virtual offsets, and the format is terminated by a fixed 28-byte
// empty EOF block.
//
// A virtual offset ([VirtualOffset]) addresses a position in the uncompressed
// stream by combining the compressed file offset of the enclosing block (upper
// 48 bits) with the uncompressed offset within that block (lower 16 bits).
// See [NewVirtualOffset], [VirtualOffset.BlockOffset], and
// [VirtualOffset.WithinBlock].
//
// The package provides three primary types:
//
//   - [Reader] decompresses a BGZF stream sequentially. It implements
//     io.Reader and io.ByteReader.
//   - [Writer] compresses data into BGZF blocks. It implements io.WriteCloser
//     and can compress blocks in parallel (see NewParallelWriter).
//   - [IndexedReader] reads BGZF data with random access by virtual offset over
//     an io.ReadSeeker, backed by an LRU cache of decompressed blocks
//     (DefaultCacheSize blocks by default). With a [GZIndex] loaded, it also
//     supports Seek by uncompressed byte position.
package bgzf
