// Package bbi reads local UCSC BBI files — bigWig and bigBed — with random
// access by genomic region.
//
// A BBI file is self-indexed: a chromosome B+ tree maps reference names to ids
// and an R-tree spatially indexes the compressed data blocks, so no sidecar
// index (unlike tabix) is needed. Only base-resolution data is read; the
// display zoom-level summaries are ignored, so queried values are exact.
//
// The reader surface mirrors htsio/tabix.Reader:
//
//	r, err := bbi.Open("scores.bw")
//	seq, err := r.Query("chr1", 999, 1000) // 0-based half-open
//	for rec, err := range seq { ... }       // rec.Value (bigWig) or rec.Line (bigBed)
//
// Byte order is detected from the file magic (BBI files may be little- or
// big-endian). The package depends only on the standard library.
package bbi
