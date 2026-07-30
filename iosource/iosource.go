// Package iosource provides pluggable, random-access byte sources for genomic
// data files.
//
// A [ByteSource] abstracts over local files and remote transports behind a
// concurrency-safe io.ReaderAt, so that index-driven readers can seek into a
// large file (BAM/CRAM/VCF/BigWig, ...) and fetch only the byte ranges they
// need without downloading the whole file. This package ships local-file and
// HTTP(S)-Range implementations, both using only the standard library. Callers
// that speak other transports (SFTP, S3, ...) implement [ByteSource] directly.
package iosource

import "io"

// ByteSource is a random-access, concurrency-safe source of bytes.
//
// ReadAt must satisfy the io.ReaderAt contract, including being safe for
// concurrent use: a header scan and one or more indexed region queries may
// share a single ByteSource. In particular a short read (n < len(p)) must be
// accompanied by a non-nil error, and io.EOF is returned when a read reaches
// the end of the source.
type ByteSource interface {
	io.ReaderAt

	// Size reports the total length of the source in bytes.
	Size() (int64, error)

	// Close releases any resources held by the source.
	Close() error
}

// ReadSeeker adapts a ByteSource to the io.ReadSeeker that block-compressed
// readers expect.
//
// It exists because the two halves of the abstraction meet awkwardly: a remote
// source is naturally random-access (ReadAt over byte ranges), while BGZF
// decompression is naturally sequential-with-seeks. io.SectionReader bridges
// exactly that, given a length — which is why ByteSource reports Size.
func ReadSeeker(src ByteSource) (io.ReadSeeker, error) {
	size, err := src.Size()
	if err != nil {
		return nil, err
	}
	return io.NewSectionReader(src, 0, size), nil
}
