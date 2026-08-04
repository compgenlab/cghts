package vcf

import (
	"context"

	"github.com/compgenlab/cghts/htsio/tabix"
	"github.com/compgenlab/cghts/iosource"
)

// OpenVcfFile opens a VCF for streaming reads from any locator: a filesystem
// path, an http(s):// URL, or any scheme registered with iosource such as
// s3://. gzip and BGZF are detected from the magic bytes, as in [NewVcfFile].
//
// A plain path is handed to NewVcfFile unchanged, so local behaviour is exactly
// what it always was.
//
// Note what a streaming read of a remote object costs: the whole file crosses
// the wire, because a stream has no index to skip with. Use
// [OpenIndexedVcfReader] when the question is about a region.
func OpenVcfFile(ctx context.Context, locator string) (*VcfReader, error) {
	if !iosource.IsRemote(locator) {
		return NewVcfFile(locator)
	}
	rc, err := iosource.OpenReader(ctx, locator)
	if err != nil {
		return nil, err
	}
	return newVcfFrom(rc, rc)
}

// OpenIndexedVcfReader opens a tabix-indexed VCF for random access from any
// locator, resolving the .tbi or .csi over the same transport as the data.
//
// A plain path is handed to [NewIndexedVcfReader] unchanged.
func OpenIndexedVcfReader(ctx context.Context, locator string) (*IndexedVcfReader, error) {
	if !iosource.IsRemote(locator) {
		return NewIndexedVcfReader(locator)
	}
	tr, err := tabix.Open(ctx, locator)
	if err != nil {
		return nil, err
	}
	// Ownership of tr transfers; the returned reader's Close releases it.
	return NewIndexedVcfReaderFrom(tr, locator), nil
}
