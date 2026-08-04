package tabix

import (
	"context"
	"fmt"

	"github.com/compgenlab/cghts/iosource"
)

// IndexSuffixes are the sidecar index names tried, in order.
var IndexSuffixes = []string{".tbi", ".csi"}

// Open opens a tabix-indexed file from any locator: a filesystem path, an
// http(s):// URL, or any scheme registered with iosource such as s3://. The
// index sidecar is resolved over the same transport as the data, so a remote
// file finds its remote index.
//
// A plain path is handed to [NewReaderSize] unchanged rather than routed through
// iosource. That is deliberate: the local constructor is the behaviour every
// existing caller already has, and re-deriving it here would be a second
// implementation free to drift from the first.
func Open(ctx context.Context, locator string, opts ...ReaderOption) (*Reader, error) {
	if !iosource.IsRemote(locator) {
		o := readerOpts{cacheWindows: DefaultRecordCacheWindows}
		for _, opt := range opts {
			opt(&o)
		}
		return NewReaderSize(locator, o.cacheWindows)
	}

	src, err := iosource.Open(ctx, locator)
	if err != nil {
		return nil, err
	}
	data, err := iosource.ReadSeeker(src)
	if err != nil {
		src.Close()
		return nil, err
	}
	idx, _, err := iosource.ResolveSibling(locator, IndexSuffixes, iosource.Sibling(ctx))
	if err != nil {
		src.Close()
		return nil, fmt.Errorf("tabix: %s needs an index: %w", locator, err)
	}
	defer idx.Close() // loadIndexFrom consumes it entirely

	// The source is handed to WithCloser so closing the Reader releases the
	// transport; nothing else holds it.
	r, err := NewReaderFromSource(data, idx, append(opts, WithCloser(src))...)
	if err != nil {
		src.Close()
		return nil, err
	}
	return r, nil
}
