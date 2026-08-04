package varstore

import (
	"context"
	"fmt"
	"io"
	"os"

	"github.com/compgenlab/cghts/iosource"
)

// member is one file of a store — calls, sites or regions — held open for the
// life of the store.
//
// Holding it open rather than reopening per query matters for two reasons.
// Remotely, every open costs a request and a footer parse, and a store is
// queried many times. Locally it was simply wasted work: the previous code
// reopened the file and re-parsed the footer for each scan.
type member struct {
	// ra is an *io.SectionReader rather than the ByteSource itself: parquet-go
	// discovers a reader's length by asserting for `Size() int64`, and
	// ByteSource.Size returns (int64, error). SectionReader has the shape it
	// wants and costs nothing.
	ra   *io.SectionReader
	size int64
	name string // locator, for error messages
	src  iosource.ByteSource
}

// openMember opens a store member from any locator: a filesystem path, an
// http(s):// URL, or any registered scheme such as s3://.
func openMember(ctx context.Context, locator string) (*member, error) {
	src, err := iosource.Open(ctx, locator)
	if err != nil {
		return nil, err
	}
	size, err := src.Size()
	if err != nil {
		src.Close()
		return nil, fmt.Errorf("%s: %w", locator, err)
	}
	return &member{ra: io.NewSectionReader(src, 0, size), size: size, name: locator, src: src}, nil
}

func (m *member) Close() error {
	if m == nil || m.src == nil {
		return nil
	}
	return m.src.Close()
}

// memberExists reports whether a store member is present, for any locator kind.
//
// Absence is not an error: a store written with --no-callable has no regions
// member, and Classify's behaviour depends on telling that apart from a failure.
func memberExists(ctx context.Context, locator string) bool {
	if !iosource.IsRemote(locator) {
		_, err := os.Stat(locator)
		return err == nil
	}
	src, err := iosource.Open(ctx, locator)
	if err != nil {
		return false
	}
	src.Close()
	return true
}
