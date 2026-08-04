package varstore

import (
	"context"
	"fmt"
	"io"
	"os"

	"github.com/compgenlab/cghts/iosource"
	"github.com/parquet-go/parquet-go"
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

// MemberShape reports a store member's row count and size in bytes, reading
// only its parquet footer.
//
// It exists for diagnosing a store that will not open. Since a missing manifest
// has no escape hatch, a caller holding an unreadable store needs some way to
// learn what is in it, and row counts are footer metadata rather than a scan --
// so this answers cheaply without answering any genotype question. Diagnosis is
// deliberately not access.
func MemberShape(ctx context.Context, locator string) (rows, size int64, err error) {
	m, err := openMember(ctx, locator)
	if err != nil {
		return 0, 0, err
	}
	defer m.Close()
	f, err := parquet.OpenFile(m.ra, m.size)
	if err != nil {
		return 0, 0, err
	}
	return f.NumRows(), m.size, nil
}
