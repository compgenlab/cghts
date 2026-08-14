package varstore

import (
	"context"
	"fmt"
	"io"
	"os"

	"github.com/compgenlab/cghts/iosource"
	"github.com/parquet-go/parquet-go"
	"sync"
)

// table is one file of a store — calls, sites or regions — held open for the
// life of the store.
//
// Holding it open rather than reopening per query matters for two reasons.
// Remotely, every open costs a request and a footer parse, and a store is
// queried many times. Locally it was simply wasted work: the previous code
// reopened the file and re-parsed the footer for each scan.
type table struct {
	// ra is an *io.SectionReader rather than the ByteSource itself: parquet-go
	// discovers a reader's length by asserting for `Size() int64`, and
	// ByteSource.Size returns (int64, error). SectionReader has the shape it
	// wants and costs nothing.
	ra   *io.SectionReader
	size int64
	name string // locator, for error messages
	src  iosource.ByteSource

	// once guards file and parseErr. A table is read concurrently once shards
	// are scanned in parallel, and the lazy parse below is the only mutable
	// state on the read path.
	once     sync.Once
	parseErr error

	// file is the parsed footer, kept because parsing it is the expensive half
	// of holding the table open and it does not change for the life of the
	// store. Populated lazily by parsed(); a store is opened for many queries
	// but a table it never reads should not pay for one.
	file *parquet.File
}

// parsed returns the table's parquet footer, parsing it once.
//
// Everything that reads a table wants this, and before it existed each caller
// parsed its own: ParquetVolume parsed the calls footer at open and discarded it,
// scanParquetPruned parsed one per scan and then handed the raw ReaderAt to
// NewGenericReader, which parsed a second. Remotely each parse is a request.
// parsed returns the table's footer, parsing it once.
//
// ONCE, AND SAFELY. This was a bare check-then-assign, which is a data race the
// moment two goroutines read the same table -- and shards exist precisely so
// that several can be read at the same time. Two racing parses would also both
// pay for the footer, which over object storage is a request each.
func (m *table) parsed() (*parquet.File, error) {
	m.once.Do(func() {
		m.file, m.parseErr = parquet.OpenFile(m.ra, m.size)
	})
	return m.file, m.parseErr
}

// openTable opens a store table from any locator: a filesystem path, an
// http(s):// URL, or any registered scheme such as s3://.
func openTable(ctx context.Context, locator string) (*table, error) {
	src, err := iosource.Open(ctx, locator)
	if err != nil {
		return nil, err
	}
	size, err := src.Size()
	if err != nil {
		src.Close()
		return nil, fmt.Errorf("%s: %w", locator, err)
	}
	return &table{ra: io.NewSectionReader(src, 0, size), size: size, name: locator, src: src}, nil
}

func (m *table) Close() error {
	if m == nil || m.src == nil {
		return nil
	}
	return m.src.Close()
}

// tableExists reports whether a store table is present, for any locator kind.
//
// Absence is not by itself an error, so this answers rather than failing.
// Whether a *particular* absence is legitimate is a separate question the
// manifest settles, since it records how many rows each table held.
//
// Note that the writer creates all three tables regardless, so --no-callable
// yields a present, zero-row regions file rather than a missing one -- the
// comments here used to say otherwise, and any logic keyed on that absence
// would never have fired.
func tableExists(ctx context.Context, locator string) bool {
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

// TableShape reports a store table's size in bytes and, if it has a readable
// parquet footer, its row count. A rows of -1 means the table is present but
// was never finalized.
//
// That distinction is the whole point. This exists for diagnosing a store that
// will not open, and since a missing manifest has no escape hatch, a caller left
// holding one needs some way to learn what is in it. A footer is written only by
// the writer's Close, so "present, N bytes, no footer" says precisely that the
// conversion died while writing this table -- which reporting it as simply
// absent would hide. An error means the table is not there at all.
//
// Row counts are footer metadata rather than a scan, so this is cheap, and it
// answers nothing about genotypes: diagnosis is deliberately not access.
func TableShape(ctx context.Context, locator string) (rows, size int64, err error) {
	m, err := openTable(ctx, locator)
	if err != nil {
		return 0, 0, err
	}
	defer m.Close()
	f, err := parquet.OpenFile(m.ra, m.size)
	if err != nil {
		return -1, m.size, nil
	}
	return f.NumRows(), m.size, nil
}
