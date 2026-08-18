package varstore

import (
	"errors"
	"fmt"
	"io"
	"io/fs"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"sync"
)

// Where a store's tables are written.
//
// A store is four files that only mean anything together, and until now they
// could only be local ones. Nothing about the format required that: Parquet
// emits row groups forward and writes its footer last, so a table needs an
// io.Writer and never seeks. The filesystem was simply what the writer happened
// to be spelled in.
//
// Sink is that spelling made replaceable. The local implementation is the
// behaviour that was already here; a remote one is registered by importing the
// package that provides it, the same way iosource gains a read transport:
//
//	import _ "github.com/compgenlab/cghts/varstore/sinks3"
//
// which keeps the AWS SDK out of this package, and out of every program that
// only ever writes local stores.
//
// WHY THE SINK OWNS REMOVAL as well as creation: aborting a conversion has to
// undo what it made, and undoing differs by more than the call used. Locally a
// half-written table is a file to unlink. On an object store nothing is
// visible until the upload completes, so there is no half-written table to
// remove and abandoning is what has to happen instead. A caller cannot pick
// correctly between those, so it does not choose.

// Sink creates and removes the tables of one store.
//
// Names are table file names -- "calls.parquet", "manifest.json.gz" -- not
// paths. Where the store lives is the sink's business, which is what lets the
// writer be written once.
type Sink interface {
	// Create returns a writer for one table, replacing any previous content.
	Create(name string) (io.WriteCloser, error)

	// Remove deletes a table. A table that is not there is not an error.
	Remove(name string) error

	// Stat reports a table's size, and whether it exists.
	Stat(name string) (size int64, ok bool, err error)

	// Describe names the destination for an error message.
	Describe() string
}

// Aborter is a Sink whose tables are not visible until they are finished, and
// which therefore abandons rather than deletes.
//
// Implemented by object-store sinks. The writer prefers it when tearing down a
// failed conversion: calling Remove on an upload that was never completed would
// delete nothing and leave the parts behind.
type Aborter interface {
	Abort(name string) error
}

var (
	sinksMu sync.RWMutex
	sinks   = map[string]func(base string) (Sink, error){}
)

// RegisterSink associates a scheme with a way of writing stores to it.
//
// Registering the same scheme twice panics, matching iosource and the standard
// library's treatment of duplicate driver registration: two writers disagreeing
// about who owns a scheme is a programming error, not a runtime one.
func RegisterSink(scheme string, open func(base string) (Sink, error)) {
	sinksMu.Lock()
	defer sinksMu.Unlock()
	if _, dup := sinks[scheme]; dup {
		panic("varstore: sink for scheme " + scheme + " registered twice")
	}
	sinks[scheme] = open
}

// SinkSchemes lists the registered remote schemes, sorted. Local paths are
// always writable and are not listed.
func SinkSchemes() []string {
	sinksMu.RLock()
	defer sinksMu.RUnlock()
	out := make([]string, 0, len(sinks))
	for s := range sinks {
		out = append(out, s)
	}
	sort.Strings(out)
	return out
}

// CanWrite reports whether a store can be written to base.
//
// It exists so a command can refuse a locator with a message naming what is
// missing -- an unregistered scheme otherwise looks like a directory with a
// very strange name.
func CanWrite(base string) bool {
	scheme := schemeOf(base)
	if scheme == "" {
		return true
	}
	sinksMu.RLock()
	defer sinksMu.RUnlock()
	_, ok := sinks[scheme]
	return ok
}

// OpenSink returns the sink for a store base, local or remote.
func OpenSink(base string) (Sink, error) {
	scheme := schemeOf(base)
	if scheme == "" {
		return fileSink{dir: trimStoreDir(base)}, nil
	}
	sinksMu.RLock()
	open, ok := sinks[scheme]
	sinksMu.RUnlock()
	if !ok {
		known := SinkSchemes()
		if len(known) == 0 {
			return nil, fmt.Errorf(
				"varstore: cannot write to %s: no transport for %q is linked into this program", base, scheme)
		}
		return nil, fmt.Errorf(
			"varstore: cannot write to %s: no transport for %q (registered: %s)",
			base, scheme, strings.Join(known, ", "))
	}
	return open(trimStoreDir(base))
}

// schemeOf returns the URI scheme of a locator, or "" for a filesystem path.
//
// A Windows drive letter is a path, not a scheme: "C:\store" has a colon in
// the same place a scheme does, so the separator that follows is what tells
// them apart.
func schemeOf(locator string) string {
	i := strings.Index(locator, "://")
	if i <= 0 {
		return ""
	}
	return locator[:i]
}

// fileSink writes a store into a directory, which is what every store was
// before this interface existed.
type fileSink struct{ dir string }

func (f fileSink) Create(name string) (io.WriteCloser, error) {
	if err := os.MkdirAll(f.dir, 0o755); err != nil {
		return nil, err
	}
	// A split table's name carries a directory -- "calls/00007.parquet" -- and
	// on an object store that is just a key with a slash in it. On a filesystem
	// it is a directory that has to exist first, which is the one place the two
	// kinds of sink genuinely differ.
	if dir := filepath.Dir(filepath.Join(f.dir, name)); dir != f.dir {
		if err := os.MkdirAll(dir, 0o755); err != nil {
			return nil, err
		}
	}
	return os.Create(joinStore(f.dir, name))
}

func (f fileSink) Remove(name string) error {
	err := os.Remove(joinStore(f.dir, name))
	if err != nil && !errors.Is(err, fs.ErrNotExist) {
		return err
	}
	return nil
}

func (f fileSink) Stat(name string) (int64, bool, error) {
	st, err := os.Stat(joinStore(f.dir, name))
	switch {
	case errors.Is(err, fs.ErrNotExist):
		return 0, false, nil
	case err != nil:
		return 0, false, err
	}
	return st.Size(), true, nil
}

func (f fileSink) Describe() string { return f.dir }
