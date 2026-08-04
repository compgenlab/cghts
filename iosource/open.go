package iosource

import (
	"context"
	"fmt"
	"io"
	"sort"
	"strings"
	"sync"
)

// Opener builds a ByteSource for a locator of one scheme.
type Opener func(ctx context.Context, locator string) (ByteSource, error)

var (
	schemesMu sync.RWMutex
	schemes   = map[string]Opener{}
)

// Register associates an Opener with a URI scheme ("s3", "gs", …).
//
// Transports live outside this package when they carry heavy dependencies —
// an S3 client pulls in the AWS SDK, which not every consumer wants linked —
// so they register themselves and callers enable one with a blank import:
//
//	import _ "github.com/compgenlab/cghts/iosource/s3"
//
// Registering the same scheme twice panics, matching the standard library's
// treatment of duplicate driver registration: it means two transports disagree
// about who owns a scheme, which is a programming error, not a runtime one.
func Register(scheme string, open Opener) {
	schemesMu.Lock()
	defer schemesMu.Unlock()
	if _, dup := schemes[scheme]; dup {
		panic("iosource: scheme " + scheme + " registered twice")
	}
	schemes[scheme] = open
}

// Schemes lists the registered schemes, sorted. "file" and "http(s)" are always
// available and are not listed here.
func Schemes() []string {
	schemesMu.RLock()
	defer schemesMu.RUnlock()
	out := make([]string, 0, len(schemes))
	for s := range schemes {
		out = append(out, s)
	}
	sort.Strings(out)
	return out
}

// Open returns a random-access source for a locator, dispatching on its scheme.
//
// A plain path opens a local file, http:// and https:// use Range requests, and
// any other scheme must have been registered. This is what lets an indexed
// reader accept "the data is over there" without knowing where "there" is.
func Open(ctx context.Context, locator string) (ByteSource, error) {
	scheme := Scheme(locator)
	switch scheme {
	case "":
		return OpenFile(locator)
	case "http", "https":
		return NewHTTPRange(locator)
	}
	schemesMu.RLock()
	open, ok := schemes[scheme]
	schemesMu.RUnlock()
	if !ok {
		known := Schemes()
		if len(known) == 0 {
			return nil, fmt.Errorf("iosource: no transport registered for %q (import a transport package, e.g. _ \"github.com/compgenlab/cghts/iosource/s3\")", scheme+"://")
		}
		return nil, fmt.Errorf("iosource: no transport registered for %q (registered: %s)", scheme+"://", strings.Join(known, ", "))
	}
	return open(ctx, locator)
}

// Sibling opens locator+suffix using the transport the locator names, so index
// resolution works the same way whatever the data source is.
func Sibling(ctx context.Context) SiblingOpener {
	return func(locator, suffix string) (io.ReadCloser, error) {
		switch Scheme(locator) {
		case "":
			return FileSibling(locator, suffix)
		case "http", "https":
			return HTTPSibling(locator, suffix)
		}
		src, err := Open(ctx, locator+suffix)
		if err != nil {
			return nil, err
		}
		size, err := src.Size()
		if err != nil {
			src.Close()
			return nil, err
		}
		return readCloser{Reader: io.NewSectionReader(src, 0, size), Closer: src}, nil
	}
}

type readCloser struct {
	io.Reader
	io.Closer
}

// OpenReader opens a locator for sequential reading.
//
// The bytes are returned exactly as stored; nothing is decompressed. Use it
// where a whole small file is wanted -- an index, a sample list, a manifest --
// rather than the random access [Open] provides.
func OpenReader(ctx context.Context, locator string) (io.ReadCloser, error) {
	src, err := Open(ctx, locator)
	if err != nil {
		return nil, err
	}
	size, err := src.Size()
	if err != nil {
		src.Close()
		return nil, err
	}
	return readCloser{Reader: io.NewSectionReader(src, 0, size), Closer: src}, nil
}

// Scheme returns the URI scheme of a locator, or "" for a filesystem path.
//
// A Windows drive letter ("C:\...") is not a scheme, so a single character does
// not count.
//
// It is exported because a caller has to classify a locator exactly the way
// [Open] will dispatch it. Code that decides "is this a file to stat, or a URL
// to fetch" and disagrees with Open by one edge case does not fail loudly -- it
// opens the wrong thing, or reads a URL as a name.
func Scheme(locator string) string {
	i := strings.Index(locator, "://")
	if i <= 1 {
		return ""
	}
	return strings.ToLower(locator[:i])
}

// IsRemote reports whether a locator names something other than a local file.
func IsRemote(locator string) bool { return Scheme(locator) != "" }
