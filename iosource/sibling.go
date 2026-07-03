package iosource

import (
	"fmt"
	"io"
	"net/http"
	"os"
)

// SiblingOpener opens an index "sibling" of a data resource — for example the
// ".bai" beside a ".bam", or the ".tbi"/".csi" beside a bgzipped VCF. It
// receives the data locator (path or URL) and the sibling suffix, and returns
// a reader over the sibling's bytes. Callers speaking a custom transport
// (SFTP, S3, ...) supply their own opener.
type SiblingOpener func(locator, suffix string) (io.ReadCloser, error)

// FileSibling opens locator+suffix from the local filesystem.
func FileSibling(locator, suffix string) (io.ReadCloser, error) {
	return os.Open(locator + suffix)
}

// HTTPSibling fetches locator+suffix over HTTP(S) using [DefaultClient].
func HTTPSibling(locator, suffix string) (io.ReadCloser, error) {
	resp, err := DefaultClient.Get(locator + suffix)
	if err != nil {
		return nil, err
	}
	if resp.StatusCode != http.StatusOK {
		resp.Body.Close()
		return nil, fmt.Errorf("http sibling %s%s: status %d", locator, suffix, resp.StatusCode)
	}
	return resp.Body, nil
}

// ResolveSibling tries each suffix in order and returns a reader for the first
// that opens, along with the matched suffix. It is used to locate an index
// whose exact extension is not known ahead of time (e.g. ".tbi" vs ".csi").
// The caller owns closing the returned reader.
func ResolveSibling(locator string, suffixes []string, open SiblingOpener) (io.ReadCloser, string, error) {
	var lastErr error
	for _, suffix := range suffixes {
		rc, err := open(locator, suffix)
		if err == nil {
			return rc, suffix, nil
		}
		lastErr = err
	}
	if lastErr == nil {
		lastErr = fmt.Errorf("no sibling found for %s (tried %v)", locator, suffixes)
	}
	return nil, "", lastErr
}
