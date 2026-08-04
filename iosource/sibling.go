package iosource

import (
	"fmt"
	"io"
	"net/http"
	"os"
	"strings"
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
// The error names every suffix tried, not just the last failure. Returning only
// the last one meant a missing index reported a 404 for ".csi" and never
// mentioned ".tbi" -- which reads as though the wrong index kind was expected,
// when in fact neither was there. That is the most common failure for a remote
// file, so it is worth spelling out.
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
		return nil, "", fmt.Errorf("no sibling found for %s (tried %s)",
			locator, strings.Join(suffixes, ", "))
	}
	return nil, "", fmt.Errorf("no sibling found for %s (tried %s; last error: %w)",
		locator, strings.Join(suffixes, ", "), lastErr)
}
