package iosource

import (
	"fmt"
	"io"
	"net"
	"net/http"
	"strconv"
	"strings"
	"sync/atomic"
	"time"
)

// DefaultClient is the shared HTTP client for remote byte-source and reference
// access.
//
// It deliberately sets connection-level timeouts (dial, TLS handshake, and
// time-to-first-response-byte) rather than an overall Client.Timeout: genomic
// reads can be large and slow but legitimate, so we want to fail fast on a
// stalled or unresponsive server without capping total body-transfer time.
var DefaultClient = &http.Client{
	Transport: &http.Transport{
		Proxy: http.ProxyFromEnvironment,
		DialContext: (&net.Dialer{
			Timeout:   30 * time.Second,
			KeepAlive: 30 * time.Second,
		}).DialContext,
		ForceAttemptHTTP2:     true,
		MaxIdleConns:          100,
		IdleConnTimeout:       90 * time.Second,
		TLSHandshakeTimeout:   30 * time.Second,
		ExpectContinueTimeout: 5 * time.Second,
		ResponseHeaderTimeout: 60 * time.Second,
	},
}

// HTTPOption configures an [HTTPRange].
type HTTPOption func(*HTTPRange)

// WithHeader adds a request header (e.g. "Authorization") sent on every
// request the source makes.
func WithHeader(key, value string) HTTPOption {
	return func(h *HTTPRange) { h.header.Add(key, value) }
}

// WithClient overrides the HTTP client used by the source.
func WithClient(c *http.Client) HTTPOption {
	return func(h *HTTPRange) { h.client = c }
}

// HTTPRange is a [ByteSource] backed by an HTTP(S) server that supports Range
// requests. It is safe for concurrent use.
type HTTPRange struct {
	url    string
	client *http.Client
	header http.Header
	size   atomic.Int64 // total length, or -1 until learned
}

// NewHTTPRange opens url as a random-access byte source. It issues a HEAD to
// learn the content length up front; if the server does not report a usable
// length, the size is learned lazily from the Content-Range of the first
// ranged read (see [HTTPRange.Size]).
func NewHTTPRange(url string, opts ...HTTPOption) (*HTTPRange, error) {
	h := &HTTPRange{url: url, client: DefaultClient, header: make(http.Header)}
	h.size.Store(-1)
	for _, opt := range opts {
		opt(h)
	}

	// Best-effort size probe via HEAD; failures are tolerated because size can
	// still be recovered from the first ranged GET's Content-Range header.
	req, err := http.NewRequest(http.MethodHead, url, nil)
	if err != nil {
		return nil, err
	}
	h.applyHeader(req)
	if resp, err := h.client.Do(req); err == nil {
		resp.Body.Close()
		if resp.StatusCode == http.StatusOK && resp.ContentLength >= 0 {
			h.size.Store(resp.ContentLength)
		}
	}
	return h, nil
}

func (h *HTTPRange) applyHeader(req *http.Request) {
	for k, vs := range h.header {
		for _, v := range vs {
			req.Header.Add(k, v)
		}
	}
}

// ReadAt implements io.ReaderAt using an HTTP Range request. It fills p from
// offset off and returns io.EOF when the read reaches the end of the resource.
func (h *HTTPRange) ReadAt(p []byte, off int64) (int, error) {
	if len(p) == 0 {
		return 0, nil
	}
	end := off + int64(len(p)) - 1

	req, err := http.NewRequest(http.MethodGet, h.url, nil)
	if err != nil {
		return 0, err
	}
	h.applyHeader(req)
	req.Header.Set("Range", fmt.Sprintf("bytes=%d-%d", off, end))

	resp, err := h.client.Do(req)
	if err != nil {
		return 0, err
	}
	defer resp.Body.Close()

	switch resp.StatusCode {
	case http.StatusPartialContent:
		if h.size.Load() < 0 {
			if total := parseContentRangeTotal(resp.Header.Get("Content-Range")); total >= 0 {
				h.size.Store(total)
			}
		}
	case http.StatusOK:
		// Server ignored the Range header and returned the whole body; skip to
		// the requested offset so the source still behaves correctly (just less
		// efficiently). resp.ContentLength, if known, is the full size.
		if resp.ContentLength >= 0 {
			h.size.Store(resp.ContentLength)
		}
		if _, err := io.CopyN(io.Discard, resp.Body, off); err != nil {
			if err == io.EOF {
				return 0, io.EOF
			}
			return 0, err
		}
	case http.StatusRequestedRangeNotSatisfiable:
		return 0, io.EOF
	default:
		return 0, fmt.Errorf("http range %s: unexpected status %d", h.url, resp.StatusCode)
	}

	n, err := io.ReadFull(resp.Body, p)
	if err == io.ErrUnexpectedEOF || err == io.EOF {
		// Reached the end of the resource before filling p.
		return n, io.EOF
	}
	return n, err
}

// Size reports the total length of the resource, probing the server with a
// small ranged read if the length was not learned at open time.
func (h *HTTPRange) Size() (int64, error) {
	if s := h.size.Load(); s >= 0 {
		return s, nil
	}
	var b [1]byte
	if _, err := h.ReadAt(b[:], 0); err != nil && err != io.EOF {
		return 0, err
	}
	if s := h.size.Load(); s >= 0 {
		return s, nil
	}
	return 0, fmt.Errorf("http range %s: server did not report content length", h.url)
}

// Close releases resources. For an HTTP source this is a no-op.
func (h *HTTPRange) Close() error { return nil }

// parseContentRangeTotal extracts TOTAL from a "bytes START-END/TOTAL" header,
// returning -1 when the total is unknown ("*") or unparseable.
func parseContentRangeTotal(v string) int64 {
	i := strings.LastIndex(v, "/")
	if i < 0 {
		return -1
	}
	total := strings.TrimSpace(v[i+1:])
	if total == "" || total == "*" {
		return -1
	}
	n, err := strconv.ParseInt(total, 10, 64)
	if err != nil {
		return -1
	}
	return n
}
