package seqio

import (
	"net"
	"net/http"
	"time"
)

// httpClient is the shared HTTP client for all remote reference access
// (remote FASTA, refget, and the on-disk cache fetches).
//
// It deliberately sets connection-level timeouts (dial, TLS handshake, and
// time-to-first-response-byte) rather than an overall Client.Timeout: reference
// downloads can be large and slow but legitimate, so we want to fail fast on a
// stalled or unresponsive server without capping the total body-transfer time.
var httpClient = &http.Client{
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
