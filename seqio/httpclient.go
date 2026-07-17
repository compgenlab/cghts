package seqio

import "github.com/compgenlab/cghts/iosource"

// httpClient is the shared HTTP client for all remote reference access
// (remote FASTA, refget, and the on-disk cache fetches). It is defined in the
// iosource package so remote byte sources and reference readers share one
// tuned transport; see [iosource.DefaultClient] for the timeout rationale.
var httpClient = iosource.DefaultClient
