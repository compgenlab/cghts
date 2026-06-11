package htsio

// SamWriter is the interface for writing SAM/BAM/CRAM records.
//
// Implementations are safe for concurrent use, but records are written in call
// order, so a single producer is the normal pattern. After Close, Write may
// return an error and must not be relied upon.
type SamWriter interface {
	// Write writes a SamRecord to the output. It does not take ownership of rec.
	Write(rec *SamRecord) error
	// Close flushes any buffered records and writes the format trailer (e.g. the
	// BGZF/CRAM EOF marker), then releases resources. It must be called exactly
	// once after the final Write to produce a valid file; the returned error
	// reports any failure from that final flush. Close is idempotent: a second
	// call is a no-op that returns the first call's result.
	Close() error
}
