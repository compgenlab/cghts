package htsio

// SamWriter is the interface for writing SAM/BAM/CRAM records.
type SamWriter interface {
	// Write writes a SamRecord to the output.
	Write(rec *SamRecord) error
	// Close flushes remaining data and releases resources.
	Close() error
}
