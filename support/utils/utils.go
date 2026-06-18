package utils

import (
	"fmt"
	"io"
	"strings"
)

// TrimFloat formats x with prec digits after the decimal point and then
// removes any trailing zeros, along with a trailing decimal point if no
// fractional digits remain. For example TrimFloat(1.2500, 4) returns "1.25"
// and TrimFloat(3.0, 2) returns "3".
func TrimFloat(x float64, prec int) string {
	s := fmt.Sprintf("%.*f", prec, x)
	s = strings.TrimRight(s, "0")
	s = strings.TrimRight(s, ".")
	return s
}

// PositionTrackingReader wraps an [io.Reader] and counts the total number of
// bytes that have been read through it. This is useful for reporting progress
// or recording a byte offset in a stream that does not support seeking.
type PositionTrackingReader struct {
	r io.Reader
	n int64
}

// NewPositionTrackingReader returns a [PositionTrackingReader] that reads from
// r and starts its byte count at zero.
func NewPositionTrackingReader(r io.Reader) *PositionTrackingReader {
	return &PositionTrackingReader{r, 0}
}

// Position returns the total number of bytes read so far through the reader.
func (c *PositionTrackingReader) Position() int64 {
	return c.n
}

// Read reads from the underlying reader, adds the number of bytes read to the
// running total, and returns the same values as the underlying reader. It
// implements [io.Reader].
func (c *PositionTrackingReader) Read(p []byte) (int, error) {
	k, err := c.r.Read(p)
	c.n += int64(k)
	return k, err
}

// Semaphore is a counting semaphore backed by a buffered channel. Acquire and
// Release must be balanced (each goroutine that Acquires must Release).
type Semaphore chan struct{}

// NewSemaphore returns a [Semaphore] that allows up to n concurrent holders.
// At most n calls to Acquire can succeed before a matching Release is made;
// further callers block until a slot is freed.
func NewSemaphore(n int) Semaphore {
	return make(Semaphore, n)
}

// Acquire takes one slot from the semaphore, blocking until a slot is available
// when the semaphore is already full.
func (s Semaphore) Acquire() {
	s <- struct{}{}
}

// Release returns one slot to the semaphore, allowing a goroutine blocked in
// Acquire to proceed. Each Acquire must be balanced by exactly one Release.
func (s Semaphore) Release() {
	<-s
}

// Close releases the semaphore's resources. It must only be called once all
// goroutines have finished using it — closing while a goroutine is blocked in
// Acquire would panic with "send on closed channel". Close is optional; a
// semaphore that is simply abandoned is reclaimed by the garbage collector.
func (s Semaphore) Close() {
	close(s)
}
