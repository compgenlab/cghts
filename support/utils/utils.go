package utils

import (
	"fmt"
	"io"
	"strings"
)

func TrimFloat(x float64, prec int) string {
	s := fmt.Sprintf("%.*f", prec, x)
	s = strings.TrimRight(s, "0")
	s = strings.TrimRight(s, ".")
	return s
}

type PositionTrackingReader struct {
	r io.Reader
	n int64
}

func NewPositionTrackingReader(r io.Reader) *PositionTrackingReader {
	return &PositionTrackingReader{r, 0}
}
func (c *PositionTrackingReader) Position() int64 {
	return c.n
}
func (c *PositionTrackingReader) Read(p []byte) (int, error) {
	k, err := c.r.Read(p)
	c.n += int64(k)
	return k, err
}

// Semaphore is a counting semaphore backed by a buffered channel. Acquire and
// Release must be balanced (each goroutine that Acquires must Release).
type Semaphore chan struct{}

func NewSemaphore(n int) Semaphore {
	return make(Semaphore, n)
}

func (s Semaphore) Acquire() {
	s <- struct{}{}
}

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
