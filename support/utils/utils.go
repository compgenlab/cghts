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

func (s Semaphore) Close() {
	close(s)
}
