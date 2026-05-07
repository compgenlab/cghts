// Package codec provides general-purpose byte-level codecs used by CRAM
// and potentially other HTS formats.
//
// Each codec supports both a simple []byte API and streaming io.Reader/io.Writer
// interfaces. The streaming interfaces read/write complete compressed blocks —
// they do not support incremental compression across multiple blocks.
//
// Decoder usage:
//
//	decoder := codec.NewRans4x8Decoder(compressedReader)
//	decoded, err := io.ReadAll(decoder)
//
// Encoder usage:
//
//	encoder := codec.NewRans4x8Encoder(outputWriter, codec.Order0)
//	encoder.Write(data)
//	encoder.Close() // compresses and flushes to outputWriter
package codec

import (
	"bytes"
	"fmt"
	"io"
)

// Codec method constants matching CRAM block method IDs.
const (
	MethodRans4x8  = 4
	MethodRansNx16 = 5
	MethodArith    = 6
	MethodFqzcomp  = 7
	MethodNameTok  = 8
)

// Order constants for rANS encoding.
const (
	Order0 = 0
	Order1 = 1
)

// Decoder reads compressed data from an io.Reader, decodes a complete block,
// and provides the decoded bytes via the Read method.
type Decoder struct {
	r       io.Reader
	decoded []byte
	pos     int
	err     error
	method  int
}

// NewRans4x8Decoder creates a decoder that reads rANS 4x8 compressed data.
func NewRans4x8Decoder(r io.Reader) *Decoder {
	return &Decoder{r: r, method: MethodRans4x8}
}

// NewRansNx16Decoder creates a decoder that reads rANS Nx16 compressed data.
func NewRansNx16Decoder(r io.Reader) *Decoder {
	return &Decoder{r: r, method: MethodRansNx16}
}

// NewArithDecoder creates a decoder that reads adaptive arithmetic coded data.
func NewArithDecoder(r io.Reader) *Decoder {
	return &Decoder{r: r, method: MethodArith}
}

func (d *Decoder) init() {
	if d.decoded != nil || d.err != nil {
		return
	}
	data, err := io.ReadAll(d.r)
	if err != nil {
		d.err = err
		return
	}
	switch d.method {
	case MethodRans4x8:
		d.decoded, d.err = DecodeRans4x8(data)
	case MethodRansNx16:
		d.decoded, d.err = DecodeRansNx16(data)
	case MethodArith:
		d.decoded, d.err = DecodeArithDynamic(data)
	default:
		d.err = fmt.Errorf("codec: unsupported method %d", d.method)
	}
}

// Read implements io.Reader. On the first call, it reads all compressed data
// from the underlying reader and decodes the block.
func (d *Decoder) Read(p []byte) (int, error) {
	d.init()
	if d.err != nil {
		return 0, d.err
	}
	if d.pos >= len(d.decoded) {
		return 0, io.EOF
	}
	n := copy(p, d.decoded[d.pos:])
	d.pos += n
	return n, nil
}

// Encoder collects data via Write, then compresses and writes to the
// underlying io.Writer on Close.
type Encoder struct {
	w      io.Writer
	buf    bytes.Buffer
	method int
	order  int
	closed bool
}

// NewRans4x8Encoder creates an encoder that writes rANS 4x8 compressed data.
// order is Order0 or Order1.
func NewRans4x8Encoder(w io.Writer, order int) *Encoder {
	return &Encoder{w: w, method: MethodRans4x8, order: order}
}

// NewRansNx16Encoder creates an encoder that writes rANS Nx16 compressed data.
func NewRansNx16Encoder(w io.Writer) *Encoder {
	return &Encoder{w: w, method: MethodRansNx16}
}

// Write buffers data for compression. The actual encoding happens on Close.
func (e *Encoder) Write(p []byte) (int, error) {
	if e.closed {
		return 0, fmt.Errorf("codec: write after close")
	}
	return e.buf.Write(p)
}

// Close compresses the buffered data and writes the compressed block to
// the underlying writer.
func (e *Encoder) Close() error {
	if e.closed {
		return nil
	}
	e.closed = true

	data := e.buf.Bytes()
	var compressed []byte

	switch e.method {
	case MethodRans4x8:
		compressed = EncodeRans4x8(data, e.order)
	case MethodRansNx16:
		compressed = EncodeRansNx16(data)
	default:
		return fmt.Errorf("codec: unsupported encode method %d", e.method)
	}

	_, err := e.w.Write(compressed)
	return err
}
