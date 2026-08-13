// Package sinks3 lets a varstore be written straight to S3 (and S3-compatible
// endpoints), so a conversion puts its output where it will be read instead of
// filling a local disk that somebody then has to upload and clear.
//
// Importing this package registers the "s3" scheme with varstore, so a blank
// import is enough:
//
//	import _ "github.com/compgenlab/cghts/varstore/sinks3"
//
//	cgkit vcf-toparquet --out s3://bucket/cohort input.vcf.gz
//
// It is a separate package for the same reason iosource/s3 is: it pulls in the
// AWS SDK, and a program that only ever writes local stores should not carry
// that. Credentials come from the standard AWS chain, and AWS_ENDPOINT_URL
// points at a non-AWS target -- the same settings the read side already uses,
// because it is the same client.
//
// WHAT THIS COSTS, and it is worth knowing before pointing a long conversion at
// it: a member is uploaded in parts and does not exist until its last part
// lands. A run that dies without abandoning its uploads leaves those parts
// invisible to a bucket listing and billed until something removes them. The
// writer abandons them on every failure it can see; set a lifecycle rule with
// AbortIncompleteMultipartUpload for the failures it cannot -- a killed
// process, a lost node.
package sinks3

import (
	"context"
	"fmt"
	"io"
	"strings"
	"sync"

	"github.com/compgenlab/cghts/iosource/s3"
	"github.com/compgenlab/cghts/varstore"
)

func init() {
	varstore.RegisterSink("s3", func(base string) (varstore.Sink, error) {
		return New(context.Background(), base)
	})
}

// Sink writes a store's members to one S3 prefix.
type Sink struct {
	ctx    context.Context
	client *s3.Client
	base   string

	// Writers in flight, so a failed conversion can abandon their uploads
	// rather than leave the parts behind.
	mu      sync.Mutex
	writing map[string]*s3.Writer
}

// New returns a Sink writing under an s3:// prefix.
func New(ctx context.Context, base string) (*Sink, error) {
	if !strings.HasPrefix(base, "s3://") {
		return nil, fmt.Errorf("sinks3: %s is not an s3:// locator", base)
	}
	c, err := s3.Shared(ctx)
	if err != nil {
		return nil, err
	}
	return &Sink{
		ctx: ctx, client: c, base: strings.TrimSuffix(base, "/"),
		writing: map[string]*s3.Writer{},
	}, nil
}

func (s *Sink) uri(name string) string { return s.base + "/" + name }

func (s *Sink) Create(name string) (io.WriteCloser, error) {
	w, err := s.client.Create(s.ctx, s.uri(name))
	if err != nil {
		return nil, err
	}
	s.mu.Lock()
	s.writing[name] = w
	s.mu.Unlock()
	return &tracked{Sink: s, name: name, w: w}, nil
}

func (s *Sink) Remove(name string) error { return s.client.Delete(s.ctx, s.uri(name)) }

func (s *Sink) Stat(name string) (int64, bool, error) {
	return s.client.Stat(s.ctx, s.uri(name))
}

func (s *Sink) Describe() string { return s.base }

// Abort abandons a member's upload.
//
// This is what makes the Sink an [varstore.Aborter], and why the writer asks
// for one: on an object store there is no half-written member to delete, only
// an upload in progress. Deleting would remove nothing and leave the parts.
func (s *Sink) Abort(name string) error {
	s.mu.Lock()
	w := s.writing[name]
	delete(s.writing, name)
	s.mu.Unlock()
	if w == nil {
		// Already finished, so there is a real object to remove instead.
		return s.Remove(name)
	}
	return w.Abort()
}

// tracked forgets a writer once it has completed, so a later abort removes the
// finished object rather than trying to abandon an upload that has ended.
type tracked struct {
	*Sink
	name string
	w    *s3.Writer
}

func (t *tracked) Write(p []byte) (int, error) { return t.w.Write(p) }

func (t *tracked) Close() error {
	err := t.w.Close()
	t.Sink.mu.Lock()
	delete(t.Sink.writing, t.name)
	t.Sink.mu.Unlock()
	return err
}

var (
	_ varstore.Sink    = (*Sink)(nil)
	_ varstore.Aborter = (*Sink)(nil)
)
