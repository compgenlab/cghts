package s3

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"io"

	"github.com/aws/aws-sdk-go-v2/aws"
	awss3 "github.com/aws/aws-sdk-go-v2/service/s3"
	"github.com/aws/aws-sdk-go-v2/service/s3/types"
	"github.com/aws/smithy-go"
)

// Writing objects, so a producer can put its output where it sits rather than
// writing a local file somebody then has to upload.
//
// This is the other half of what this package already does for reading, and it
// exists because the formats being written do not need it to be hard: Parquet
// emits row groups forward and puts its footer last, so a writer needs an
// io.Writer and never seeks. The same is true of gzip, of BGZF, and of every
// other stream this library produces.
//
// MULTIPART ONLY WHEN IT HAS TO BE. Nothing is sent until PartSize has
// accumulated, so an object smaller than that -- a manifest, a small index --
// goes out as a single PutObject and never creates a multipart upload at all.
// That matters beyond tidiness: an abandoned multipart upload is invisible to a
// bucket listing and is billed until something removes it, so the cheapest way
// to handle one is not to start it.

// PartSize is how much is buffered before a part is sent.
//
// S3 requires every part except the last to be at least 5 MiB; 8 MiB leaves
// room above that floor without holding much. The 10,000-part limit puts the
// largest writable object at 80 GB, and a caller writing more than that should
// say so with WithPartSize rather than have it chosen for them.
const PartSize = 8 << 20

// maxParts is S3's limit. Exceeding it fails at CompleteMultipartUpload, after
// everything has been transferred, which is the worst possible moment to find
// out -- so it is checked as the parts are made.
const maxParts = 10_000

// Writer streams an object into S3.
//
// Not safe for concurrent use: it is a stream, and the order bytes arrive in is
// the order they are stored.
type Writer struct {
	ctx    context.Context
	client *awss3.Client
	bucket string
	key    string

	partSize int
	buf      []byte

	uploadID string
	parts    []types.CompletedPart
	written  int64
	closed   bool
}

// WriterOption adjusts a Writer before it is used.
type WriterOption func(*Writer)

// WithPartSize sets the buffer that decides when a part is sent. Values below
// S3's 5 MiB floor are raised to it, since a smaller part is refused on upload
// rather than at the point it was configured.
func WithPartSize(n int) WriterOption {
	return func(w *Writer) {
		if n < 5<<20 {
			n = 5 << 20
		}
		w.partSize = n
	}
}

// Create returns a Writer for an s3:// URI.
//
// No request is made here. The object comes into existence when the Writer is
// closed -- or, once it has grown past a part, when the multipart upload
// completes -- so a Writer that is created and abandoned leaves nothing behind.
func (c *Client) Create(ctx context.Context, uri string, opts ...WriterOption) (*Writer, error) {
	bucket, key, err := ParseURI(uri)
	if err != nil {
		return nil, err
	}
	if key == "" {
		return nil, fmt.Errorf("s3: %s names a bucket, not an object", uri)
	}
	w := &Writer{
		ctx: ctx, client: c.api, bucket: bucket, key: key,
		partSize: PartSize,
	}
	for _, o := range opts {
		o(w)
	}
	w.buf = make([]byte, 0, w.partSize)
	return w, nil
}

func (w *Writer) Write(p []byte) (int, error) {
	if w.closed {
		return 0, errors.New("s3: write to a closed object")
	}
	n := len(p)
	w.written += int64(n)
	for len(p) > 0 {
		room := w.partSize - len(w.buf)
		if room > len(p) {
			room = len(p)
		}
		w.buf = append(w.buf, p[:room]...)
		p = p[room:]
		if len(w.buf) >= w.partSize {
			if err := w.flush(); err != nil {
				return n - len(p), err
			}
		}
	}
	return n, nil
}

// flush sends the buffer as a part, starting the multipart upload if this is
// the first one.
func (w *Writer) flush() error {
	if len(w.buf) == 0 {
		return nil
	}
	if w.uploadID == "" {
		out, err := w.client.CreateMultipartUpload(w.ctx, &awss3.CreateMultipartUploadInput{
			Bucket: aws.String(w.bucket), Key: aws.String(w.key),
		})
		if err != nil {
			return fmt.Errorf("s3: starting an upload of %s: %w", w.uri(), err)
		}
		w.uploadID = aws.ToString(out.UploadId)
	}
	if len(w.parts) >= maxParts {
		return fmt.Errorf("s3: %s needs more than %d parts at %d bytes each; use WithPartSize",
			w.uri(), maxParts, w.partSize)
	}

	n := int32(len(w.parts) + 1)
	out, err := w.client.UploadPart(w.ctx, &awss3.UploadPartInput{
		Bucket: aws.String(w.bucket), Key: aws.String(w.key),
		UploadId: aws.String(w.uploadID), PartNumber: aws.Int32(n),
		Body: bytes.NewReader(w.buf),
	})
	if err != nil {
		return fmt.Errorf("s3: uploading part %d of %s: %w", n, w.uri(), err)
	}
	w.parts = append(w.parts, types.CompletedPart{ETag: out.ETag, PartNumber: aws.Int32(n)})
	w.buf = w.buf[:0]
	return nil
}

// Close finishes the object, which is the point at which it becomes visible.
//
// An object that never grew past one part is sent here as a single PutObject,
// so the common small case costs one request and leaves no multipart upload to
// clean up.
func (w *Writer) Close() error {
	if w.closed {
		return nil
	}
	w.closed = true

	if w.uploadID == "" {
		_, err := w.client.PutObject(w.ctx, &awss3.PutObjectInput{
			Bucket: aws.String(w.bucket), Key: aws.String(w.key),
			Body: bytes.NewReader(w.buf),
		})
		if err != nil {
			return fmt.Errorf("s3: writing %s: %w", w.uri(), err)
		}
		w.buf = nil
		return nil
	}

	if err := w.flush(); err != nil {
		return errors.Join(err, w.abortLocked())
	}
	_, err := w.client.CompleteMultipartUpload(w.ctx, &awss3.CompleteMultipartUploadInput{
		Bucket: aws.String(w.bucket), Key: aws.String(w.key),
		UploadId:        aws.String(w.uploadID),
		MultipartUpload: &types.CompletedMultipartUpload{Parts: w.parts},
	})
	if err != nil {
		return errors.Join(fmt.Errorf("s3: completing %s: %w", w.uri(), err), w.abortLocked())
	}
	return nil
}

// Abort abandons the object without creating it.
//
// Worth calling on any failure path. Parts already uploaded do not appear in a
// bucket listing and are billed until they are removed, so an upload nobody
// aborts is storage nobody can see. A bucket lifecycle rule with
// AbortIncompleteMultipartUpload is the backstop for the case where the process
// dies before it can call this.
func (w *Writer) Abort() error {
	w.closed = true
	return w.abortLocked()
}

func (w *Writer) abortLocked() error {
	if w.uploadID == "" {
		return nil
	}
	_, err := w.client.AbortMultipartUpload(w.ctx, &awss3.AbortMultipartUploadInput{
		Bucket: aws.String(w.bucket), Key: aws.String(w.key),
		UploadId: aws.String(w.uploadID),
	})
	w.uploadID = ""
	if err != nil {
		return fmt.Errorf("s3: abandoning %s: %w", w.uri(), err)
	}
	return nil
}

// Written reports how many bytes have been handed to this Writer.
func (w *Writer) Written() int64 { return w.written }

func (w *Writer) uri() string { return "s3://" + w.bucket + "/" + w.key }

// Delete removes an object. A key that is not there is not an error: the
// caller wanted it gone.
func (c *Client) Delete(ctx context.Context, uri string) error {
	bucket, key, err := ParseURI(uri)
	if err != nil {
		return err
	}
	if _, err := c.api.DeleteObject(ctx, &awss3.DeleteObjectInput{
		Bucket: aws.String(bucket), Key: aws.String(key),
	}); err != nil && !isNotFound(err) {
		return fmt.Errorf("s3: removing %s: %w", uri, err)
	}
	return nil
}

// Stat returns an object's size, and whether it is there at all.
func (c *Client) Stat(ctx context.Context, uri string) (size int64, ok bool, err error) {
	bucket, key, err := ParseURI(uri)
	if err != nil {
		return 0, false, err
	}
	out, err := c.api.HeadObject(ctx, &awss3.HeadObjectInput{
		Bucket: aws.String(bucket), Key: aws.String(key),
	})
	if err != nil {
		if isNotFound(err) {
			return 0, false, nil
		}
		return 0, false, fmt.Errorf("s3: reading %s: %w", uri, err)
	}
	return aws.ToInt64(out.ContentLength), true, nil
}

func isNotFound(err error) bool {
	var nf *types.NotFound
	if errors.As(err, &nf) {
		return true
	}
	var nsk *types.NoSuchKey
	if errors.As(err, &nsk) {
		return true
	}
	// A HeadObject against a missing key answers 404 with no body, so the
	// typed error above is not always what arrives.
	var api smithy.APIError
	if errors.As(err, &api) {
		switch api.ErrorCode() {
		case "NotFound", "NoSuchKey", "404":
			return true
		}
	}
	return false
}

var _ io.WriteCloser = (*Writer)(nil)
