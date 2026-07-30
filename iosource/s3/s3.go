// Package s3 reads objects from S3 (and S3-compatible endpoints) as
// [iosource.ByteSource]s, so any cghts reader that accepts a source can query a
// file where it sits instead of downloading it.
//
// Importing this package registers the "s3" scheme with iosource, so a blank
// import is enough to make s3:// locators work everywhere:
//
//	import _ "github.com/compgenlab/cghts/iosource/s3"
//
//	src, err := iosource.Open(ctx, "s3://bucket/clinvar.vcf.gz")
//
// Credentials come from the standard AWS chain — environment, shared config,
// then an instance or container role — so a deployed pod can use an assumed
// role rather than static keys. AWS_ENDPOINT_URL (or AWS_ENDPOINT_URL_S3)
// points at a non-AWS, S3-compatible target and switches on path-style
// addressing, which a local gateway needs because it has no wildcard DNS.
package s3

import (
	"context"
	"errors"
	"fmt"
	"io"
	"os"
	"strings"
	"sync"

	"github.com/aws/aws-sdk-go-v2/aws"
	awsconfig "github.com/aws/aws-sdk-go-v2/config"
	awss3 "github.com/aws/aws-sdk-go-v2/service/s3"
	"github.com/compgenlab/cghts/iosource"
)

func init() {
	iosource.Register("s3", func(ctx context.Context, locator string) (iosource.ByteSource, error) {
		c, err := Shared(ctx)
		if err != nil {
			return nil, err
		}
		return c.Open(ctx, locator)
	})
}

// Client reads from an S3-compatible endpoint.
type Client struct{ api *awss3.Client }

// Option configures a Client.
type Option func(*settings)

type settings struct {
	profile  string
	region   string
	endpoint string
	apply    []func(*awss3.Options)
}

// WithProfile selects a named profile from the shared config files, the same
// one `aws --profile` would use.
//
// Without it the SDK honours AWS_PROFILE, falling back to "default". This
// exists so a program can offer a flag or a config setting rather than making
// callers set an environment variable.
func WithProfile(name string) Option {
	return func(s *settings) { s.profile = name }
}

// WithRegion overrides the region, which otherwise comes from AWS_REGION,
// AWS_DEFAULT_REGION or the selected profile.
func WithRegion(region string) Option {
	return func(s *settings) { s.region = region }
}

// WithEndpoint points at a non-AWS, S3-compatible endpoint, overriding
// AWS_ENDPOINT_URL. Setting one also switches on path-style addressing.
func WithEndpoint(url string) Option {
	return func(s *settings) { s.endpoint = url }
}

// WithHTTPClient overrides the HTTP client. Mainly for tests that need to
// observe what actually went over the wire.
func WithHTTPClient(h aws.HTTPClient) Option {
	return func(s *settings) {
		s.apply = append(s.apply, func(o *awss3.Options) { o.HTTPClient = h })
	}
}

// New builds a client from the ambient AWS configuration.
//
// Credentials come from the standard chain, in order: environment variables;
// the shared credentials and config files (~/.aws/credentials, ~/.aws/config,
// honouring AWS_PROFILE and relocatable with AWS_SHARED_CREDENTIALS_FILE /
// AWS_CONFIG_FILE); then a container or EC2 instance role. Profile forms the
// SDK resolves — SSO, role_arn with source_profile, credential_process — all
// work, because none of that is reimplemented here.
func New(ctx context.Context, opts ...Option) (*Client, error) {
	var set settings
	for _, opt := range opts {
		opt(&set)
	}

	loadOpts := []func(*awsconfig.LoadOptions) error{}
	if set.profile != "" {
		loadOpts = append(loadOpts, awsconfig.WithSharedConfigProfile(set.profile))
	}
	if set.region != "" {
		loadOpts = append(loadOpts, awsconfig.WithRegion(set.region))
	}
	cfg, err := awsconfig.LoadDefaultConfig(ctx, loadOpts...)
	if err != nil {
		if set.profile != "" {
			return nil, fmt.Errorf("s3: aws config for profile %q: %w", set.profile, err)
		}
		return nil, fmt.Errorf("s3: aws config: %w", err)
	}
	if cfg.Region == "" {
		// A region is required to sign. Any value works against a compatible
		// gateway, and us-east-1 is the conventional default.
		cfg.Region = "us-east-1"
	}
	endpoint := set.endpoint
	if endpoint == "" {
		endpoint = firstEnv("AWS_ENDPOINT_URL_S3", "AWS_ENDPOINT_URL")
	}
	api := awss3.NewFromConfig(cfg, func(o *awss3.Options) {
		if endpoint != "" {
			o.BaseEndpoint = aws.String(endpoint)
			o.UsePathStyle = true
		}
		for _, f := range set.apply {
			f(o)
		}
	})
	return &Client{api: api}, nil
}

var (
	sharedOnce sync.Once
	sharedC    *Client
	sharedErr  error
)

// Shared returns a process-wide client, built on first use.
//
// Lazy so that merely importing this package — which is how the scheme gets
// registered — never touches the credential chain. A program that only ever
// reads local files must not fail because it has no AWS configuration.
func Shared(ctx context.Context) (*Client, error) {
	sharedOnce.Do(func() { sharedC, sharedErr = New(ctx) })
	return sharedC, sharedErr
}

func firstEnv(names ...string) string {
	for _, n := range names {
		if v := strings.TrimSpace(os.Getenv(n)); v != "" {
			return v
		}
	}
	return ""
}

// ParseURI splits an s3://bucket/key locator.
func ParseURI(uri string) (bucket, key string, err error) {
	if !strings.HasPrefix(uri, "s3://") {
		return "", "", fmt.Errorf("s3: not an s3:// locator: %q", uri)
	}
	bucket, key, _ = strings.Cut(strings.TrimPrefix(uri, "s3://"), "/")
	if bucket == "" {
		return "", "", fmt.Errorf("s3: no bucket in %q", uri)
	}
	return bucket, strings.Trim(key, "/"), nil
}

// Source is a ByteSource backed by ranged GETs.
type Source struct {
	api         *awss3.Client
	bucket, key string
	uri         string
	mu          sync.Mutex
	size        int64 // -1 until learned
}

var _ iosource.ByteSource = (*Source)(nil)

// Open returns a source for an s3:// locator.
func (c *Client) Open(ctx context.Context, uri string) (*Source, error) {
	bucket, key, err := ParseURI(uri)
	if err != nil {
		return nil, err
	}
	s := &Source{api: c.api, bucket: bucket, key: key, uri: uri, size: -1}
	// Learn the length now: it bounds the section reader that block-compressed
	// readers seek within, and it fails early and clearly when the object is
	// absent rather than midway through a query.
	if _, err := s.Size(); err != nil {
		return nil, err
	}
	return s, nil
}

// OpenReader streams a whole object. Used for index sidecars, which are read
// start to finish and are small enough not to warrant ranging.
func (c *Client) OpenReader(ctx context.Context, uri string) (io.ReadCloser, error) {
	bucket, key, err := ParseURI(uri)
	if err != nil {
		return nil, err
	}
	out, err := c.api.GetObject(ctx, &awss3.GetObjectInput{
		Bucket: aws.String(bucket), Key: aws.String(key),
	})
	if err != nil {
		return nil, fmt.Errorf("s3: get %s: %w", uri, err)
	}
	return out.Body, nil
}

// Size reports the object's length.
func (s *Source) Size() (int64, error) {
	s.mu.Lock()
	defer s.mu.Unlock()
	if s.size >= 0 {
		return s.size, nil
	}
	out, err := s.api.HeadObject(context.Background(), &awss3.HeadObjectInput{
		Bucket: aws.String(s.bucket), Key: aws.String(s.key),
	})
	if err != nil {
		return 0, fmt.Errorf("s3: head %s: %w", s.uri, err)
	}
	if out.ContentLength == nil {
		return 0, fmt.Errorf("s3: head %s: no content length", s.uri)
	}
	s.size = *out.ContentLength
	return s.size, nil
}

// ReadAt implements io.ReaderAt with a ranged GET. Safe for concurrent use, as
// ByteSource requires: each call is an independent request holding no state.
func (s *Source) ReadAt(p []byte, off int64) (int, error) {
	if len(p) == 0 {
		return 0, nil
	}
	if off < 0 {
		return 0, fmt.Errorf("s3: read %s: negative offset %d", s.uri, off)
	}
	out, err := s.api.GetObject(context.Background(), &awss3.GetObjectInput{
		Bucket: aws.String(s.bucket),
		Key:    aws.String(s.key),
		Range:  aws.String(fmt.Sprintf("bytes=%d-%d", off, off+int64(len(p))-1)),
	})
	if err != nil {
		// Reading at or past the end is EOF, not a failure — block-compressed
		// readers probe past the end legitimately.
		if isRangeUnsatisfiable(err) {
			return 0, io.EOF
		}
		return 0, fmt.Errorf("s3: read %s at %d: %w", s.uri, off, err)
	}
	defer out.Body.Close()

	n, err := io.ReadFull(out.Body, p)
	if errors.Is(err, io.ErrUnexpectedEOF) || errors.Is(err, io.EOF) {
		// io.ReaderAt requires a non-nil error on a short read.
		return n, io.EOF
	}
	return n, err
}

// Close releases the source. Each read is its own request, so there is nothing
// held open; this exists to satisfy ByteSource.
func (s *Source) Close() error { return nil }

func isRangeUnsatisfiable(err error) bool {
	var re interface{ HTTPStatusCode() int }
	if errors.As(err, &re) && re.HTTPStatusCode() == 416 {
		return true
	}
	return strings.Contains(err.Error(), "InvalidRange")
}

// PutForTest uploads an object. It exists so packages in this repository can
// stage their own fixtures in integration tests without taking on a write-side
// S3 dependency of their own; provisioning proper — multipart, checksums,
// abort-on-failure — belongs to whatever tool writes the data.
func PutForTest(ctx context.Context, bucket, key string, body io.Reader) error {
	c, err := Shared(ctx)
	if err != nil {
		return err
	}
	_, err = c.api.PutObject(ctx, &awss3.PutObjectInput{
		Bucket: aws.String(bucket), Key: aws.String(key), Body: body,
	})
	return err
}
