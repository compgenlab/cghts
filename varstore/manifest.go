package varstore

import (
	"compress/gzip"
	"context"
	"encoding/json"
	"fmt"
	"io"
	"time"

	"github.com/compgenlab/cghts/iosource"
)

// ManifestMember is the store member that marks a conversion complete.
const ManifestMember = "manifest"

// ManifestFile is the manifest's name inside the store directory.
const ManifestFile = ManifestMember + ".json.gz"

// ManifestVersion is the schema version of the manifest document itself.
// It exists so a future reader can refuse a layout it does not understand
// rather than silently misreading one; nothing about the parquet members is
// versioned by it.
const ManifestVersion = 1

// ManifestPath returns the manifest file for a store.
func ManifestPath(base string) string { return joinStore(base, ManifestFile) }

// Manifest records what a conversion actually wrote.
//
// It exists because the parquet members cannot answer one question about
// themselves: whether the run that produced them reached the end. A footer is
// written only by the writer's Close, so a member that parses is a member that
// was finished -- but a *set* of finished members says nothing about how much of
// the intended input went into them. The metadata already in the calls file is
// worse than silent on this: MetaSource and MetaContigs are stamped at
// construction, before a single record is read, so a store holding three
// chromosomes of a twenty-two-input conversion names all twenty-two inputs and
// declares all twenty-two contigs.
//
// So the manifest is written last, after every member is closed, and only then.
// Its presence is the claim "this conversion ran to completion"; the counts in
// it are what a reader checks that claim against.
//
// Chromosomes is the part that earns the file. It turns "this store is smaller
// than I expected" into "chr4 through chr22 have no sites", which is the only
// thing that can contradict the intent recorded at construction.
type Manifest struct {
	FormatVersion int       `json:"format_version"`
	Complete      bool      `json:"complete"`
	Created       time.Time `json:"created"`

	Program string   `json:"program,omitempty"`
	Command string   `json:"command,omitempty"`
	Sources []string `json:"sources,omitempty"`

	// Meta is what the caller said the store *is*, where the fields above record
	// how it was made. See ReservedMetaKeys. Omitted entirely when empty, so a
	// store converted without any reads back identical to one written before
	// this field existed -- absent means "not stated", never "stated as nothing".
	Meta map[string]string `json:"meta,omitempty"`

	Params  ManifestParams        `json:"params"`
	Samples []string              `json:"samples"`
	Counts  ManifestCounts        `json:"counts"`
	Members map[string]MemberInfo `json:"members"`

	Chromosomes     []ChromCensus `json:"chromosomes"`
	ContigsDeclared []string      `json:"contigs_declared,omitempty"`
}

// ManifestParams records the conversion settings that change what the store
// means, as opposed to how it is encoded.
type ManifestParams struct {
	MinDP         int32         `json:"min_dp"`
	NoCallable    bool          `json:"no_callable"`
	RowGroupSize  int64         `json:"row_group_size"`
	SpanSemantics SpanSemantics `json:"span_semantics"`

	// Info are the INFO fields captured into sites.parquet, if any.
	//
	// THIS IS THE ONLY PLACE ABSENCE IS ANSWERABLE. A typed reader gets a zero
	// both for a column that is not in the file and for one holding zero, and
	// this store's own schema records the same trap for RefEnd -- "zero in
	// stores written before this column existed; treat that as unknown". A
	// consumer asking "does this store know its imputation quality" must be
	// able to get "no" rather than "0.0".
	Info []InfoField `json:"info,omitempty"`
}

// ManifestCounts are store-wide totals.
type ManifestCounts struct {
	Samples int   `json:"samples"`
	Calls   int64 `json:"calls"`
	Sites   int64 `json:"sites"`
	Regions int64 `json:"regions"`
}

// MemberInfo is one member's size on disk, as written.
//
// Rows is checked against the member's own parquet footer at open. That is the
// check which catches a member that is present and well-formed but does not
// belong to this store -- one copied in from another conversion, say -- which
// nothing else can see, since sites and regions carry no metadata of their own.
type MemberInfo struct {
	Rows  int64 `json:"rows"`
	Bytes int64 `json:"bytes"`
}

// ChromCensus is what one chromosome contributed, as counted by the writer
// while it wrote rather than as reported by the caller.
type ChromCensus struct {
	Name     string `json:"name"`
	Sites    int64  `json:"sites"`
	Calls    int64  `json:"calls"`
	FirstPos int32  `json:"first_pos"`
	LastPos  int32  `json:"last_pos"`
}

// WriteManifest writes the manifest for the store at base.
//
// The write is atomic: the document goes to a temporary file in the store
// directory and is renamed into place. A manifest is a completion marker, so a
// half-written one would be worse than none at all -- it would claim the store
// is finished while being unreadable proof of the opposite.
func WriteManifest(base string, m Manifest) error {
	sink, err := OpenSink(base)
	if err != nil {
		return err
	}
	return WriteManifestTo(sink, m)
}

// WriteManifestTo writes the manifest through a sink.
//
// SIMPLER THAN THE LOCAL PATH IT REPLACED, and deliberately so. Writing a
// manifest used to mean a temporary file, an fsync and a rename, because a
// half-written local file is visible while it is being written and a manifest
// is a completion marker: one that appeared before its own contents would claim
// a store is finished while being proof of the opposite.
//
// A sink does not have that problem to solve. A local member is created and
// closed in one call here, and an object store's member does not exist at all
// until its upload completes -- so "appears whole or not at all" is the sink's
// property rather than something this function has to arrange.
func WriteManifestTo(sink Sink, m Manifest) error {
	f, err := sink.Create(ManifestFile)
	if err != nil {
		return fmt.Errorf("creating the manifest in %s: %w", sink.Describe(), err)
	}
	zw := gzip.NewWriter(f)
	enc := json.NewEncoder(zw)
	enc.SetIndent("", "  ")
	if err := enc.Encode(m); err != nil {
		f.Close()
		return fmt.Errorf("encoding manifest: %w", err)
	}
	if err := zw.Close(); err != nil {
		f.Close()
		return fmt.Errorf("compressing manifest: %w", err)
	}
	if err := f.Close(); err != nil {
		return fmt.Errorf("writing the manifest to %s: %w", sink.Describe(), err)
	}
	return nil
}

func ReadManifest(base string) (*Manifest, error) {
	return ReadManifestContext(context.Background(), base)
}

// ReadManifestContext reads the manifest from any locator.
func ReadManifestContext(ctx context.Context, base string) (*Manifest, error) {
	locator := ManifestPath(TrimStoreSuffix(base))
	src, err := iosource.Open(ctx, locator)
	if err != nil {
		return nil, err
	}
	defer src.Close()
	size, err := src.Size()
	if err != nil {
		return nil, err
	}
	zr, err := gzip.NewReader(io.NewSectionReader(src, 0, size))
	if err != nil {
		return nil, fmt.Errorf("reading %s: %w", locator, err)
	}
	defer zr.Close()
	var m Manifest
	if err := json.NewDecoder(zr).Decode(&m); err != nil {
		return nil, fmt.Errorf("reading %s: %w", locator, err)
	}
	return &m, nil
}
