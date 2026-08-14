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

// VolumeManifestName is the store table that marks a conversion complete.
const VolumeManifestName = "manifest"

// VolumeManifestFile is the manifest's name inside the store directory.
const VolumeManifestFile = VolumeManifestName + ".json.gz"

// VolumeManifestVersion is the schema version of the manifest document itself.
// It exists so a future reader can refuse a layout it does not understand
// rather than silently misreading one; nothing about the parquet tables is
// versioned by it.
const VolumeManifestVersion = 1

// VolumeManifestPath returns the manifest file for a store.
func VolumeManifestPath(base string) string { return joinStore(base, VolumeManifestFile) }

// VolumeManifest records what a conversion actually wrote.
//
// It exists because the parquet tables cannot answer one question about
// themselves: whether the run that produced them reached the end. A footer is
// written only by the writer's Close, so a table that parses is a table that
// was finished -- but a *set* of finished tables says nothing about how much of
// the intended input went into them. The metadata already in the calls file is
// worse than silent on this: MetaSource and MetaContigs are stamped at
// construction, before a single record is read, so a store holding three
// chromosomes of a twenty-two-input conversion names all twenty-two inputs and
// declares all twenty-two contigs.
//
// So the manifest is written last, after every table is closed, and only then.
// Its presence is the claim "this conversion ran to completion"; the counts in
// it are what a reader checks that claim against.
//
// Chromosomes is the part that earns the file. It turns "this store is smaller
// than I expected" into "chr4 through chr22 have no sites", which is the only
// thing that can contradict the intent recorded at construction.
type VolumeManifest struct {
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

	Params  ManifestParams       `json:"params"`
	Samples []string             `json:"samples"`
	Counts  ManifestCounts       `json:"counts"`
	Tables  map[string]TableInfo `json:"tables"`

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

	// DepthBands are the boundaries at which callable runs were broken, empty
	// when they were not banded.
	//
	// Recorded because two stores banded differently do not mean the same thing
	// by a run: the same MinDP describes a tight class in one and a whole
	// chromosome arm in the other, and a consumer comparing them across parts
	// has no other way to find out.
	DepthBands []int32 `json:"depth_bands,omitempty"`

	// Format are the FORMAT fields captured onto the ALT calls, if any. The
	// same reasoning as Info: absence and a zero read alike from a typed
	// reader, and only this can tell them apart.
	Format []FormatField `json:"format,omitempty"`

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

// UnmarshalJSON reads a volume manifest, accepting the pre-rename key for the
// table index.
//
// The three parquet tables used to be "members", and the word had to give: it
// named one of the tables inside a volume AND one of the volumes inside a
// store, which are different axes entirely. Volumes kept it; tables got the
// name they always were.
//
// Stores written before the rename carry the table index under "members", and
// reading it is not cosmetic. The row counts there are what catch a table that
// parses but belongs to a different conversion. Ignoring the old key would not
// fail loudly -- it would leave an empty map, and the check would pass over
// nothing at all.
//
// Written only as "tables". This reads both.
func (m *VolumeManifest) UnmarshalJSON(b []byte) error {
	type raw VolumeManifest
	var v struct {
		raw
		LegacyTables map[string]TableInfo `json:"members"`
	}
	if err := json.Unmarshal(b, &v); err != nil {
		return err
	}
	*m = VolumeManifest(v.raw)
	if len(m.Tables) == 0 {
		m.Tables = v.LegacyTables
	}
	return nil
}

// TableInfo is one table's size on disk, as written.
//
// Rows is checked against the table's own parquet footer at open. That is the
// check which catches a table that is present and well-formed but does not
// belong to this store -- one copied in from another conversion, say -- which
// nothing else can see, since sites and regions carry no metadata of their own.
type TableInfo struct {
	Rows  int64 `json:"rows"`
	Bytes int64 `json:"bytes"`

	// Shards, when this table is split by coordinate. Empty means one file
	// under the table's own name, which is every store written before
	// splitting existed and every store small enough not to want it.
	//
	// THE INDEX IS THE POINT. Parquet already prunes row groups by their
	// statistics, but the footer carrying those statistics is at the end of the
	// file -- so a locus query against a whole-genome table must fetch and
	// parse a footer describing hundreds of gigabytes before it can decide to
	// read none of it. A shard index in the manifest answers the same question
	// having read only the manifest, and over object storage that is the
	// difference between one small GET and a large one.
	Shards []ShardInfo `json:"shards,omitempty"`
}

// ShardInfo is one coordinate range of a split table.
//
// Shards are ALIGNED ACROSS MEMBERS: shard k of calls, sites and regions cover
// the same interval, so a locus query reads the k'th of each and never has to
// reconcile two different partitionings. A shard never spans chromosomes, which
// is what lets First and Last be compared without a chromosome check.
type ShardInfo struct {
	// Name is the file, relative to the store: "calls/00007.parquet".
	Name  string `json:"name"`
	Chrom string `json:"chrom"`

	// First and Last bound the SITE coordinates this shard covers, inclusive.
	//
	// For regions the bound is on the interval, not on any single run: a run is
	// broken at shard boundaries when it is written, so every run lies wholly
	// inside one shard. That is what makes a shard answerable on its own --
	// without it, a query would have to read every earlier shard in case a run
	// started there and reached in.
	First int32 `json:"first"`
	Last  int32 `json:"last"`

	Rows  int64 `json:"rows"`
	Bytes int64 `json:"bytes"`
}

// Covers reports whether a position falls in this shard.
func (s ShardInfo) Covers(chrom string, pos int32) bool {
	return SameChrom(s.Chrom, chrom) && pos >= s.First && pos <= s.Last
}

// Overlaps reports whether a coordinate span touches this shard.
func (s ShardInfo) Overlaps(chrom string, first, last int32) bool {
	return SameChrom(s.Chrom, chrom) && last >= s.First && first <= s.Last
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

// WriteVolumeManifest writes the manifest for the store at base.
//
// The write is atomic: the document goes to a temporary file in the store
// directory and is renamed into place. A manifest is a completion marker, so a
// half-written one would be worse than none at all -- it would claim the store
// is finished while being unreadable proof of the opposite.
func WriteVolumeManifest(base string, m VolumeManifest) error {
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
// A sink does not have that problem to solve. A local table is created and
// closed in one call here, and an object store's table does not exist at all
// until its upload completes -- so "appears whole or not at all" is the sink's
// property rather than something this function has to arrange.
func WriteManifestTo(sink Sink, m VolumeManifest) error {
	f, err := sink.Create(VolumeManifestFile)
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

func ReadVolumeManifest(base string) (*VolumeManifest, error) {
	return ReadVolumeManifestContext(context.Background(), base)
}

// ReadVolumeManifestContext reads the manifest from any locator.
func ReadVolumeManifestContext(ctx context.Context, base string) (*VolumeManifest, error) {
	locator := VolumeManifestPath(TrimStoreSuffix(base))
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
	var m VolumeManifest
	if err := json.NewDecoder(zr).Decode(&m); err != nil {
		return nil, fmt.Errorf("reading %s: %w", locator, err)
	}
	return &m, nil
}
