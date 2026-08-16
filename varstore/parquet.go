package varstore

import (
	"context"
	"errors"
	"fmt"
	"io"
	"io/fs"
	"iter"
	"os"
	"sort"
	"strconv"
	"strings"
	"time"

	"github.com/parquet-go/parquet-go"
	"github.com/parquet-go/parquet-go/compress"
	"github.com/parquet-go/parquet-go/compress/snappy"
	"github.com/parquet-go/parquet-go/compress/uncompressed"
	"github.com/parquet-go/parquet-go/compress/zstd"
)

// Parquet file metadata keys. The sample roster lives here rather than in a
// fourth file: Classify needs every sample, but the calls file only ever
// mentions carriers, so a sample carrying nothing anywhere would otherwise be
// invisible and get reported as not-assayed everywhere.
const (
	MetaSamples = "cgkit.samples"        // newline-separated sample ids, source order
	MetaMinDP   = "cgkit.min_dp"         // callable threshold used at conversion
	MetaProgram = "cgkit.program"        // cgkit version
	MetaCommand = "cgkit.command"        // full command line
	MetaSource  = "cgkit.source"         // input filename
	MetaNoCall  = "cgkit.nocallable"     // "1" when regions are absent by request
	MetaSpans   = "cgkit.span_semantics" // see SpanSemantics

	// MetaContigs holds the source's ##contig header lines verbatim, newline
	// separated. Kept because a store is expected to be exported back to VCF, and
	// those lines are how a VCF says which reference it was called against -- without
	// them the export is not self-describing, and cannot carry contig lengths at all.
	// Derived from the calls alone, the best available is a bare ID for whichever
	// contigs a query happened to name.
	MetaContigs = "cgkit.contigs"
)

// MetaPrefix namespaces caller-supplied metadata inside the parquet key/value
// metadata, so a WriterOpts.Meta key can never shadow one of the keys above --
// "cgkit.meta.source" and "cgkit.source" are different keys, and only the
// latter is provenance the writer vouches for.
const MetaPrefix = "cgkit.meta."

// Reserved metadata keys. Nothing rejects a key outside this set; they exist so
// that the writer, a reader, and a CLI exposing them as flags agree on spelling
// rather than each inventing "ref" or "genome" for the same fact.
//
// What earns a place here is a fact the store cannot recover about itself.
// Sources, sample roster, creation time and the per-chromosome census are all
// recorded already; the assembly is not, and neither is the name of the release
// a set of per-chromosome inputs came from.
const (
	// MetaKeyDataset names the release or callset the store was built from --
	// the thing a directory of 24 per-chromosome VCFs is collectively called.
	MetaKeyDataset = "dataset"

	// MetaKeyReference names the assembly the calls were made against. The
	// ##contig lines in MetaContigs carry lengths that imply it, but implying is
	// not declaring, and a store read against the wrong assembly does not fail --
	// it answers with coordinates that mean something else.
	MetaKeyReference = "reference"

	// MetaKeyCaller names the variant caller and version. It bears on queries
	// rather than just on provenance: DP, GQ and RGQ mean different things
	// between callers, and those are the fields Gate acts on.
	MetaKeyCaller = "caller"

	// MetaKeyAccession is the study or dataset accession (phs000000, EGAD...).
	MetaKeyAccession = "accession"

	// MetaKeyURL is where the data was retrieved from. Kept separate from
	// MetaKeyAccession because an accession identifies and a URL locates: a
	// store can have either without the other, and folding them loses whichever
	// was not written.
	MetaKeyURL = "url"

	// MetaKeyVersion is the dataset's own release version, as opposed to the
	// cgkit version that did the conversion -- which is already in MetaProgram.
	MetaKeyVersion = "version"

	// MetaKeyDescription is free text.
	MetaKeyDescription = "description"
)

// ReservedMetaKeys lists the reserved keys in the order a report should present
// them, most identifying first. A CLI generating one flag per key should range
// over this rather than repeat the list, so the two cannot drift.
var ReservedMetaKeys = []string{
	MetaKeyDataset,
	MetaKeyReference,
	MetaKeyCaller,
	MetaKeyAccession,
	MetaKeyURL,
	MetaKeyVersion,
	MetaKeyDescription,
}

// ValidMetaKey reports whether k is usable as a metadata key.
//
// Keys are constrained where values are not. A value is a claim by the caller
// and is recorded verbatim -- the writer cannot know whether "GRCh38" is true,
// and normalizing it would turn a caller's assertion into the library's. A key
// is an identifier that has to survive a round trip through a parquet metadata
// key and a JSON object table, and be greppable afterwards, so it is held to
// lowercase [a-z0-9_-]. The dot is excluded specifically: MetaPrefix already
// uses it as a separator, and a key containing one would read back ambiguously.
func ValidMetaKey(k string) bool {
	if k == "" {
		return false
	}
	for _, r := range k {
		switch {
		case r >= 'a' && r <= 'z', r >= '0' && r <= '9', r == '_', r == '-':
		default:
			return false
		}
	}
	return true
}

// copyMeta returns an independent copy of m, or nil when m is empty so the
// manifest's omitempty elides the field rather than writing "meta": {}. The copy
// matters because the manifest outlives the WriterOpts the caller passed in.
func copyMeta(m map[string]string) map[string]string {
	if len(m) == 0 {
		return nil
	}
	out := make(map[string]string, len(m))
	for k, v := range m {
		out[k] = v
	}
	return out
}

// sortedKeys returns m's keys in sorted order, so that metadata is stamped
// deterministically and two conversions given the same inputs produce the same
// bytes. Go map iteration would not.
func sortedKeys(m map[string]string) []string {
	if len(m) == 0 {
		return nil
	}
	ks := make([]string, 0, len(m))
	for k := range m {
		ks = append(ks, k)
	}
	sort.Strings(ks)
	return ks
}

// validateMeta checks every key, naming all the bad ones rather than only the
// first: a caller passing several flags wants one round trip, not one per typo.
func validateMeta(m map[string]string) error {
	var bad []string
	for k := range m {
		if !ValidMetaKey(k) {
			bad = append(bad, strconv.Quote(k))
		}
	}
	if len(bad) == 0 {
		return nil
	}
	sort.Strings(bad)
	return fmt.Errorf("invalid metadata key(s) %s: keys must be non-empty and use only "+
		"lowercase letters, digits, underscore and hyphen", strings.Join(bad, ", "))
}

// SpanSemantics records what the intervals in the regions file are entitled to
// claim, which depends entirely on what the source format asserted.
type SpanSemantics string

const (
	// SpansSites means the intervals only mark catalog sites at which a sample
	// was called. Nothing is claimed about the bases in between. This is all a
	// plain VCF supports, and it confines queries to the sites catalog.
	SpansSites SpanSemantics = "sites"

	// SpansBlocks means the intervals came from gVCF reference blocks, which
	// are positive statements about whole spans. Only such a store may answer
	// for positions absent from the catalog. Not yet produced by any converter.
	SpansBlocks SpanSemantics = "blocks"
)

// CodecFor maps a --compression flag value to a parquet codec.
func CodecFor(name string) (compress.Codec, error) {
	switch strings.ToLower(name) {
	case "zstd", "":
		return &zstd.Codec{}, nil
	case "snappy":
		return &snappy.Codec{}, nil
	case "none", "uncompressed":
		return &uncompressed.Codec{}, nil
	}
	return nil, fmt.Errorf("unknown compression %q (use zstd, snappy, or none)", name)
}

// WriterOpts configures a Parquet store writer.
type WriterOpts struct {
	// Sink is where the tables are written. Nil means "work it out from the
	// base", which is what every caller wants; it is settable so a test can
	// supply one without a filesystem or a bucket.
	Sink Sink

	Codec        compress.Codec
	RowGroupSize int64
	Samples      []string
	MinDP        int32
	NoCallable   bool
	Program      string

	// Coverage opens the optional coverage table, for genomic block spans that
	// answer OFF the sites catalog. Off by default: a store carrying an empty
	// one asserts "covered nowhere", where an absent one asserts nothing.
	Coverage bool

	// MaxGap is the largest uncovered stretch a coverage block was permitted to
	// span when it was built, recorded so a reader knows what "covered" meant
	// here. Zero means the blocks are the source's own, unmerged.
	//
	// Recorded for the same reason MinDP is: two stores built with different
	// tolerances do not mean the same thing by "covered", and a roll-up
	// spanning both has to be able to see that they differ.
	MaxGap  int32
	Command string

	// ShardSites splits each table every N sites, so a locus query reads one
	// small file rather than pruning row groups inside a large one. Zero writes
	// the store unsplit, which is what every store before this is.
	//
	// The win is not the pruning -- Parquet already prunes row groups -- it is
	// that the statistics doing the pruning live in the FOOTER, so deciding not
	// to read a whole-genome table still costs fetching and parsing a footer
	// describing it. A shard index in the manifest decides the same thing having
	// read only the manifest.
	ShardSites int64

	// BeforeRotate is called just before a shard is closed, so a caller holding
	// state that spans sites can flush it into the shard it belongs to.
	//
	// The converter's open callable runs are exactly that. A run must lie wholly
	// inside one shard -- otherwise a query would have to read every earlier
	// shard in case a run started there and reached in -- so runs are broken
	// here, exactly as they are already broken at depth-band boundaries.
	BeforeRotate func() error

	// DepthBands are the boundaries at which a callable run is broken, so each
	// run spans one depth class and its MinDP is a tight bound rather than the
	// worst moment in an arbitrary stretch. Empty leaves runs unbanded, which is
	// what every store written before this existed is.
	DepthBands []int32

	// Format are FORMAT fields to capture onto the ALT calls, each as its own
	// typed column in calls.parquet.
	//
	// Empty by default and deliberately so: calls is the large table, so a
	// column here costs roughly a hundred times what the same column costs on
	// the sites catalog. See format.go.
	Format []FormatField

	// Info are INFO fields to capture from the source into sites.parquet, each
	// as its own typed column. Empty -- the usual case -- writes exactly the
	// schema this store has always had. See info.go for why these are columns
	// and why the source key is preserved rather than normalized.
	Info []InfoField

	// Sources are the inputs the conversion was asked to consume, in order.
	// The manifest keeps the list rather than a joined string, because "which
	// of these actually went in" is the question a partial store raises.
	Sources []string

	// Contigs are the source's ##contig lines, verbatim, in header order. Callers
	// converting several inputs should pass the union: a per-chromosome callset
	// declares only its own contig in each file, so taking one file's would lose
	// the rest.
	Contigs []string

	// Spans declares what the run intervals may claim. Defaults to SpansSites,
	// which is all a plain VCF can support.
	Spans SpanSemantics

	// Meta is caller-supplied metadata describing what the store *is*, as
	// opposed to how it was built: see ReservedMetaKeys for the keys with an
	// agreed meaning. The map is open, so a caller may record anything; keys are
	// validated by NewWriter, values never are.
	Meta map[string]string
}

// Writer builds the three files of a Parquet store. Rows are buffered and
// flushed in batches so memory stays bounded no matter how large the input is.
type Writer struct {
	calls    *parquet.GenericWriter[Call]
	sites    *parquet.GenericWriter[Site]
	regions  *parquet.GenericWriter[CalledSiteRun]
	coverage *parquet.GenericWriter[CoverageBlock]

	// The dynamic sites writer, used only when INFO fields are captured. Exactly
	// one of sites and sitesAny is ever non-nil.
	sitesAny   *parquet.GenericWriter[any]
	siteSchema *parquet.Schema

	// The same, for calls, when FORMAT fields are captured.
	callsAny   *parquet.GenericWriter[any]
	callSchema *parquet.Schema
	callRows   []any
	fmtScratch []infoSlots

	callBuf     []Call
	siteBuf     []Site
	regionBuf   []CalledSiteRun
	coverageBuf []CoverageBlock

	// Rows for the dynamic path, allocated once and refilled. A whole-genome
	// catalog is tens of millions of sites and a map per row would be the
	// dominant cost of writing one.
	siteRows []any

	// Sharding: the coordinate ranges each table is split into.
	//
	// A shard is cut every ShardSites sites and at every chromosome change,
	// which is what keeps First and Last comparable without a chromosome test.
	// The three tables are cut TOGETHER on the same site boundaries, so shard
	// k of calls, sites and regions cover the same interval and a locus query
	// reads the k'th of each rather than reconciling two partitionings.
	shardIdx   int
	shardChrom string
	shardFirst int32
	shardLast  int32
	shardSites int64
	shardRows  map[string]int64
	shards     map[string][]ShardInfo

	// Backing storage for captured values, one slot per buffered row per field.
	//
	// Optional columns are written through POINTERS, because parquet-go encodes
	// an optional field's zero value as null -- so a plain float64(0) would
	// store R2=0 as "no R2 here", which is a different claim and the one this
	// column exists to avoid making. Pointers into a preallocated batch keep
	// that correct without an allocation per value.
	infoScratch []infoSlots

	// The tables being written, in creation order, with the names they will
	// have. Held so Close can finish them and abort can undo them -- and named
	// rather than handled, because a remote table has no file to name itself.
	sink   Sink
	tables []tableSink

	// base and opts are kept so Finish can describe the store without the
	// caller having to hand back what it already passed to NewWriter.
	base string
	opts WriterOpts

	// The per-chromosome census is accumulated here rather than by the caller,
	// so it records what was actually written. A count supplied by the
	// converter would be a second opinion that could drift from the rows; this
	// one cannot.
	chroms   []ChromCensus
	chromIdx map[string]int

	NCalls    int64
	NSites    int64
	NRegions  int64
	NCoverage int64
}

// census returns the running tally for chrom, creating it on first sight.
// Chromosome order follows the order they were written, which is the store's
// own order and not lexicographic.
// Keyed on the canonical chromosome, not the spelling: a conversion whose
// inputs mix "chr1" and "1" was producing two entries for one chromosome, in
// the field the manifest doc calls the part that earns the file. Name keeps the
// first spelling seen, because the store records the source's own naming rather
// than rewriting it.
func (w *Writer) census(chrom string, pos int32) *ChromCensus {
	key := CanonKey(chrom)
	if i, ok := w.chromIdx[key]; ok {
		c := &w.chroms[i]
		if pos < c.FirstPos {
			c.FirstPos = pos
		}
		if pos > c.LastPos {
			c.LastPos = pos
		}
		return c
	}
	if w.chromIdx == nil {
		w.chromIdx = map[string]int{}
	}
	w.chromIdx[key] = len(w.chroms)
	w.chroms = append(w.chroms, ChromCensus{Name: chrom, FirstPos: pos, LastPos: pos})
	return &w.chroms[len(w.chroms)-1]
}

const batchSize = 8192

// tableSink is one store table being written: the writer it is fed through
// and the file name it will have. Named around the sink because `table` and
// `openTable` both already mean something on the reading side.
type tableSink struct {
	name string
	sw   io.WriteCloser

	// closed once, by whoever finished the shard it belongs to. Without this a
	// split store closed each shard's sink at the rotation AND again at the end,
	// and the second close is an error the writer reported as a failed
	// conversion.
	closed bool
}

// NewWriter creates the store directory at base and the table files inside it.
//
// The directory is created here rather than left to the caller because a store
// *is* its directory: there is no longer a spelling of base that means anything
// else, so requiring a separate mkdir would only be a way to get it wrong.
// EnsureStoreDir remains for callers that want the directory to exist earlier,
// and is idempotent.
func NewWriter(base string, o WriterOpts) (*Writer, error) {
	// Before anything on disk is touched. The next step removes the previous
	// run's manifest, so a bad key caught later would already have unmade a
	// readable store in order to reject the run that was meant to replace it.
	if err := validateMeta(o.Meta); err != nil {
		return nil, err
	}
	// Same reason as validateMeta: refuse before the previous manifest is
	// cleared, so a bad --info spelling cannot unmake a readable store on its
	// way to rejecting the run meant to replace it.
	if err := ValidateInfo(o.Info); err != nil {
		return nil, err
	}
	o.Info = normalizeInfo(o.Info)
	if err := ValidateFormat(o.Format); err != nil {
		return nil, err
	}
	o.Format = normalizeFormat(o.Format)
	if o.Codec == nil {
		o.Codec = &zstd.Codec{}
	}
	if o.RowGroupSize <= 0 {
		o.RowGroupSize = 250_000
	}
	if o.Spans == "" {
		o.Spans = SpansSites
	}
	sink := o.Sink
	if sink == nil {
		var err error
		if sink, err = OpenSink(base); err != nil {
			return nil, err
		}
	}
	// Remove any manifest before touching the tables. From here until Finish
	// the store is under construction and must not carry a completion marker --
	// otherwise a --force re-run that dies partway would leave the previous
	// run's manifest vouching for this run's half-written tables.
	if err := sink.Remove(VolumeManifestFile); err != nil {
		return nil, fmt.Errorf("clearing previous manifest: %w", err)
	}
	w := &Writer{base: base, opts: o, sink: sink}

	w.shards = map[string][]ShardInfo{}
	w.shardRows = map[string]int64{}

	if err := w.openShard(); err != nil {
		return nil, err
	}
	return w, nil
}

// tableName is where a table's rows go: its own file when unsplit, or the
// current shard's file when split.
func (w *Writer) tableName(table string) string {
	if w.opts.ShardSites <= 0 {
		return TableFile(table)
	}
	return ShardFile(table, w.shardIdx)
}

// openShard creates the three table writers for the current shard.
//
// Called once for an unsplit store and once per shard for a split one, so the
// two paths cannot drift: a split store's first shard is written by exactly the
// code that writes an unsplit store's only file.
func (w *Writer) openShard() error {
	o := w.opts
	open := func(name string) (io.Writer, error) {
		f, err := w.sink.Create(name)
		if err != nil {
			_ = w.abort()
			return nil, err
		}
		w.tables = append(w.tables, tableSink{name: name, sw: f})
		return f, nil
	}

	opts := []parquet.WriterOption{
		parquet.Compression(o.Codec),
		parquet.MaxRowsPerRowGroup(o.RowGroupSize),
	}
	// Declare the order the rows are already in. Input is coordinate sorted and
	// written in stream order, so saying so lets a reader trust the per-group
	// min/max on pos, which is what makes locus lookups skip whole groups.
	sortedByLocus := append(append([]parquet.WriterOption{}, opts...),
		parquet.SortingWriterConfig(parquet.SortingColumns(
			parquet.Ascending("chrom"), parquet.Ascending("pos"))))
	// sample_id is high-cardinality and unsorted, so statistics cannot bound it;
	// a bloom filter answers "is this sample absent from this group" exactly,
	// which is what a --sample query needs.
	callOpts := append(append([]parquet.WriterOption{}, sortedByLocus...),
		parquet.BloomFilters(parquet.SplitBlockFilter(10, "sample_id")))

	cf, err := open(w.tableName(CallsTable))
	if err != nil {
		return err
	}
	if len(o.Format) > 0 {
		w.callSchema = callSchemaWith(o.Format)
		w.fmtScratch = newInfoSlots(formatAsInfo(o.Format), batchSize)
		w.callsAny = parquet.NewGenericWriter[any](cf,
			append(append([]parquet.WriterOption{}, callOpts...), w.callSchema)...)
		w.calls = nil
	} else {
		w.calls = parquet.NewGenericWriter[Call](cf, callOpts...)
		w.callsAny = nil
	}
	// THE STORE'S METADATA GOES ON EVERY SHARD, not only the first. A shard is a
	// parquet file somebody may open on its own -- with DuckDB, with pyarrow,
	// with anything -- and one that could not say which samples it holds or what
	// depth gate it was written under would be a file with no provenance.
	w.setCallMeta(MetaSamples, strings.Join(o.Samples, "\n"))
	w.setCallMeta(MetaMinDP, fmt.Sprint(o.MinDP))
	w.setCallMeta(MetaProgram, o.Program)
	w.setCallMeta(MetaCommand, o.Command)
	w.setCallMeta(MetaSource, strings.Join(o.Sources, ", "))
	if len(o.Contigs) > 0 {
		w.setCallMeta(MetaContigs, strings.Join(o.Contigs, "\n"))
	}
	w.setCallMeta(MetaSpans, string(o.Spans))
	for _, k := range sortedKeys(o.Meta) {
		w.setCallMeta(MetaPrefix+k, o.Meta[k])
	}
	if o.NoCallable {
		w.setCallMeta(MetaNoCall, "1")
	}

	sf, err := open(w.tableName(SitesTable))
	if err != nil {
		return err
	}
	if len(o.Info) > 0 {
		w.siteSchema = siteSchemaWith(o.Info)
		w.infoScratch = newInfoSlots(o.Info, batchSize)
		w.sitesAny = parquet.NewGenericWriter[any](sf,
			append(append([]parquet.WriterOption{}, sortedByLocus...), w.siteSchema)...)
		w.sites = nil
	} else {
		w.sites = parquet.NewGenericWriter[Site](sf, sortedByLocus...)
		w.sitesAny = nil
	}

	rf, err := open(w.tableName(RegionsTable))
	if err != nil {
		return err
	}
	w.regions = parquet.NewGenericWriter[CalledSiteRun](rf, append(append([]parquet.WriterOption{}, opts...),
		parquet.SortingWriterConfig(parquet.SortingColumns(
			parquet.Ascending("chrom"), parquet.Ascending("start"))),
		parquet.BloomFilters(parquet.SplitBlockFilter(10, "sample_id")))...)

	// OPENED ONLY WHEN ASKED FOR, so an absent coverage table means "nobody
	// said" rather than "covered nowhere". Creating an empty one by default
	// would assert the second, which is the claim that turns an unknown into a
	// wrong answer -- every off-catalog position would read as never covered
	// for every sample, which is exactly what the table exists to fix.
	if w.opts.Coverage {
		cf, err := open(w.tableName(CoverageTable))
		if err != nil {
			return err
		}
		w.coverage = parquet.NewGenericWriter[CoverageBlock](cf, append(append([]parquet.WriterOption{}, opts...),
			parquet.SortingWriterConfig(parquet.SortingColumns(
				parquet.Ascending("chrom"), parquet.Ascending("start"))),
			parquet.BloomFilters(parquet.SplitBlockFilter(10, "sample_id")))...)
	}

	w.shardSites = 0
	w.shardChrom = ""
	w.shardRows = map[string]int64{}
	return nil
}

// abort closes and removes any files opened so far, used when construction
// fails partway; a partial store is worse than none, since the set is meant to
// be inseparable.
//
// The removal error is reported rather than swallowed. If the unlink fails --
// a read-only mount, a sticky-bit directory, a stale NFS handle -- what is left
// behind is a partial store that the overwrite guard will then treat as a real
// one and refuse to replace. That is worth saying out loud rather than
// returning as though nothing happened.
func (w *Writer) abort() error {
	var errs []error
	// The manifest is written last and cleared first, so it should not be here.
	// A --force re-run over a previous store is the case where a stale one could
	// survive, and a completion marker sitting beside discarded tables is the
	// worst outcome available.
	if err := w.sink.Remove(VolumeManifestFile); err != nil {
		errs = append(errs, fmt.Errorf("removing %s: %w", VolumeManifestFile, err))
	}
	aborter, remote := w.sink.(Aborter)
	for _, m := range w.tables {
		// ABANDON where a table only becomes visible when it is finished, and
		// UNLINK where it is visible as it is written. On an object store there
		// is no half-written table to delete -- what exists is an upload in
		// progress, and Remove would delete nothing while leaving its parts
		// behind, invisible to a listing and billed.
		if remote {
			if err := aborter.Abort(m.name); err != nil {
				errs = append(errs, err)
			}
			continue
		}
		// The handle may already be closed by a failed Close; the unlink is the
		// part that matters, so the close error is not interesting here.
		m.sw.Close()
		if err := w.sink.Remove(m.name); err != nil {
			errs = append(errs, fmt.Errorf("removing %s: %w", m.name, err))
		}
	}
	w.tables = nil
	return errors.Join(errs...)
}

// Discard abandons a conversion, leaving nothing behind.
//
// A failure partway through must not leave the tables on disk. They would look
// like a store, they would be truncated or incomplete, and -- because
// conversion refuses to overwrite an existing store -- their presence would
// block the retry.
//
// The parquet writers are deliberately *not* closed first. Closing one writes a
// complete, valid footer into the file, so the old sequence -- footer, footer,
// footer, then three unlinks -- meant a process killed partway through Discard
// left behind exactly the well-formed partial store this exists to prevent.
// There is no reason to finalize a file that is about to be removed.
func (w *Writer) Discard() error {
	return w.abort()
}

// Finish closes the tables and then writes the manifest that marks the store
// complete and readable.
//
// This is the ordinary way to end a conversion; Close alone leaves a store
// without a manifest, which readers refuse. The order is the whole point: every
// table is finalized first, and only a run that got that far writes the marker
// saying so.
// A failure after Close discards the store rather than leaving it. By then the
// three tables carry valid footers, so what is on disk looks exactly like a
// finished conversion minus its marker -- and a reader has no way to tell that
// from an interrupted one. It would be refused at open forever, with no way to
// retry that the overwrite guard did not also have to be talked past. Better to
// leave nothing.
func (w *Writer) Finish() error {
	if err := w.Close(); err != nil {
		return err
	}
	m, err := w.manifest()
	if err != nil {
		return errors.Join(err, w.Discard())
	}
	if err := WriteManifestTo(w.sink, m); err != nil {
		return errors.Join(err, w.Discard())
	}
	return nil
}

// manifest describes the store as written. Table sizes are read back from disk
// rather than accumulated, so they describe the files rather than the intent.
func (w *Writer) manifest() (VolumeManifest, error) {
	tables := map[string]TableInfo{}
	counted := map[string]int64{
		CallsTable:   w.NCalls,
		SitesTable:   w.NSites,
		RegionsTable: w.NRegions,
	}
	// Listed only when the store actually carries one. A coverage entry with
	// zero rows would read as "covered nowhere" to every later query, where no
	// entry at all reads as "nobody said" -- and only the second is true of a
	// conversion that was never given a gVCF.
	if w.opts.Coverage {
		counted[CoverageTable] = w.NCoverage
	}
	for table, rows := range counted {
		// A SPLIT MEMBER IS SIZED SHARD BY SHARD, and each shard's size is read
		// back from the sink rather than accumulated -- the manifest describes
		// the files that exist, not the rows the writer believes it wrote, which
		// is what makes it able to catch a table that does not belong.
		if shards := w.shards[table]; len(shards) > 0 {
			var total int64
			sized := make([]ShardInfo, 0, len(shards))
			for _, si := range shards {
				size, ok, err := w.sink.Stat(si.Name)
				if err != nil {
					return VolumeManifest{}, fmt.Errorf("sizing %s: %w", si.Name, err)
				}
				if !ok {
					return VolumeManifest{}, fmt.Errorf("%s is missing from %s", si.Name, w.sink.Describe())
				}
				si.Bytes = size
				total += size
				sized = append(sized, si)
			}
			tables[table] = TableInfo{Rows: rows, Bytes: total, Shards: sized}
			continue
		}

		size, ok, err := w.sink.Stat(TableFile(table))
		if err != nil {
			return VolumeManifest{}, fmt.Errorf("sizing %s table: %w", table, err)
		}
		if !ok {
			return VolumeManifest{}, fmt.Errorf("the %s table is missing from %s", table, w.sink.Describe())
		}
		tables[table] = TableInfo{Rows: rows, Bytes: size}
	}
	return VolumeManifest{
		FormatVersion: VolumeManifestVersion,
		Complete:      true,
		Created:       time.Now().UTC(),
		Program:       w.opts.Program,
		Command:       w.opts.Command,
		Sources:       w.opts.Sources,
		Meta:          copyMeta(w.opts.Meta),
		Params: ManifestParams{
			MinDP:         w.opts.MinDP,
			NoCallable:    w.opts.NoCallable,
			RowGroupSize:  w.opts.RowGroupSize,
			SpanSemantics: w.opts.Spans,
			DepthBands:    w.opts.DepthBands,
			MaxGap:        w.opts.MaxGap,
			Coverage:      w.opts.Coverage,
			Info:          w.opts.Info,
			Format:        w.opts.Format,
		},
		Samples: w.opts.Samples,
		Counts: ManifestCounts{
			Samples: len(w.opts.Samples),
			Calls:   w.NCalls,
			Sites:   w.NSites,
			Regions: w.NRegions,
		},
		Tables:          tables,
		Chromosomes:     w.chroms,
		ContigsDeclared: w.opts.Contigs,
	}, nil
}

// WriteCall buffers one ALT-carrying genotype.
func (w *Writer) WriteCall(c Call) error {
	return w.WriteCallFormat(c, nil)
}

// WriteCallFormat buffers one ALT call along with the FORMAT values captured
// for it, keyed by the source VCF's own FORMAT key ("PID", not "fmt_pid").
//
// A key with no value for this sample is LEFT OUT rather than zeroed: the
// columns are optional so that "this sample had no PID here" and "its PID is 0"
// stay different claims. The map may be reused between calls.
func (w *Writer) WriteCallFormat(c Call, values map[string]any) error {
	w.shardRows[CallsTable]++
	c.RefEnd = defaultRefEnd(c.Pos, c.Ref, c.RefEnd)
	w.census(c.Chrom, c.Pos).Calls++
	w.NCalls++

	if w.callsAny == nil {
		if len(values) > 0 {
			return fmt.Errorf("call %s at %s:%d carries FORMAT values but the store captures none",
				c.SampleID, c.Chrom, c.Pos)
		}
		w.callBuf = append(w.callBuf, c)
		if len(w.callBuf) >= batchSize {
			return w.flushCalls()
		}
		return nil
	}

	row, err := w.callRow(c, values)
	if err != nil {
		return err
	}
	w.callRows = append(w.callRows, row)
	if len(w.callRows) >= batchSize {
		return w.flushCalls()
	}
	return nil
}

// callRow renders one call as a map matching the dynamic schema, reusing the
// buffer's own slot so a full batch allocates its maps once. Keys absent from
// values are DELETED rather than inherited from the slot's previous occupant.
func (w *Writer) callRow(c Call, values map[string]any) (map[string]any, error) {
	var m map[string]any
	if n := len(w.callRows); cap(w.callRows) > n {
		if prev, ok := w.callRows[:cap(w.callRows)][n].(map[string]any); ok && prev != nil {
			m = prev
			clear(m)
		}
	}
	if m == nil {
		m = make(map[string]any, len(w.opts.Format)+12)
	}
	m["sample_id"] = c.SampleID
	m["chrom"] = c.Chrom
	m["pos"] = c.Pos
	m["ref"] = c.Ref
	m["alt"] = c.Alt
	m["ref_end"] = c.RefEnd
	m["gt"] = c.GT
	m["dp"] = c.DP
	m["ad_ref"] = c.ADRef
	m["ad_alt"] = c.ADAlt
	m["gq"] = c.GQ
	m["min_dp"] = c.MinDP

	slot := len(w.callRows)
	for i, f := range w.opts.Format {
		raw, ok := values[f.Name]
		if !ok {
			continue
		}
		v, err := w.fmtScratch[i].store(
			InfoField{Name: f.Name, Type: f.Type}, slot, raw)
		if err != nil {
			return nil, fmt.Errorf("call %s at %s:%d: %w", c.SampleID, c.Chrom, c.Pos, err)
		}
		m[f.Column] = v
	}
	return m, nil
}

// WriteSite buffers one catalog entry.
func (w *Writer) WriteSite(s Site) error {
	return w.WriteSiteInfo(s, nil)
}

// WriteSiteInfo buffers one catalog entry along with the INFO values captured at
// it, keyed by the source VCF's own INFO key ("R2", not "info_r2").
//
// A key with no value at this site is simply LEFT OUT of the map rather than
// given a zero: the columns are optional precisely so that "the program emitted
// no R2 here" and "R2 is 0 here" stay different claims. Passing a nil map is
// what WriteSite does, and is correct for a store capturing nothing.
//
// The map may be reused between calls -- values are converted on the spot.
func (w *Writer) WriteSiteInfo(s Site, info map[string]any) error {
	s.RefEnd = defaultRefEnd(s.Pos, s.Ref, s.RefEnd)
	// Before anything is buffered, so a rotation lands between sites.
	if err := w.noteSite(s.Chrom, s.Pos); err != nil {
		return err
	}
	w.census(s.Chrom, s.Pos).Sites++
	w.NSites++
	w.shardRows[SitesTable]++

	if w.sitesAny == nil {
		// Captured nothing: values, if any were offered, have nowhere to go.
		// Silently dropping them would be worse than the error, since the store
		// would look like it holds them.
		if len(info) > 0 {
			return fmt.Errorf("site %s:%d carries INFO values but the store captures none", s.Chrom, s.Pos)
		}
		w.siteBuf = append(w.siteBuf, s)
		if len(w.siteBuf) >= batchSize {
			return w.flushSites()
		}
		return nil
	}

	row, err := w.siteRow(s, info)
	if err != nil {
		return err
	}
	w.siteRows = append(w.siteRows, row)
	if len(w.siteRows) >= batchSize {
		return w.flushSites()
	}
	return nil
}

// siteRow renders one site as a map matching the dynamic schema.
//
// The map is taken from the buffer's own slot and refilled, so a full batch
// allocates batchSize maps once for the life of the writer rather than one per
// site. Keys absent from info are DELETED rather than left over from the
// previous occupant of this slot -- otherwise a site would inherit the R2 of
// whatever site used the map 8,192 rows ago, which is the kind of wrong that
// looks entirely plausible.
func (w *Writer) siteRow(s Site, info map[string]any) (map[string]any, error) {
	var m map[string]any
	if n := len(w.siteRows); cap(w.siteRows) > n {
		if prev, ok := w.siteRows[:cap(w.siteRows)][n].(map[string]any); ok && prev != nil {
			m = prev
			clear(m)
		}
	}
	if m == nil {
		m = make(map[string]any, len(w.opts.Info)+8)
	}

	m["chrom"] = s.Chrom
	m["pos"] = s.Pos
	m["ref"] = s.Ref
	m["alt"] = s.Alt
	m["ref_end"] = s.RefEnd
	m["ac"] = s.AC
	m["an"] = s.AN
	m["n_carriers"] = s.NCarriers
	m["n_called"] = s.NCalled
	m["n_lowdp"] = s.NLowDP

	slot := len(w.siteRows)
	for i, f := range w.opts.Info {
		raw, ok := info[f.Name]
		if !ok {
			if f.Type == InfoFlag {
				m[f.Column] = false
			}
			continue
		}
		v, err := w.infoScratch[i].store(f, slot, raw)
		if err != nil {
			return nil, fmt.Errorf("site %s:%d: %w", s.Chrom, s.Pos, err)
		}
		m[f.Column] = v
	}
	return m, nil
}

// infoSlots is one captured field's backing storage for a batch of rows.
//
// Only the slice its declared type needs is allocated; the others stay nil.
type infoSlots struct {
	f64 []float64
	i32 []int32
	str []string
}

func newInfoSlots(fields []InfoField, n int) []infoSlots {
	if len(fields) == 0 {
		return nil
	}
	out := make([]infoSlots, len(fields))
	for i, f := range fields {
		switch f.Type {
		case InfoFloat:
			out[i].f64 = make([]float64, n)
		case InfoInteger:
			out[i].i32 = make([]int32, n)
		case InfoString:
			out[i].str = make([]string, n)
		}
	}
	return out
}

// store puts a value in this row's slot and returns what the parquet writer
// should receive for it.
//
// A POINTER for the optional columns, and the reason is not style: parquet-go
// encodes an optional field's zero value as null, so handing it a plain
// float64(0) would record a genuine R2 of 0 as "no R2 at this site". Those are
// different claims -- "worthless" and "unknown" -- and a filter for R2 >= 0.3
// only treats them alike by accident. A Flag is required rather than optional,
// where false is a real value and needs no pointer.
//
// Deliberately narrow about types: a Float column takes a float, not a string
// that looks like one. Parsing belongs to whoever read the VCF header and knows
// the declared type; accepting both here would let a caller's mistake become a
// column of zeros nobody questions.
func (s *infoSlots) store(f InfoField, slot int, raw any) (any, error) {
	switch f.Type {
	case InfoInteger:
		var n int32
		switch x := raw.(type) {
		case int32:
			n = x
		case int:
			n = int32(x)
		case int64:
			n = int32(x)
		default:
			return nil, fmt.Errorf("info %s declared Integer, given %T", f.Name, raw)
		}
		s.i32[slot] = n
		return &s.i32[slot], nil

	case InfoFloat:
		var x float64
		switch v := raw.(type) {
		case float64:
			x = v
		case float32:
			x = float64(v)
		default:
			return nil, fmt.Errorf("info %s declared Float, given %T", f.Name, raw)
		}
		s.f64[slot] = x
		return &s.f64[slot], nil

	case InfoFlag:
		b, ok := raw.(bool)
		if !ok {
			return nil, fmt.Errorf("info %s declared Flag, given %T", f.Name, raw)
		}
		return b, nil

	default:
		v, ok := raw.(string)
		if !ok {
			return nil, fmt.Errorf("info %s declared String, given %T", f.Name, raw)
		}
		s.str[slot] = v
		return &s.str[slot], nil
	}
}

// defaultRefEnd fills in a reference end the caller left unset.
//
// len(REF) is recoverable from the row itself, so deriving it here means every
// writer gets overlap for plain records without having to think about it, and the
// column is never a field of zeros that makes HasRefSpans lie. A caller that knows
// better -- a symbolic ALT with INFO/END, a gVCF block -- passes its own wider
// value and this leaves it alone.
func defaultRefEnd(pos int32, ref string, refEnd int32) int32 {
	derived := pos - 1 + int32(len(ref))
	if derived <= pos-1 {
		derived = pos
	}
	if refEnd > derived {
		return refEnd
	}
	return derived
}

// WriteRegion buffers one callable run.
func (w *Writer) WriteRegion(r CalledSiteRun) error {
	// A RUN MUST LIE WHOLLY INSIDE ONE SHARD, and a caller that lets one span a
	// boundary is refused rather than obeyed.
	//
	// The reader prunes region shards by coordinate, so a run filed in shard 3
	// but starting in shard 1 is a run no query for shard 1 will ever see --
	// and every locus it covered there reads as NEVER ASSAYED instead of
	// reference. That is silent, plausible, and moves a denominator, which is
	// exactly the failure this package is arranged against. Callers break their
	// runs in BeforeRotate; this is what catches one that forgot.
	if w.opts.ShardSites > 0 && w.shardSites > 0 && SameChrom(r.Chrom, w.shardChrom) &&
		(r.Start < w.shardFirst || r.Start > w.shardLast) {
		// BOTH DIRECTIONS, because the two mistakes look nothing alike and both
		// are silent. A run beginning EARLIER than the shard was never broken at
		// the boundary. One beginning LATER was created after the boundary
		// passed and filed in the shard that is closing -- which is what happens
		// when a caller extends a run into the next shard's first site and only
		// then writes that site, letting the rotation discover itself too late.
		//
		// Either way the run lands in a shard whose range does not contain it,
		// so no query for those positions will ever see it and every locus it
		// covered reads as never assayed.
		return fmt.Errorf(
			"run %s %s:%d-%d lies outside the open shard (%s:%d-%d): ask WouldRotate before "+
				"extending a run into a new shard, close the runs, then Rotate -- otherwise every "+
				"locus this run covers reads as never assayed",
			r.SampleID, r.Chrom, r.Start, r.End, w.shardChrom, w.shardFirst, w.shardLast)
	}
	w.shardRows[RegionsTable]++
	w.regionBuf = append(w.regionBuf, r)
	w.NRegions++
	if len(w.regionBuf) >= batchSize {
		return w.flushRegions()
	}
	return nil
}

func (w *Writer) flushCalls() error {
	if w.callsAny != nil {
		if len(w.callRows) == 0 {
			return nil
		}
		if _, err := w.callsAny.Write(w.callRows); err != nil {
			return fmt.Errorf("writing calls: %w", err)
		}
		// Length only: the maps beyond it stay allocated for callRow to reclaim.
		w.callRows = w.callRows[:0]
		return nil
	}
	if len(w.callBuf) == 0 {
		return nil
	}
	if _, err := w.calls.Write(w.callBuf); err != nil {
		return fmt.Errorf("writing calls: %w", err)
	}
	w.callBuf = w.callBuf[:0]
	return nil
}

func (w *Writer) flushSites() error {
	if w.sitesAny != nil {
		if len(w.siteRows) == 0 {
			return nil
		}
		if _, err := w.sitesAny.Write(w.siteRows); err != nil {
			return fmt.Errorf("writing sites: %w", err)
		}
		// Length only: the maps beyond it stay allocated for siteRow to reclaim.
		w.siteRows = w.siteRows[:0]
		return nil
	}
	if len(w.siteBuf) == 0 {
		return nil
	}
	if _, err := w.sites.Write(w.siteBuf); err != nil {
		return fmt.Errorf("writing sites: %w", err)
	}
	w.siteBuf = w.siteBuf[:0]
	return nil
}

func (w *Writer) flushRegions() error {
	if len(w.regionBuf) == 0 {
		return nil
	}
	if _, err := w.regions.Write(w.regionBuf); err != nil {
		return fmt.Errorf("writing regions: %w", err)
	}
	w.regionBuf = w.regionBuf[:0]
	return nil
}

// WriteCoverage buffers one coverage block.
//
// Refused unless the store was opened with Coverage set, because a block
// arriving at a store that will not record one is a caller who believes they
// are capturing something they are not -- and the loss is unrecoverable, since
// the gVCF that could have answered is not read again.
func (w *Writer) WriteCoverage(b CoverageBlock) error {
	if !w.opts.Coverage {
		return fmt.Errorf(
			"this store was not opened to hold coverage blocks; set WriterOpts.Coverage, "+
				"or the span %s %s:%d-%d is discarded and cannot be recovered without the source",
			b.SampleID, b.Chrom, b.Start, b.End)
	}
	if b.End < b.Start {
		return fmt.Errorf("coverage block %s %s:%d-%d ends before it starts",
			b.SampleID, b.Chrom, b.Start, b.End)
	}
	w.shardRows[CoverageTable]++
	w.coverageBuf = append(w.coverageBuf, b)
	w.NCoverage++
	if len(w.coverageBuf) >= batchSize {
		return w.flushCoverage()
	}
	return nil
}

func (w *Writer) flushCoverage() error {
	if w.coverage == nil || len(w.coverageBuf) == 0 {
		return nil
	}
	if _, err := w.coverage.Write(w.coverageBuf); err != nil {
		return fmt.Errorf("writing coverage: %w", err)
	}
	w.coverageBuf = w.coverageBuf[:0]
	return nil
}

// Close flushes and finalizes all three files.
//
// It stops at the first failure instead of pressing on. Continuing would finish
// the other tables, and a parquet footer is written precisely by the Close that
// would then run -- so a failed flush used to yield three structurally valid
// files of which one was silently short by up to a batch. Three well-formed
// tables that disagree is the one outcome a reader cannot detect, so the write
// stops as soon as it cannot be completed honestly.
//
// A non-nil error means nothing on disk should be trusted; the caller is
// expected to Discard. The file handles are released either way, since leaking
// them helps nobody, but the tables are left in place for Discard to remove.
func (w *Writer) Close() error {
	// The last shard is closed the same way every other one was, so its bounds
	// and row counts are recorded by the same code -- a final shard described
	// differently from its siblings is a shard a reader would prune wrongly.
	if err := w.closeShard(); err != nil {
		w.closeFiles()
		return err
	}
	return w.closeFiles()
}

// closeFiles releases the underlying handles without removing anything.
func (w *Writer) closeFiles() error {
	var errs []error
	for i := range w.tables {
		if w.tables[i].closed {
			continue
		}
		if err := w.tables[i].sw.Close(); err != nil {
			errs = append(errs, fmt.Errorf("finishing %s: %w", w.tables[i].name, err))
		}
		w.tables[i].closed = true
	}
	return errors.Join(errs...)
}

// scanParquet streams path in batches, calling fn per row. fn returns false to
// stop early. Batching keeps a whole-genome store from having to be resident.
func scanParquet[T any](set *shardSet, fn func(T) bool) error {
	return scanParquetPruned(set, scanFilter{shard: keepAllShards, rowGroup: keepAll}, fn)
}

// scanWholeTable streams one file, for callers holding a table directly.
func scanWholeTable[T any](m *table, fn func(T) bool) error {
	pf, err := m.parsed()
	if err != nil {
		return err
	}
	// The parsed file, not the raw ReaderAt: see scanParquetPruned.
	r := parquet.NewGenericReader[T](pf)
	defer r.Close()

	buf := make([]T, 1024)
	for {
		n, err := r.Read(buf)
		for i := 0; i < n; i++ {
			if !fn(buf[i]) {
				return nil
			}
		}
		if err == io.EOF {
			return nil
		}
		if err != nil {
			return err
		}
		if n == 0 {
			return nil
		}
	}
}

// hasColumn reports whether a leaf column is present in a file's schema.
func hasColumn(pf *parquet.File, name string) bool {
	for _, f := range pf.Schema().Fields() {
		if f.Name() == name {
			return true
		}
	}
	return false
}

// ParquetVolume is a Store backed by the three-file Parquet set.
type ParquetVolume struct {
	base        string
	calls       *shardSet
	hasRefSpans bool
	sites       *shardSet
	regions     *shardSet
	coverage    *shardSet
	samples     []string
	hasSites    bool
	hasRegions  bool

	// hasCoverage distinguishes a store that was given genomic block spans from
	// one that was not, which is the difference between "not covered there" and
	// "nobody said". Only the first may be reported as an answer.
	hasCoverage bool
	noCallable  bool
	spans       SpanSemantics
	minDP       int32
	meta        map[string]string
	manifest    *VolumeManifest
}

// HasCoverage reports whether this store carries genomic block spans, and so
// whether it can answer about positions absent from the sites catalog.
func (s *ParquetVolume) HasCoverage() bool { return s.hasCoverage }

// MaxGap is the largest uncovered stretch this store's coverage blocks were
// permitted to span. Zero when there are none, or when they are unmerged.
func (s *ParquetVolume) MaxGap() int32 {
	if s.manifest == nil {
		return 0
	}
	return s.manifest.Params.MaxGap
}

// Coverage streams the coverage blocks, calling fn for each until it returns
// false. Reports nothing at all for a store that carries none.
func (s *ParquetVolume) Coverage(fn func(CoverageBlock) bool) error {
	if !s.hasCoverage {
		return nil
	}
	return scanParquet(s.coverage, fn)
}

// SpanSemantics reports what this store's run intervals may claim. A store
// written from a plain VCF reports SpansSites, confining answers to the sites
// catalog.
func (s *ParquetVolume) SpanSemantics() SpanSemantics { return s.spans }

// Provenance is what a store records about how it was built. It matters at
// query time chiefly because of MinDP: a store baked that threshold into its
// called-site runs, so a query gating at a different value is not asking the
// same question the store can answer.
type Provenance struct {
	Source     string
	Program    string
	Command    string
	MinDP      int32
	NoCallable bool
	Spans      SpanSemantics

	NumSamples int

	// Meta is the caller-supplied metadata recorded at conversion, read back
	// from the calls file rather than the manifest -- so a query-time caller
	// that already has the store open need not open a second file to learn what
	// dataset it is holding. Nil when none was recorded.
	Meta map[string]string
}

// Provenance returns the conversion metadata recorded in the calls file.
func (s *ParquetVolume) Provenance() Provenance {
	return Provenance{
		Source:     s.meta[MetaSource],
		Program:    s.meta[MetaProgram],
		Command:    s.meta[MetaCommand],
		MinDP:      s.minDP,
		NoCallable: s.noCallable,
		Spans:      s.spans,
		NumSamples: len(s.samples),
		Meta:       s.metaFields(),
	}
}

// metaFields recovers the caller-supplied metadata from the calls file's
// key/value metadata by stripping MetaPrefix. Returns nil rather than an empty
// map when none was recorded, so a caller can tell "not stated" from "empty".
func (s *ParquetVolume) metaFields() map[string]string {
	var out map[string]string
	for k, v := range s.meta {
		key, ok := strings.CutPrefix(k, MetaPrefix)
		if !ok {
			continue
		}
		if out == nil {
			out = map[string]string{}
		}
		out[key] = v
	}
	return out
}

// OpenParquet opens a Parquet store. base may be given either as the base name
// or as the path to any one of the three files.
func OpenParquet(base string) (*ParquetVolume, error) {
	return OpenParquetContext(context.Background(), base)
}

// OpenParquetContext opens a store from any locator: a filesystem path, an
// http(s):// URL, or any scheme registered with iosource, such as s3://.
//
// Parquet suits remote reading unusually well. The footer carries per-row-group
// statistics, and this package already prunes groups by them, so a locus query
// against a remote store skips the pruned groups without transferring them at
// all — the same mechanism that makes it fast locally makes it cheap remotely.
//
// The tables stay open for the life of the store, so Close matters more than
// it used to: it now releases them rather than being a no-op.
func OpenParquetContext(ctx context.Context, base string) (*ParquetVolume, error) {
	base = TrimStoreSuffix(base)
	man, err := requireManifest(ctx, base)
	if err != nil {
		return nil, err
	}
	calls, err := openShardSet(ctx, base, CallsTable, man.Tables[CallsTable])
	if err != nil {
		return nil, fmt.Errorf("opening parquet store %s: %w", base, err)
	}
	// The key/value metadata lives on the first shard: a split table writes it
	// once rather than repeating it, and every shard of one store agrees about
	// what the store is.
	first := calls.single()
	pf, err := parquet.OpenFile(first.ra, first.size)
	if err != nil {
		calls.Close()
		return nil, fmt.Errorf("reading %s: %w", first.name, err)
	}
	s := &ParquetVolume{base: base, calls: calls, meta: map[string]string{}}
	for _, k := range []string{MetaSource, MetaProgram, MetaCommand, MetaMinDP, MetaContigs} {
		if v, ok := pf.Lookup(k); ok {
			s.meta[k] = v
		}
	}
	// Caller metadata is open-ended, so it cannot be read by an allowlist the way
	// the keys above are; sweep the file's metadata for the namespace instead.
	for _, kv := range pf.Metadata().KeyValueMetadata {
		if strings.HasPrefix(kv.Key, MetaPrefix) {
			s.meta[kv.Key] = kv.Value
		}
	}
	if v, err := strconv.Atoi(s.meta[MetaMinDP]); err == nil {
		s.minDP = int32(v)
	}
	if v, ok := pf.Lookup(MetaSamples); ok && v != "" {
		s.samples = strings.Split(v, "\n")
	}
	if v, ok := pf.Lookup(MetaNoCall); ok && v == "1" {
		s.noCallable = true
	}
	s.spans = SpansSites // the conservative reading for stores predating the key
	if v, ok := pf.Lookup(MetaSpans); ok && v != "" {
		s.spans = SpanSemantics(v)
	}
	// Whether this store can answer overlap queries is read from the schema rather
	// than from a metadata key: the column either exists or it does not, and a key
	// could disagree with the data. Absent means the store predates ref_end, so
	// span selection there matches on position alone -- a different answer, not
	// merely a slower one, which is why HasRefSpans is exported for callers to
	// report.
	s.hasRefSpans = hasColumn(pf, "ref_end")
	// The optional tables. Absence is legitimate only where the manifest
	// recorded nothing in them, which verifyAgainstManifest below enforces;
	// note that the writer creates all three regardless, so --no-callable
	// produces a present, zero-row regions table rather than a missing one.
	sites, hasSites, err := openOptionalShardSet(ctx, base, SitesTable, man.Tables[SitesTable])
	if err != nil {
		calls.Close()
		return nil, err
	}
	s.sites, s.hasSites = sites, hasSites
	regions, hasRegions, err := openOptionalShardSet(ctx, base, RegionsTable, man.Tables[RegionsTable])
	if err != nil {
		calls.Close()
		sites.Close()
		return nil, err
	}
	s.regions, s.hasRegions = regions, hasRegions

	// OPTIONAL, and absent is the common case: only a conversion given a gVCF
	// or a callable mask has one. Absent means nothing was claimed off the
	// catalog, never that nothing was covered.
	coverage, hasCoverage, err := openOptionalShardSet(ctx, base, CoverageTable, man.Tables[CoverageTable])
	if err != nil {
		calls.Close()
		sites.Close()
		regions.Close()
		return nil, err
	}
	s.coverage, s.hasCoverage = coverage, hasCoverage
	s.manifest = man

	if err := verifyAgainstManifest(man, pf, calls, sites, regions); err != nil {
		s.Close()
		return nil, fmt.Errorf("%s: %w", base, err)
	}
	return s, nil
}

// requireManifest loads the completion marker, refusing the store without it.
//
// A store with no manifest was not finished -- or predates them -- and either
// way nothing here can tell how much of the intended input it holds. The
// tables would open and answer queries perfectly well; that is the problem.
// An unfinished store reports "not assayed" for everything it never got to,
// which is indistinguishable from the honest answer for a position the source
// genuinely never reported.
func requireManifest(ctx context.Context, base string) (*VolumeManifest, error) {
	man, err := ReadVolumeManifestContext(ctx, base)
	if err != nil {
		return nil, fmt.Errorf("%s has no readable %s, so it cannot be shown to be "+
			"a complete store: %w\n\tit is from an interrupted conversion, or predates "+
			"manifests; re-convert it, or inspect it with vcf-varsummary",
			base, VolumeManifestFile, err)
	}
	if !man.Complete {
		return nil, fmt.Errorf("%s is marked incomplete by its %s; re-convert it",
			base, VolumeManifestFile)
	}
	if man.FormatVersion > VolumeManifestVersion {
		return nil, fmt.Errorf("%s was written by a newer cgkit (manifest format %d, "+
			"this build understands %d); upgrade to read it",
			base, man.FormatVersion, VolumeManifestVersion)
	}
	return man, nil
}

// verifyAgainstManifest checks each table against the row count recorded for
// it when the store was written.
//
// The footers say each table was finished; the manifest says which tables
// finishing this store produced. Comparing them is what catches a table that
// is well-formed but does not belong -- one copied in from another conversion,
// or left over from a run whose replacement died. Nothing else can see that,
// because sites and regions carry no metadata of their own.
//
// It costs nothing: every count here is footer metadata that was already read.
func verifyAgainstManifest(man *VolumeManifest, callsFile *parquet.File, calls, sites, regions *shardSet) error {
	type check struct {
		name string
		file *parquet.File
		set  *shardSet
	}
	for _, c := range []check{
		// The parsed file is the FIRST SHARD when calls are split, so it is
		// passed alongside the set rather than instead of it -- checking it
		// against the table's total would compare one shard's rows with all of
		// them, which is how this first failed.
		{CallsTable, callsFile, calls},
		{SitesTable, nil, sites},
		{RegionsTable, nil, regions},
	} {
		want, recorded := man.Tables[c.name]
		if !recorded {
			continue // a manifest that never described this table claims nothing
		}

		// A SPLIT MEMBER IS COUNTED ACROSS ITS SHARDS, and every shard's own
		// footer is checked against the count the index recorded for it. That
		// is strictly stronger than the whole-table total: a store missing one
		// shard and gaining rows in another would still total correctly.
		if c.set != nil && c.set.split() {
			c.file = nil
			var total int64
			for _, sh := range c.set.shards {
				f, err := parquet.OpenFile(sh.m.ra, sh.m.size)
				if err != nil {
					return fmt.Errorf("reading %s: %w", sh.info.Name, err)
				}
				got := f.NumRows()
				if got != sh.info.Rows {
					return fmt.Errorf("%s holds %d rows but the manifest records %d; "+
						"this shard does not belong to this store, or the store is damaged",
						sh.info.Name, got, sh.info.Rows)
				}
				total += got
			}
			if total != want.Rows {
				return fmt.Errorf("%s totals %d rows across its shards but the manifest records %d",
					c.name, total, want.Rows)
			}
			continue
		}

		if c.file == nil {
			m := c.set.single()
			if m == nil {
				if want.Rows == 0 {
					continue // written empty and since removed; it claimed nothing
				}
				return fmt.Errorf("%s.parquet is missing, but the manifest records %d rows in it",
					c.name, want.Rows)
			}
			f, err := parquet.OpenFile(m.ra, m.size)
			if err != nil {
				return fmt.Errorf("reading %s.parquet: %w", c.name, err)
			}
			c.file = f
		}
		if got := c.file.NumRows(); got != want.Rows {
			return fmt.Errorf("%s.parquet holds %d rows but the manifest records %d; "+
				"this table does not belong to this store, or the store is damaged",
				c.name, got, want.Rows)
		}
	}
	return nil
}

// openOptionalTable opens a table a store may legitimately lack, and verifies
// that what it finds is a readable parquet file.
//
// The two outcomes have to stay distinct. Absence is normal: a --no-callable
// store records no callable runs. A table that is *present but unreadable* is
// not, and until this checked, the two were indistinguishable -- the open error
// was discarded either way, so a truncated or half-written sites file read
// exactly like a deliberate omission and the store went on to answer queries as
// though those runs had never been asked for.
//
// Parsing the footer here is what makes the check meaningful, and it is nearly
// free: it is footer metadata, not a data scan, and every query would have
// parsed it anyway. A parquet footer is written only by the writer's Close, so
// a table that parses is a table that was finished.
//
// One limit worth naming: for a remote locator "absent" and "unreachable" are
// not cleanly separable here, because a 404 surfaces as a failed size probe
// rather than as fs.ErrNotExist. Such a table is treated as absent, which is
// the pre-existing behaviour. The manifest is what makes this exact, since it
// records which tables the conversion actually wrote.
func openOptionalTable(ctx context.Context, locator string) (*table, bool, error) {
	m, err := openTable(ctx, locator)
	if err != nil {
		// Absent is an answer; anything else is a failure. This used to discard
		// every error alike, so a table that existed but could not be read --
		// bad permissions, an I/O error, a symlink loop -- was reported as one
		// that was never written, and the store opened with a quietly missing
		// half. Remote absence reads the same way as local now: iosource wraps
		// a 404 in fs.ErrNotExist.
		if errors.Is(err, fs.ErrNotExist) {
			return nil, false, nil
		}
		return nil, false, fmt.Errorf("opening %s: %w", locator, err)
	}
	if _, err := parquet.OpenFile(m.ra, m.size); err != nil {
		m.Close()
		return nil, false, fmt.Errorf("reading %s: %w", locator, err)
	}
	return m, true, nil
}

// TrimStoreSuffix reduces any spelling of a store to its directory.
//
//	cohort                    the store
//	cohort/                   the same store; the separator is optional
//	cohort/calls.parquet      a table of it
//	cohort/manifest.json.gz   likewise
//
// Pointing at a table is worth accepting because tab completion lands there,
// and because a locator naming a file is unambiguous in a way a bare directory
// name is not. The manifest spelling is what makes `stores/*/manifest.json.gz`
// a useful glob: the shell filters to the completed stores without opening
// anything.
//
// This no longer consults the filesystem. It used to, to tell a directory base
// from a filename prefix -- a question that only existed while both forms did,
// and one that could not be answered at all for a remote locator, where the
// os.Stat simply failed and the path fell through unchanged.
func TrimStoreSuffix(p string) string {
	for _, name := range []string{
		CallsTable + ".parquet", SitesTable + ".parquet", RegionsTable + ".parquet",
		VolumeManifestFile,
	} {
		if base, ok := cutTableSuffix(p, name); ok {
			return base
		}
	}
	return trimStoreDir(p)
}

// cutTableSuffix removes a trailing "/name" (or the platform separator) from p,
// reporting whether it was there. A bare "name" with no separator is a table of
// the current directory, whose store base is "".
func cutTableSuffix(p, name string) (string, bool) {
	if p == name {
		return "", true
	}
	for _, sep := range []string{"/", string(os.PathSeparator)} {
		if base, ok := strings.CutSuffix(p, sep+name); ok {
			return trimStoreDir(base), true
		}
	}
	return "", false
}

func fileExists(p string) bool {
	st, err := os.Stat(p)
	return err == nil && !st.IsDir()
}

// Contigs returns the source's ##contig header lines, or nothing for a store
// written before they were recorded. A caller exporting to VCF should emit these
// rather than synthesizing lines from the loci it happens to have seen.
func (s *ParquetVolume) Contigs() []string {
	v := s.meta[MetaContigs]
	if v == "" {
		return nil
	}
	return strings.Split(v, "\n")
}

// Samples returns the roster recorded at conversion time.
func (s *ParquetVolume) Samples() ([]string, error) {
	if len(s.samples) == 0 {
		return nil, fmt.Errorf("%s records no sample list", CallsPath(s.base))
	}
	return s.samples, nil
}

// Calls streams the genotypes a query selects, in the store's own order.
func (s *ParquetVolume) Calls(q Query) (iter.Seq2[Call, error], error) {
	// Fail before iterating: a caller asking for reference calls from a store that
	// cannot reconstruct them should learn so now, not receive a silently
	// ALT-only stream.
	if q.IncludeRef {
		if err := s.classifiable(); err != nil {
			return nil, err
		}
	}
	p := q.plan()
	keep := callsFilter(q)

	if !q.IncludeRef {
		// The ALT-only case genuinely streams. calls.parquet is written in store
		// order, so rows are yielded as they decode with nothing buffered.
		return func(yield func(Call, error) bool) {
			err := scanParquetPruned(s.calls, keep, func(c Call) bool {
				if !p.wantsCall(c) {
					return true
				}
				return yield(c, nil)
			})
			if err != nil {
				yield(Call{}, err)
			}
		}, nil
	}

	// Reference calls have to be interleaved with the ALT calls in one order, and
	// they are derived from two other files, so this path assembles before it
	// yields. Ordering comes from walking sites.parquet -- itself in store order --
	// rather than from sorting, so contig order is preserved for free.
	//
	// It is bounded by the query, not by the store: naming samples or loci bounds
	// what is held. An unrestricted IncludeRef query over a whole store is the case
	// that still wants a streaming merge of the two scans.
	return func(yield func(Call, error) bool) {
		rows, err := s.callsWithRef(p, keep)
		if err != nil {
			yield(Call{}, err)
			return
		}
		for _, c := range rows {
			if !yield(c, nil) {
				return
			}
		}
	}, nil
}

// callsWithRef assembles ALT and reference calls for a query, in store order.
func (s *ParquetVolume) callsWithRef(p *plan, keep scanFilter) ([]Call, error) {
	roster, err := s.Samples()
	if err != nil {
		return nil, err
	}
	wanted := make([]string, 0, len(roster))
	for _, name := range roster {
		if p.wantsSample(name) {
			wanted = append(wanted, name)
		}
	}

	// One pass over the calls. Two things come out of it: the ALT rows to emit,
	// grouped by site, and -- keyed by record rather than locus, and deliberately
	// UNGATED -- which samples carry some alternate. A below-gate ALT makes a
	// sample an uncertain carrier, never a reference observation, and any alternate
	// of a multiallelic record disqualifies every site split out of it.
	altBySite := map[Locus][]Call{}
	carries := map[string]map[RecordKey]bool{}
	if err := scanParquetPruned(s.calls, keep, func(c Call) bool {
		if !p.wantsSample(c.SampleID) {
			return true
		}
		loc := c.Locus()
		// Record level, not site level: a call at a SIBLING locus of the same
		// record still disqualifies the sample from being reference here, and
		// testing wantsSite would throw that evidence away.
		if !p.wantsRecord(loc, c.RefEnd) {
			return true
		}
		if carries[c.SampleID] == nil {
			carries[c.SampleID] = map[RecordKey]bool{}
		}
		carries[c.SampleID][loc.Record()] = true
		if p.wantsSite(loc, c.RefEnd) && p.q.Gate.Admits(c) {
			k := canonLocus(loc)
			altBySite[k] = append(altBySite[k], c)
		}
		return true
	}); err != nil {
		return nil, err
	}

	runs, err := s.runsFor(wanted, p)
	if err != nil {
		return nil, err
	}

	var out []Call
	err = scanParquetPruned(s.sites, callsFilter(p.q), func(site Site) bool {
		loc := site.Locus()
		if !p.wantsSite(loc, site.RefEnd) {
			return true
		}
		k := canonLocus(loc)
		// One pass over the roster, emitting each sample's ALT call or its reference
		// call. Interleaved rather than ALT-block-then-reference-block, so rows are
		// genuinely ordered by sample within a locus -- which is what the Calls
		// contract promises. Emitting the blocks separately put an ALT carrier ahead
		// of a lower-numbered sample's reference call.
		alt := make(map[string]Call, len(altBySite[k]))
		for _, c := range altBySite[k] {
			alt[c.SampleID] = c
		}
		delete(altBySite, k)
		for _, name := range wanted {
			if c, ok := alt[name]; ok {
				out = append(out, c)
				continue
			}
			if carries[name][loc.Record()] {
				continue // carries this or another alternate of the record
			}
			if !runs.covers(name, site.Chrom, site.Pos) {
				continue // never called here at adequate depth
			}
			out = append(out, HomRefCall(name, loc))
		}
		return true
	})
	if err != nil {
		return nil, err
	}
	// A call whose site is missing from the catalog should be impossible, but
	// dropping it silently would be worse than emitting it late. Sorted rather
	// than drained in map order: Go randomizes that, so the same query against
	// the same store returned these rows in a different order on every run, and
	// a caller diffing two runs saw a difference that was not in the data.
	if len(altBySite) > 0 {
		keys := make([]Locus, 0, len(altBySite))
		for k := range altBySite {
			keys = append(keys, k)
		}
		// Not the store's contig order -- these rows are by definition absent
		// from the catalog that defines it. Any total order will do; what
		// matters is that it is the same one twice.
		sort.Slice(keys, func(i, j int) bool {
			a, b := keys[i], keys[j]
			switch {
			case a.Chrom != b.Chrom:
				return a.Chrom < b.Chrom
			case a.Pos != b.Pos:
				return a.Pos < b.Pos
			case a.Ref != b.Ref:
				return a.Ref < b.Ref
			default:
				return a.Alt < b.Alt
			}
		})
		for _, k := range keys {
			calls := altBySite[k]
			sort.SliceStable(calls, func(i, j int) bool {
				return calls[i].SampleID < calls[j].SampleID
			})
			out = append(out, calls...)
		}
	}
	return out, nil
}

// Classify resolves every sample at a locus.
//
// It refuses rather than guesses: without the sites catalog there is no way to
// know the position was interrogated at all, and without run information no way
// to know a silent sample was called there. Returning "non-carrier" in either
// case would be a fabricated observation.
//
// The sites catalog is a hard gate. For a store whose spans are SpansSites --
// everything a plain VCF can produce -- a locus absent from the catalog is
// StateNotAssayed for every sample, no matter that run intervals appear to
// bracket it. Those intervals only mark catalog sites; the bases between them
// were never reported, and treating a run as coverage would invent reference
// observations. Only a gVCF-derived store (SpansBlocks) could answer here.
func (s *ParquetVolume) Classify(l Locus, g Gate) ([]SampleState, error) {
	if err := s.classifiable(); err != nil {
		return nil, err
	}
	samples, err := s.Samples()
	if err != nil {
		return nil, err
	}

	// Was this position interrogated at all? This is a hard gate rather than one
	// input among several: if the source never reported the locus, the run
	// intervals must not be consulted, so return before they are even opened.
	interrogated, err := s.SiteKnown(l)
	if err != nil {
		return nil, err
	}
	if !interrogated && s.spans != SpansBlocks {
		out := make([]SampleState, 0, len(samples))
		for _, name := range samples {
			out = append(out, SampleState{SampleID: name, State: StateNotAssayed})
		}
		return out, nil
	}

	calls := map[string]Call{}
	if err := scanParquetPruned(s.calls, locusFilter(l), func(c Call) bool {
		if SameLocus(c.Locus(), l) {
			calls[c.SampleID] = c
		}
		return true
	}); err != nil {
		return nil, err
	}

	// Reached only for a locus in the catalog (or a block-semantics store), so a
	// run bracketing the position genuinely means "called here".
	called := map[string]bool{}
	if err := scanParquetPruned(s.regions, coveringFilter(l.Chrom, l.Pos), func(r CalledSiteRun) bool {
		if SameChrom(r.Chrom, l.Chrom) && l.Pos >= r.Start && l.Pos <= r.End {
			called[r.SampleID] = true
		}
		return true
	}); err != nil {
		return nil, err
	}

	out := make([]SampleState, 0, len(samples))
	for _, name := range samples {
		st := SampleState{SampleID: name}
		if c, ok := calls[name]; ok {
			cc := c
			st.Call = &cc
			if g.Admits(c) {
				st.State = StateCarrier
			} else {
				st.State = StateUncertain
			}
		} else if called[name] {
			st.State = StateNonCarrier
		} else {
			st.State = StateNotAssayed
		}
		out = append(out, st)
	}
	return out, nil
}

// classifiable reports whether the store carries the evidence needed to tell a
// reference call apart from a position that was never assayed. Every query that
// would otherwise have to guess goes through this first.
func (s *ParquetVolume) classifiable() error {
	if !s.hasSites {
		return fmt.Errorf("%w: %s is missing", ErrNotClassifiable, SitesPath(s.base))
	}
	if !s.hasRegions {
		return fmt.Errorf("%w: %s is missing", ErrNotClassifiable, RegionsPath(s.base))
	}
	if s.noCallable {
		return fmt.Errorf("%w: %s was built with --no-callable (the source had no DP field)",
			ErrNotClassifiable, s.base)
	}
	return nil
}

// runsBySample holds callable runs for the samples a query selected, grouped by
// canonical chromosome and sorted by start.
//
// Held in memory because reconstructing reference calls tests every catalog site
// against it, and re-scanning the regions file per site would be quadratic. A run
// spans as many sites as the sample was called at consecutively, so this stays far
// smaller than the regions file -- but it grows with the number of samples
// selected, which is the bound on an unrestricted IncludeRef query.
type runsBySample map[string]map[string][]CalledSiteRun

// covers reports whether a run brackets a 1-based position for one sample. Runs
// for one sample on one chromosome are disjoint and increasing -- the converter
// closes a run before opening the next -- so the last run starting at or before
// pos is the only one that can contain it.
func (r runsBySample) covers(sample, chrom string, pos int32) bool {
	runs := r[sample][CanonKey(chrom)]
	i := sort.Search(len(runs), func(i int) bool { return runs[i].Start > pos })
	if i == 0 {
		return false
	}
	return pos <= runs[i-1].End
}

// runsFor loads the called-site runs for the named samples in one pass.
//
// The query only prunes row groups here; runs surviving that are kept whatever
// their extent, since a run beginning before a span may still reach into it.
func (s *ParquetVolume) runsFor(samples []string, p *plan) (runsBySample, error) {
	want := make(map[string]bool, len(samples))
	for _, name := range samples {
		want[name] = true
	}
	out := runsBySample{}
	err := scanParquetPruned(s.regions, runScanFilter(p.q), func(r CalledSiteRun) bool {
		if !want[r.SampleID] {
			return true
		}
		if out[r.SampleID] == nil {
			out[r.SampleID] = map[string][]CalledSiteRun{}
		}
		k := CanonKey(r.Chrom)
		out[r.SampleID][k] = append(out[r.SampleID][k], r)
		return true
	})
	if err != nil {
		return nil, err
	}
	for _, byChrom := range out {
		for _, runs := range byChrom {
			sort.Slice(runs, func(i, j int) bool { return runs[i].Start < runs[j].Start })
		}
	}
	return out, nil
}

// Sites streams the catalog, calling fn per site. fn returns false to stop.
func (s *ParquetVolume) Sites(fn func(Site) bool) error {
	if !s.hasSites {
		return fmt.Errorf("%s is missing", SitesPath(s.base))
	}
	return scanParquet(s.sites, fn)
}

// Regions walks the callable runs.
//
// The counterpart to Sites, and the way a consumer inspects what a reference
// call rests on: a run carries the span, the site count and -- since depth
// banding -- the lowest depth vouched for across the whole of it.
func (s *ParquetVolume) Regions(fn func(CalledSiteRun) bool) error {
	if !s.hasRegions {
		return fmt.Errorf("%s is missing", RegionsPath(s.base))
	}
	return scanParquet(s.regions, fn)
}

// FormatFields reports the FORMAT fields this store captured onto its calls,
// nil if none.
//
// From the manifest rather than the schema, for the same reason InfoFields is:
// the schema says which columns exist, the manifest says what they MEAN -- which
// VCF key each came from, its declared type, and its Number.
func (s *ParquetVolume) FormatFields() []FormatField {
	if s.manifest == nil {
		return nil
	}
	return s.manifest.Params.Format
}

// CallsFormat walks the ALT calls yielding each with its captured FORMAT
// values.
//
// The FormatRow is only valid for the duration of the call -- it views the
// scan's reusable buffer -- so copy what you keep. A store that captured
// nothing walks normally with an empty row, so a caller need not branch.
func (s *ParquetVolume) CallsFormat(fn func(Call, FormatRow) bool) error {
	if len(s.FormatFields()) == 0 {
		return scanParquet(s.calls, func(c Call) bool { return fn(c, FormatRow{}) })
	}
	for _, sh := range s.calls.shards {
		pf, err := sh.m.parsed()
		if err != nil {
			return err
		}
		schema := pf.Schema()
		cols := map[string]int{}
		for _, path := range schema.Columns() {
			if lc, ok := schema.Lookup(path...); ok {
				cols[strings.Join(path, ".")] = lc.ColumnIndex
			}
		}
		r := parquet.NewGenericReader[any](pf, schema)
		buf := make([]parquet.Row, 512)
		stop := false
		for !stop {
			n, err := r.ReadRows(buf)
			for i := 0; i < n; i++ {
				row := buf[i]
				if !fn(callFromRow(row, cols), FormatRow{row: row, cols: cols}) {
					stop = true
					break
				}
			}
			if err == io.EOF || n == 0 {
				break
			}
			if err != nil {
				r.Close()
				return err
			}
		}
		r.Close()
		if stop {
			return nil
		}
	}
	return nil
}

// callFromRow rebuilds a Call from a dynamically-read row, BY NAME rather than
// by position: a group-built schema orders its columns alphabetically rather
// than in struct order, so an index taken from one would read the wrong column
// in the other.
func callFromRow(row parquet.Row, cols map[string]int) Call {
	str := func(name string) string {
		if i, ok := cols[name]; ok && i < len(row) && !row[i].IsNull() {
			return row[i].Clone().String()
		}
		return ""
	}
	i32 := func(name string) int32 {
		if i, ok := cols[name]; ok && i < len(row) && !row[i].IsNull() {
			return row[i].Int32()
		}
		return 0
	}
	return Call{
		SampleID: str("sample_id"), Chrom: str("chrom"), Pos: i32("pos"),
		Ref: str("ref"), Alt: str("alt"), RefEnd: i32("ref_end"),
		GT: str("gt"), DP: i32("dp"), ADRef: i32("ad_ref"),
		ADAlt: i32("ad_alt"), GQ: i32("gq"), MinDP: i32("min_dp"),
	}
}

// InfoFields reports the INFO fields this store captured, nil if none.
//
// Read from the manifest rather than from the file's schema, because the two
// answer different questions: the schema says which columns exist, the manifest
// says what they MEAN -- which VCF key each came from, its declared type, and
// its Number. A column called info_r2 alone cannot tell a consumer whether it
// holds minimac's R2 or something a caller mapped onto that name.
func (s *ParquetVolume) InfoFields() []InfoField {
	if s.manifest == nil {
		return nil
	}
	return s.manifest.Params.Info
}

// SitesInfo walks the catalog yielding each site with its captured INFO values.
//
// The InfoRow is only valid for the duration of the call: it is a view onto the
// scan's reusable buffer, so a caller keeping one past its callback gets
// whatever site came next. Copy what you need.
//
// A store that captured nothing walks normally with an empty InfoRow, so a
// caller need not branch on whether capture was on -- Present reports false for
// every field, which is the truth.
func (s *ParquetVolume) SitesInfo(fn func(Site, InfoRow) bool) error {
	if !s.hasSites {
		return fmt.Errorf("%s is missing", SitesPath(s.base))
	}
	if len(s.InfoFields()) == 0 {
		return s.Sites(func(site Site) bool { return fn(site, InfoRow{}) })
	}

	// The captured-INFO columns are the same in every shard, so the first one
	// describes the table.
	pf, err := s.sites.single().parsed()
	if err != nil {
		return err
	}
	schema := pf.Schema()

	// Column indexes by name, resolved once. The dynamic schema is ordered by
	// name rather than by struct order, so nothing may assume a position.
	cols := map[string]int{}
	for _, path := range schema.Columns() {
		lc, ok := schema.Lookup(path...)
		if !ok {
			continue
		}
		cols[strings.Join(path, ".")] = lc.ColumnIndex
	}

	// Sites come back through the typed reader as well as the row reader: the
	// struct decoding is the part worth not reimplementing, and reading the
	// file twice would double the IO. So rows are read once and Site is
	// reconstructed from the same row.
	r := parquet.NewGenericReader[any](pf, schema)
	defer r.Close()

	buf := make([]parquet.Row, 512)
	for {
		n, err := r.ReadRows(buf)
		for i := 0; i < n; i++ {
			row := buf[i]
			site := siteFromRow(row, cols)
			if !fn(site, InfoRow{row: row, cols: cols}) {
				return nil
			}
		}
		if err == io.EOF {
			return nil
		}
		if err != nil {
			return err
		}
		if n == 0 {
			return nil
		}
	}
}

// siteFromRow rebuilds a Site from a dynamically-read row.
//
// By name, never by position: a group-built schema orders its columns
// alphabetically rather than in struct order, so the two do not agree and an
// index taken from one would silently read the wrong column in the other.
func siteFromRow(row parquet.Row, cols map[string]int) Site {
	str := func(name string) string {
		if i, ok := cols[name]; ok && i < len(row) && !row[i].IsNull() {
			return row[i].Clone().String()
		}
		return ""
	}
	i32 := func(name string) int32 {
		if i, ok := cols[name]; ok && i < len(row) && !row[i].IsNull() {
			return row[i].Int32()
		}
		return 0
	}
	return Site{
		Chrom:     str("chrom"),
		Pos:       i32("pos"),
		Ref:       str("ref"),
		Alt:       str("alt"),
		RefEnd:    i32("ref_end"),
		AC:        i32("ac"),
		AN:        i32("an"),
		NCarriers: i32("n_carriers"),
		NCalled:   i32("n_called"),
		NLowDP:    i32("n_lowdp"),
	}
}

// Site returns the catalog entry for a locus, if the source reported it.
func (s *ParquetVolume) Site(l Locus) (Site, bool, error) {
	// Every sibling guards this; Site did not, and an absent sites table is a
	// reachable state -- the manifest permits one that recorded no rows -- so
	// this was a nil dereference inside the scan rather than the
	// ErrNotClassifiable its neighbours return.
	if !s.hasSites {
		return Site{}, false, fmt.Errorf("%w: %s is missing", ErrNotClassifiable, SitesPath(s.base))
	}
	var got Site
	found := false
	err := scanParquetPruned(s.sites, locusFilter(l), func(site Site) bool {
		if SameLocus(site.Locus(), l) {
			got, found = site, true
			return false
		}
		return true
	})
	return got, found, err
}

// SiteKnown reports whether a locus appears in the sites catalog, i.e. whether
// the source actually reported it. For a store built from a plain VCF this is
// the boundary of what can be answered at all, so callers presenting results to
// a user should say so rather than let "0 carriers" read as "nobody carries it".
func (s *ParquetVolume) SiteKnown(l Locus) (bool, error) {
	// Site already does this scan, with the same guard, the same filter and the
	// same early stop; this was a second copy that discarded the payload.
	_, ok, err := s.Site(l)
	return ok, err
}

// Close releases the three tables, which are held open for the store's life.
//
// It is not optional. This comment used to say the opposite -- "a no-op;
// ParquetVolume opens files per query" -- which was true before the tables
// became long-lived, and a caller trusting it leaks a file handle per store,
// or a live connection per store when the tables are remote.
func (s *ParquetVolume) Close() error {
	var first error
	for _, m := range []*shardSet{s.calls, s.sites, s.regions} {
		if err := m.Close(); err != nil && first == nil {
			first = err
		}
	}
	return first
}

// HasRefSpans reports whether the store records each record's reference span.
//
// Without it, a span query matches on start position alone, so a deletion that
// begins before the region and reaches into it is not returned. That is a
// difference in the answer rather than in the speed, so a caller offering region
// queries should say so rather than let the omission look like a real negative.
// Stores written before the ref_end column existed report false; rewriting the
// store is what fixes it.
func (s *ParquetVolume) HasRefSpans() bool { return s.hasRefSpans }

// VolumeManifest returns the completion record this store was opened with.
//
// It is never nil for an open store: opening requires one. Callers wanting the
// per-chromosome census, the conversion parameters or the recorded counts read
// them from here rather than re-deriving them, which is the point of writing
// them down.
func (s *ParquetVolume) VolumeManifest() *VolumeManifest { return s.manifest }

// setCallMeta writes key/value metadata to whichever calls writer is live.
func (w *Writer) setCallMeta(k, v string) {
	if w.callsAny != nil {
		w.callsAny.SetKeyValueMetadata(k, v)
		return
	}
	w.calls.SetKeyValueMetadata(k, v)
}
