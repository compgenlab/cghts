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
// key and a JSON object member, and be greppable afterwards, so it is held to
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
	// Sink is where the members are written. Nil means "work it out from the
	// base", which is what every caller wants; it is settable so a test can
	// supply one without a filesystem or a bucket.
	Sink Sink

	Codec        compress.Codec
	RowGroupSize int64
	Samples      []string
	MinDP        int32
	NoCallable   bool
	Program      string
	Command      string

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
	calls   *parquet.GenericWriter[Call]
	sites   *parquet.GenericWriter[Site]
	regions *parquet.GenericWriter[CalledSiteRun]

	callBuf   []Call
	siteBuf   []Site
	regionBuf []CalledSiteRun

	// The members being written, in creation order, with the names they will
	// have. Held so Close can finish them and abort can undo them -- and named
	// rather than handled, because a remote member has no file to name itself.
	sink    Sink
	members []memberSink

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

	NCalls   int64
	NSites   int64
	NRegions int64
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

// memberSink is one store member being written: the writer it is fed through
// and the file name it will have. Named around the sink because `member` and
// `openMember` both already mean something on the reading side.
type memberSink struct {
	name string
	sw   io.WriteCloser
}

// NewWriter creates the store directory at base and the member files inside it.
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
	// Remove any manifest before touching the members. From here until Finish
	// the store is under construction and must not carry a completion marker --
	// otherwise a --force re-run that dies partway would leave the previous
	// run's manifest vouching for this run's half-written members.
	if err := sink.Remove(ManifestFile); err != nil {
		return nil, fmt.Errorf("clearing previous manifest: %w", err)
	}
	w := &Writer{base: base, opts: o, sink: sink}

	open := func(name string) (io.Writer, error) {
		f, err := sink.Create(name)
		if err != nil {
			_ = w.abort()
			return nil, err
		}
		w.members = append(w.members, memberSink{name: name, sw: f})
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

	cf, err := open(MemberFile(CallsMember))
	if err != nil {
		return nil, err
	}
	w.calls = parquet.NewGenericWriter[Call](cf, callOpts...)
	w.calls.SetKeyValueMetadata(MetaSamples, strings.Join(o.Samples, "\n"))
	w.calls.SetKeyValueMetadata(MetaMinDP, fmt.Sprint(o.MinDP))
	w.calls.SetKeyValueMetadata(MetaProgram, o.Program)
	w.calls.SetKeyValueMetadata(MetaCommand, o.Command)
	w.calls.SetKeyValueMetadata(MetaSource, strings.Join(o.Sources, ", "))
	if len(o.Contigs) > 0 {
		w.calls.SetKeyValueMetadata(MetaContigs, strings.Join(o.Contigs, "\n"))
	}
	w.calls.SetKeyValueMetadata(MetaSpans, string(o.Spans))
	for _, k := range sortedKeys(o.Meta) {
		w.calls.SetKeyValueMetadata(MetaPrefix+k, o.Meta[k])
	}
	if o.NoCallable {
		w.calls.SetKeyValueMetadata(MetaNoCall, "1")
	}

	sf, err := open(MemberFile(SitesMember))
	if err != nil {
		return nil, err
	}
	w.sites = parquet.NewGenericWriter[Site](sf, sortedByLocus...)

	rf, err := open(MemberFile(RegionsMember))
	if err != nil {
		return nil, err
	}
	w.regions = parquet.NewGenericWriter[CalledSiteRun](rf, append(append([]parquet.WriterOption{}, opts...),
		parquet.SortingWriterConfig(parquet.SortingColumns(
			parquet.Ascending("chrom"), parquet.Ascending("start"))),
		parquet.BloomFilters(parquet.SplitBlockFilter(10, "sample_id")))...)

	return w, nil
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
	// survive, and a completion marker sitting beside discarded members is the
	// worst outcome available.
	if err := w.sink.Remove(ManifestFile); err != nil {
		errs = append(errs, fmt.Errorf("removing %s: %w", ManifestFile, err))
	}
	aborter, remote := w.sink.(Aborter)
	for _, m := range w.members {
		// ABANDON where a member only becomes visible when it is finished, and
		// UNLINK where it is visible as it is written. On an object store there
		// is no half-written member to delete -- what exists is an upload in
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
	w.members = nil
	return errors.Join(errs...)
}

// Discard abandons a conversion, leaving nothing behind.
//
// A failure partway through must not leave the members on disk. They would look
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

// Finish closes the members and then writes the manifest that marks the store
// complete and readable.
//
// This is the ordinary way to end a conversion; Close alone leaves a store
// without a manifest, which readers refuse. The order is the whole point: every
// member is finalized first, and only a run that got that far writes the marker
// saying so.
// A failure after Close discards the store rather than leaving it. By then the
// three members carry valid footers, so what is on disk looks exactly like a
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

// manifest describes the store as written. Member sizes are read back from disk
// rather than accumulated, so they describe the files rather than the intent.
func (w *Writer) manifest() (Manifest, error) {
	members := map[string]MemberInfo{}
	for member, rows := range map[string]int64{
		CallsMember:   w.NCalls,
		SitesMember:   w.NSites,
		RegionsMember: w.NRegions,
	} {
		size, ok, err := w.sink.Stat(MemberFile(member))
		if err != nil {
			return Manifest{}, fmt.Errorf("sizing %s member: %w", member, err)
		}
		if !ok {
			return Manifest{}, fmt.Errorf("the %s member is missing from %s", member, w.sink.Describe())
		}
		members[member] = MemberInfo{Rows: rows, Bytes: size}
	}
	return Manifest{
		FormatVersion: ManifestVersion,
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
		},
		Samples: w.opts.Samples,
		Counts: ManifestCounts{
			Samples: len(w.opts.Samples),
			Calls:   w.NCalls,
			Sites:   w.NSites,
			Regions: w.NRegions,
		},
		Members:         members,
		Chromosomes:     w.chroms,
		ContigsDeclared: w.opts.Contigs,
	}, nil
}

// WriteCall buffers one ALT-carrying genotype.
func (w *Writer) WriteCall(c Call) error {
	c.RefEnd = defaultRefEnd(c.Pos, c.Ref, c.RefEnd)
	w.census(c.Chrom, c.Pos).Calls++
	w.callBuf = append(w.callBuf, c)
	w.NCalls++
	if len(w.callBuf) >= batchSize {
		return w.flushCalls()
	}
	return nil
}

// WriteSite buffers one catalog entry.
func (w *Writer) WriteSite(s Site) error {
	s.RefEnd = defaultRefEnd(s.Pos, s.Ref, s.RefEnd)
	w.census(s.Chrom, s.Pos).Sites++
	w.siteBuf = append(w.siteBuf, s)
	w.NSites++
	if len(w.siteBuf) >= batchSize {
		return w.flushSites()
	}
	return nil
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
	w.regionBuf = append(w.regionBuf, r)
	w.NRegions++
	if len(w.regionBuf) >= batchSize {
		return w.flushRegions()
	}
	return nil
}

func (w *Writer) flushCalls() error {
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

// Close flushes and finalizes all three files.
//
// It stops at the first failure instead of pressing on. Continuing would finish
// the other members, and a parquet footer is written precisely by the Close that
// would then run -- so a failed flush used to yield three structurally valid
// files of which one was silently short by up to a batch. Three well-formed
// members that disagree is the one outcome a reader cannot detect, so the write
// stops as soon as it cannot be completed honestly.
//
// A non-nil error means nothing on disk should be trusted; the caller is
// expected to Discard. The file handles are released either way, since leaking
// them helps nobody, but the members are left in place for Discard to remove.
func (w *Writer) Close() error {
	for _, fn := range []func() error{w.flushCalls, w.flushSites, w.flushRegions} {
		if err := fn(); err != nil {
			w.closeFiles()
			return err
		}
	}
	for _, c := range []io.Closer{w.calls, w.sites, w.regions} {
		if err := c.Close(); err != nil {
			w.closeFiles()
			return err
		}
	}
	return w.closeFiles()
}

// closeFiles releases the underlying handles without removing anything.
func (w *Writer) closeFiles() error {
	var errs []error
	for _, m := range w.members {
		if err := m.sw.Close(); err != nil {
			errs = append(errs, fmt.Errorf("finishing %s: %w", m.name, err))
		}
	}
	return errors.Join(errs...)
}

// scanParquet streams path in batches, calling fn per row. fn returns false to
// stop early. Batching keeps a whole-genome store from having to be resident.
func scanParquet[T any](m *member, fn func(T) bool) error {
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

// ParquetStore is a Store backed by the three-file Parquet set.
type ParquetStore struct {
	base        string
	calls       *member
	hasRefSpans bool
	sites       *member
	regions     *member
	samples     []string
	hasSites    bool
	hasRegions  bool
	noCallable  bool
	spans       SpanSemantics
	minDP       int32
	meta        map[string]string
	manifest    *Manifest
}

// SpanSemantics reports what this store's run intervals may claim. A store
// written from a plain VCF reports SpansSites, confining answers to the sites
// catalog.
func (s *ParquetStore) SpanSemantics() SpanSemantics { return s.spans }

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
func (s *ParquetStore) Provenance() Provenance {
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
func (s *ParquetStore) metaFields() map[string]string {
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
func OpenParquet(base string) (*ParquetStore, error) {
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
// The members stay open for the life of the store, so Close matters more than
// it used to: it now releases them rather than being a no-op.
func OpenParquetContext(ctx context.Context, base string) (*ParquetStore, error) {
	base = TrimStoreSuffix(base)
	man, err := requireManifest(ctx, base)
	if err != nil {
		return nil, err
	}
	calls, err := openMember(ctx, CallsPath(base))
	if err != nil {
		return nil, fmt.Errorf("opening parquet store %s: %w", base, err)
	}
	pf, err := parquet.OpenFile(calls.ra, calls.size)
	if err != nil {
		calls.Close()
		return nil, fmt.Errorf("reading %s: %w", calls.name, err)
	}
	s := &ParquetStore{base: base, calls: calls, meta: map[string]string{}}
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
	// The optional members. Absence is legitimate only where the manifest
	// recorded nothing in them, which verifyAgainstManifest below enforces;
	// note that the writer creates all three regardless, so --no-callable
	// produces a present, zero-row regions member rather than a missing one.
	sites, hasSites, err := openOptionalMember(ctx, SitesPath(base))
	if err != nil {
		calls.Close()
		return nil, err
	}
	s.sites, s.hasSites = sites, hasSites
	regions, hasRegions, err := openOptionalMember(ctx, RegionsPath(base))
	if err != nil {
		calls.Close()
		sites.Close()
		return nil, err
	}
	s.regions, s.hasRegions = regions, hasRegions
	s.manifest = man

	if err := verifyAgainstManifest(man, pf, sites, regions); err != nil {
		s.Close()
		return nil, fmt.Errorf("%s: %w", base, err)
	}
	return s, nil
}

// requireManifest loads the completion marker, refusing the store without it.
//
// A store with no manifest was not finished -- or predates them -- and either
// way nothing here can tell how much of the intended input it holds. The
// members would open and answer queries perfectly well; that is the problem.
// An unfinished store reports "not assayed" for everything it never got to,
// which is indistinguishable from the honest answer for a position the source
// genuinely never reported.
func requireManifest(ctx context.Context, base string) (*Manifest, error) {
	man, err := ReadManifestContext(ctx, base)
	if err != nil {
		return nil, fmt.Errorf("%s has no readable %s, so it cannot be shown to be "+
			"a complete store: %w\n\tit is from an interrupted conversion, or predates "+
			"manifests; re-convert it, or inspect it with vcf-varsummary",
			base, ManifestFile, err)
	}
	if !man.Complete {
		return nil, fmt.Errorf("%s is marked incomplete by its %s; re-convert it",
			base, ManifestFile)
	}
	if man.FormatVersion > ManifestVersion {
		return nil, fmt.Errorf("%s was written by a newer cgkit (manifest format %d, "+
			"this build understands %d); upgrade to read it",
			base, man.FormatVersion, ManifestVersion)
	}
	return man, nil
}

// verifyAgainstManifest checks each member against the row count recorded for
// it when the store was written.
//
// The footers say each member was finished; the manifest says which members
// finishing this store produced. Comparing them is what catches a member that
// is well-formed but does not belong -- one copied in from another conversion,
// or left over from a run whose replacement died. Nothing else can see that,
// because sites and regions carry no metadata of their own.
//
// It costs nothing: every count here is footer metadata that was already read.
func verifyAgainstManifest(man *Manifest, calls *parquet.File, sites, regions *member) error {
	type check struct {
		name string
		file *parquet.File
		m    *member
	}
	for _, c := range []check{
		{CallsMember, calls, nil},
		{SitesMember, nil, sites},
		{RegionsMember, nil, regions},
	} {
		want, recorded := man.Members[c.name]
		if !recorded {
			continue // a manifest that never described this member claims nothing
		}
		if c.file == nil {
			if c.m == nil {
				if want.Rows == 0 {
					continue // written empty and since removed; it claimed nothing
				}
				return fmt.Errorf("%s.parquet is missing, but the manifest records %d rows in it",
					c.name, want.Rows)
			}
			f, err := parquet.OpenFile(c.m.ra, c.m.size)
			if err != nil {
				return fmt.Errorf("reading %s.parquet: %w", c.name, err)
			}
			c.file = f
		}
		if got := c.file.NumRows(); got != want.Rows {
			return fmt.Errorf("%s.parquet holds %d rows but the manifest records %d; "+
				"this member does not belong to this store, or the store is damaged",
				c.name, got, want.Rows)
		}
	}
	return nil
}

// openOptionalMember opens a member a store may legitimately lack, and verifies
// that what it finds is a readable parquet file.
//
// The two outcomes have to stay distinct. Absence is normal: a --no-callable
// store records no callable runs. A member that is *present but unreadable* is
// not, and until this checked, the two were indistinguishable -- the open error
// was discarded either way, so a truncated or half-written sites file read
// exactly like a deliberate omission and the store went on to answer queries as
// though those runs had never been asked for.
//
// Parsing the footer here is what makes the check meaningful, and it is nearly
// free: it is footer metadata, not a data scan, and every query would have
// parsed it anyway. A parquet footer is written only by the writer's Close, so
// a member that parses is a member that was finished.
//
// One limit worth naming: for a remote locator "absent" and "unreachable" are
// not cleanly separable here, because a 404 surfaces as a failed size probe
// rather than as fs.ErrNotExist. Such a member is treated as absent, which is
// the pre-existing behaviour. The manifest is what makes this exact, since it
// records which members the conversion actually wrote.
func openOptionalMember(ctx context.Context, locator string) (*member, bool, error) {
	m, err := openMember(ctx, locator)
	if err != nil {
		// Absent is an answer; anything else is a failure. This used to discard
		// every error alike, so a member that existed but could not be read --
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
//	cohort/calls.parquet      a member of it
//	cohort/manifest.json.gz   likewise
//
// Pointing at a member is worth accepting because tab completion lands there,
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
		CallsMember + ".parquet", SitesMember + ".parquet", RegionsMember + ".parquet",
		ManifestFile,
	} {
		if base, ok := cutMemberSuffix(p, name); ok {
			return base
		}
	}
	return trimStoreDir(p)
}

// cutMemberSuffix removes a trailing "/name" (or the platform separator) from p,
// reporting whether it was there. A bare "name" with no separator is a member of
// the current directory, whose store base is "".
func cutMemberSuffix(p, name string) (string, bool) {
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
func (s *ParquetStore) Contigs() []string {
	v := s.meta[MetaContigs]
	if v == "" {
		return nil
	}
	return strings.Split(v, "\n")
}

// Samples returns the roster recorded at conversion time.
func (s *ParquetStore) Samples() ([]string, error) {
	if len(s.samples) == 0 {
		return nil, fmt.Errorf("%s records no sample list", CallsPath(s.base))
	}
	return s.samples, nil
}

// Calls streams the genotypes a query selects, in the store's own order.
func (s *ParquetStore) Calls(q Query) (iter.Seq2[Call, error], error) {
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
func (s *ParquetStore) callsWithRef(p *plan, keep rowGroupFilter) ([]Call, error) {
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
func (s *ParquetStore) Classify(l Locus, g Gate) ([]SampleState, error) {
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
func (s *ParquetStore) classifiable() error {
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
func (s *ParquetStore) runsFor(samples []string, p *plan) (runsBySample, error) {
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
func (s *ParquetStore) Sites(fn func(Site) bool) error {
	if !s.hasSites {
		return fmt.Errorf("%s is missing", SitesPath(s.base))
	}
	return scanParquet(s.sites, fn)
}

// Site returns the catalog entry for a locus, if the source reported it.
func (s *ParquetStore) Site(l Locus) (Site, bool, error) {
	// Every sibling guards this; Site did not, and an absent sites member is a
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
func (s *ParquetStore) SiteKnown(l Locus) (bool, error) {
	// Site already does this scan, with the same guard, the same filter and the
	// same early stop; this was a second copy that discarded the payload.
	_, ok, err := s.Site(l)
	return ok, err
}

// Close releases the three members, which are held open for the store's life.
//
// It is not optional. This comment used to say the opposite -- "a no-op;
// ParquetStore opens files per query" -- which was true before the members
// became long-lived, and a caller trusting it leaks a file handle per store,
// or a live connection per store when the members are remote.
func (s *ParquetStore) Close() error {
	var first error
	for _, m := range []*member{s.calls, s.sites, s.regions} {
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
func (s *ParquetStore) HasRefSpans() bool { return s.hasRefSpans }

// Manifest returns the completion record this store was opened with.
//
// It is never nil for an open store: opening requires one. Callers wanting the
// per-chromosome census, the conversion parameters or the recorded counts read
// them from here rather than re-deriving them, which is the point of writing
// them down.
func (s *ParquetStore) Manifest() *Manifest { return s.manifest }
