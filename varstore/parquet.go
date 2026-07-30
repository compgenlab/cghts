package varstore

import (
	"errors"
	"fmt"
	"io"
	"iter"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"

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
)

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
	Codec        compress.Codec
	RowGroupSize int64
	Samples      []string
	MinDP        int32
	NoCallable   bool
	Program      string
	Command      string
	Source       string

	// Spans declares what the run intervals may claim. Defaults to SpansSites,
	// which is all a plain VCF can support.
	Spans SpanSemantics
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

	files []*os.File

	NCalls   int64
	NSites   int64
	NRegions int64
}

const batchSize = 8192

// NewWriter creates the three store files for base.
func NewWriter(base string, o WriterOpts) (*Writer, error) {
	if o.Codec == nil {
		o.Codec = &zstd.Codec{}
	}
	if o.RowGroupSize <= 0 {
		o.RowGroupSize = 250_000
	}
	w := &Writer{}

	open := func(path string) (*os.File, error) {
		f, err := os.Create(path)
		if err != nil {
			w.abort()
			return nil, err
		}
		w.files = append(w.files, f)
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

	cf, err := open(CallsPath(base))
	if err != nil {
		return nil, err
	}
	w.calls = parquet.NewGenericWriter[Call](cf, callOpts...)
	w.calls.SetKeyValueMetadata(MetaSamples, strings.Join(o.Samples, "\n"))
	w.calls.SetKeyValueMetadata(MetaMinDP, fmt.Sprint(o.MinDP))
	w.calls.SetKeyValueMetadata(MetaProgram, o.Program)
	w.calls.SetKeyValueMetadata(MetaCommand, o.Command)
	w.calls.SetKeyValueMetadata(MetaSource, o.Source)
	if o.Spans == "" {
		o.Spans = SpansSites
	}
	w.calls.SetKeyValueMetadata(MetaSpans, string(o.Spans))
	if o.NoCallable {
		w.calls.SetKeyValueMetadata(MetaNoCall, "1")
	}

	sf, err := open(SitesPath(base))
	if err != nil {
		return nil, err
	}
	w.sites = parquet.NewGenericWriter[Site](sf, sortedByLocus...)

	rf, err := open(RegionsPath(base))
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
func (w *Writer) abort() {
	for _, f := range w.files {
		f.Close()
		os.Remove(f.Name())
	}
	w.files = nil
}

// Discard abandons a conversion, leaving nothing behind.
//
// A failure partway through must not leave the members on disk. They would look
// like a store, they would be truncated or incomplete, and -- because
// conversion refuses to overwrite an existing store -- their presence would
// block the retry.
func (w *Writer) Discard() {
	for _, c := range []io.Closer{w.calls, w.sites, w.regions} {
		if c != nil {
			_ = c.Close()
		}
	}
	w.abort()
}

// WriteCall buffers one ALT-carrying genotype.
func (w *Writer) WriteCall(c Call) error {
	w.callBuf = append(w.callBuf, c)
	w.NCalls++
	if len(w.callBuf) >= batchSize {
		return w.flushCalls()
	}
	return nil
}

// WriteSite buffers one catalog entry.
func (w *Writer) WriteSite(s Site) error {
	w.siteBuf = append(w.siteBuf, s)
	w.NSites++
	if len(w.siteBuf) >= batchSize {
		return w.flushSites()
	}
	return nil
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
func (w *Writer) Close() error {
	var errs []error
	for _, fn := range []func() error{w.flushCalls, w.flushSites, w.flushRegions} {
		if err := fn(); err != nil {
			errs = append(errs, err)
		}
	}
	for _, c := range []io.Closer{w.calls, w.sites, w.regions} {
		if err := c.Close(); err != nil {
			errs = append(errs, err)
		}
	}
	for _, f := range w.files {
		if err := f.Close(); err != nil {
			errs = append(errs, err)
		}
	}
	return errors.Join(errs...)
}

// scanParquet streams path in batches, calling fn per row. fn returns false to
// stop early. Batching keeps a whole-genome store from having to be resident.
func scanParquet[T any](path string, fn func(T) bool) error {
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()
	r := parquet.NewGenericReader[T](f)
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

// ParquetStore is a Store backed by the three-file Parquet set.
type ParquetStore struct {
	base       string
	samples    []string
	hasSites   bool
	hasRegions bool
	noCallable bool
	spans      SpanSemantics
	minDP      int32
	meta       map[string]string
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
	}
}

// OpenParquet opens a Parquet store. base may be given either as the base name
// or as the path to any one of the three files.
func OpenParquet(base string) (*ParquetStore, error) {
	base = TrimStoreSuffix(base)
	callsPath := CallsPath(base)
	f, err := os.Open(callsPath)
	if err != nil {
		return nil, fmt.Errorf("opening parquet store %s: %w", base, err)
	}
	defer f.Close()
	st, err := f.Stat()
	if err != nil {
		return nil, err
	}
	pf, err := parquet.OpenFile(f, st.Size())
	if err != nil {
		return nil, fmt.Errorf("reading %s: %w", callsPath, err)
	}
	s := &ParquetStore{base: base, meta: map[string]string{}}
	for _, k := range []string{MetaSource, MetaProgram, MetaCommand, MetaMinDP} {
		if v, ok := pf.Lookup(k); ok {
			s.meta[k] = v
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
	s.hasSites = fileExists(SitesPath(base))
	s.hasRegions = fileExists(RegionsPath(base))
	return s, nil
}

// TrimStoreSuffix reduces any member path of a store to its base name.
// It accepts every spelling a user might reasonably type:
//
//	cohort                  prefix form
//	cohort.calls.parquet    a member of the prefix form
//	cohort/                 directory form
//	cohort/calls.parquet    a member of the directory form
//	cohort                  a directory, written without the trailing slash
//
// The last case is why this consults the filesystem: having asked for
// "--out cohort/", nobody should have to remember the slash to read it back.
func TrimStoreSuffix(p string) string {
	// A member inside a directory-form store.
	for _, m := range []string{CallsMember, SitesMember, RegionsMember} {
		if filepath.Base(p) == m+".parquet" {
			return ensureTrailingSep(filepath.Dir(p))
		}
	}
	// A member of a prefix-form store.
	for _, sfx := range []string{CallsSuffix, SitesSuffix, RegionsSuffix} {
		if strings.HasSuffix(p, sfx) {
			return strings.TrimSuffix(p, sfx)
		}
	}
	if IsDirBase(p) {
		return p
	}
	// A bare name that is really a directory holding a store.
	if st, err := os.Stat(p); err == nil && st.IsDir() {
		if fileExists(filepath.Join(p, CallsMember+".parquet")) {
			return ensureTrailingSep(p)
		}
	}
	return p
}

// ensureTrailingSep marks a path as a directory base.
func ensureTrailingSep(p string) string {
	if p == "" {
		return "." + string(os.PathSeparator)
	}
	if IsDirBase(p) {
		return p
	}
	return p + string(os.PathSeparator)
}

func fileExists(p string) bool {
	st, err := os.Stat(p)
	return err == nil && !st.IsDir()
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
			err := scanParquetPruned(CallsPath(s.base), keep, func(c Call) bool {
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
	if err := scanParquetPruned(CallsPath(s.base), keep, func(c Call) bool {
		if !p.wantsSample(c.SampleID) {
			return true
		}
		loc := c.Locus()
		// Record level, not site level: a call at a SIBLING locus of the same
		// record still disqualifies the sample from being reference here, and
		// testing wantsSite would throw that evidence away.
		if !p.wantsRecord(loc) {
			return true
		}
		if carries[c.SampleID] == nil {
			carries[c.SampleID] = map[RecordKey]bool{}
		}
		carries[c.SampleID][loc.Record()] = true
		if p.wantsSite(loc) && p.q.Gate.Admits(c) {
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
	err = scanParquetPruned(SitesPath(s.base), siteScanFilter(p.q), func(site Site) bool {
		loc := site.Locus()
		if !p.wantsSite(loc) {
			return true
		}
		k := canonLocus(loc)
		// ALT rows first, already in roster order from the scan above.
		out = append(out, altBySite[k]...)
		delete(altBySite, k)
		for _, name := range wanted {
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
	// dropping it silently would be worse than emitting it out of order.
	for _, calls := range altBySite {
		out = append(out, calls...)
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
	if err := scanParquetPruned(CallsPath(s.base), locusFilter(l), func(c Call) bool {
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
	if err := scanParquetPruned(RegionsPath(s.base), coveringFilter(l.Chrom, l.Pos), func(r CalledSiteRun) bool {
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
	err := scanParquetPruned(RegionsPath(s.base), runScanFilter(p.q), func(r CalledSiteRun) bool {
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
	return scanParquet(SitesPath(s.base), fn)
}

// Site returns the catalog entry for a locus, if the source reported it.
func (s *ParquetStore) Site(l Locus) (Site, bool, error) {
	var got Site
	found := false
	err := scanParquetPruned(SitesPath(s.base), locusFilter(l), func(site Site) bool {
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
	if !s.hasSites {
		return false, fmt.Errorf("%w: %s is missing", ErrNotClassifiable, SitesPath(s.base))
	}
	found := false
	err := scanParquetPruned(SitesPath(s.base), locusFilter(l), func(site Site) bool {
		if SameLocus(site.Locus(), l) {
			found = true
			return false
		}
		return true
	})
	return found, err
}

// Close is a no-op; ParquetStore opens files per query.
func (s *ParquetStore) Close() error { return nil }
