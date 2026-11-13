package varstore

import (
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"time"

	"github.com/compgenlab/cghts/iosource"
	"github.com/parquet-go/parquet-go"
	"github.com/parquet-go/parquet-go/compress/zstd"
)

// A ref_blocks FRAGMENT: one sample's reference blocks, from one gVCF, waiting to
// be attached to a store.
//
// WHY A FRAGMENT RATHER THAN WRITING INTO THE STORE. A cohort's gVCFs are read
// by thousands of independent jobs, one per sample, and a store's ref_blocks is
// sorted across all of them. A job can only write what it alone knows; the merge
// into the store's order is a later, separate step (attach). So a fragment is
// the unit of work: single-sample, written in one streaming pass, and safe to
// throw away and redo.
//
// THE LAYOUT. A directory holding one Parquet file per chromosome, named by
// index ("00000.parquet") because contig names may hold characters a path should
// not, and fragment.json, which names them:
//
//	S0001.frag/
//	  00000.parquet   chr1 blocks
//	  00001.parquet   chr2 blocks
//	  fragment.json   written last; marks the fragment complete
//
// THE SIDECAR IS THE COMPLETION MARKER, exactly as manifest.json.gz is for a
// store: removed first when a fragment is (re)written, written last when it is
// finished. A job killed partway leaves no marker, so attach refuses it rather
// than merging half a genome as though it were the whole -- and a scheduler can
// rerun any job whose marker is missing without asking what state it reached.

// FragmentMetaFile is the completion marker inside a fragment directory.
const FragmentMetaFile = "fragment.json"

// FragmentFormatVersion is the version of the fragment layout this package writes.
const FragmentFormatVersion = 1

// FragmentMeta describes a fragment and the settings it was built with.
//
// Attach compares these across fragments and against the store: blocks banded
// or gap-merged differently do not mean the same thing by "covered", and a table
// mixing them would answer one gate differently for different people.
type FragmentMeta struct {
	FormatVersion int       `json:"format_version"`
	Complete      bool      `json:"complete"`
	Created       time.Time `json:"created"`
	Program       string    `json:"program,omitempty"`
	Command       string    `json:"command,omitempty"`
	Source        string    `json:"source"`

	SampleID    string  `json:"sample_id"`
	MinDP       int32   `json:"min_dp"`
	DepthBands  []int32 `json:"depth_bands,omitempty"`
	DepthRatio  float64 `json:"depth_ratio,omitempty"`
	MaxGap      int32   `json:"max_gap"`
	Passthrough bool    `json:"passthrough,omitempty"`

	// Chroms are the per-chromosome files, in the order written -- the source's
	// own contig order.
	Chroms []FragmentChrom `json:"chroms"`
}

// FragmentChrom is one chromosome's file in a fragment.
type FragmentChrom struct {
	Name  string `json:"name"`
	File  string `json:"file"`
	Rows  int64  `json:"rows"`
	First int32  `json:"first"` // start of the first block
	Last  int32  `json:"last"`  // end of the last block
}

// FragmentWriter writes one fragment.
type FragmentWriter struct {
	sink  Sink
	meta  FragmentMeta
	files []string // every file created, for Discard

	out   io.WriteCloser
	pw    *parquet.GenericWriter[RefBlock]
	buf   []RefBlock
	cur   *FragmentChrom
	seen  map[string]bool
	ended int32 // end of the last block written on the current chromosome
	done  bool
}

// NewFragmentWriter starts a fragment at dir, local or remote. meta's SampleID
// and settings are recorded; Chroms, Complete and FormatVersion are the
// writer's to fill.
//
// A previous fragment.json at dir is removed before anything else is written,
// for the same reason NewWriter clears a store's manifest first: otherwise a
// rerun that dies would leave the old marker vouching for the new run's files.
func NewFragmentWriter(dir string, meta FragmentMeta) (*FragmentWriter, error) {
	if meta.SampleID == "" {
		return nil, fmt.Errorf("a fragment needs the sample it describes")
	}
	if err := ValidateDepthBands(meta.DepthBands); err != nil {
		return nil, err
	}
	if meta.MaxGap < 0 {
		return nil, fmt.Errorf("max gap %d is negative", meta.MaxGap)
	}
	sink, err := OpenSink(dir)
	if err != nil {
		return nil, err
	}
	if err := sink.Remove(FragmentMetaFile); err != nil {
		return nil, fmt.Errorf("clearing the previous %s in %s: %w", FragmentMetaFile, sink.Describe(), err)
	}
	meta.Chroms = nil
	meta.Complete = false
	meta.FormatVersion = FragmentFormatVersion
	return &FragmentWriter{sink: sink, meta: meta, seen: map[string]bool{}}, nil
}

// Write adds one block. Blocks must arrive in coordinate order: a chromosome's
// blocks together, each starting after the last one ended.
func (w *FragmentWriter) Write(b RefBlock) error {
	if w.done {
		return fmt.Errorf("fragment already finished")
	}
	if b.SampleID != w.meta.SampleID {
		return fmt.Errorf("block for %s written to %s's fragment", b.SampleID, w.meta.SampleID)
	}
	if b.End < b.Start {
		return fmt.Errorf("block %s:%d-%d ends before it starts", b.Chrom, b.Start, b.End)
	}
	if w.cur == nil || !SameChrom(w.cur.Name, b.Chrom) {
		if w.seen[CanonKey(b.Chrom)] {
			return fmt.Errorf("%s appears twice in the input, separated by another chromosome; "+
				"a fragment is written in one pass and needs each chromosome contiguous", b.Chrom)
		}
		if err := w.openChrom(b.Chrom); err != nil {
			return err
		}
	} else if b.Start <= w.ended {
		return fmt.Errorf("block %s:%d-%d starts at or before the previous block's end %d; "+
			"blocks must be disjoint and in coordinate order", b.Chrom, b.Start, b.End, w.ended)
	}
	if w.cur.Rows == 0 {
		w.cur.First = b.Start
	}
	w.cur.Rows++
	w.cur.Last = b.End
	w.ended = b.End
	w.buf = append(w.buf, b)
	if len(w.buf) >= batchSize {
		return w.flush()
	}
	return nil
}

func (w *FragmentWriter) openChrom(chrom string) error {
	if err := w.closeChrom(); err != nil {
		return err
	}
	name := fmt.Sprintf("%05d.parquet", len(w.meta.Chroms))
	f, err := w.sink.Create(name)
	if err != nil {
		return fmt.Errorf("creating %s in %s: %w", name, w.sink.Describe(), err)
	}
	w.files = append(w.files, name)
	w.out = f
	w.pw = parquet.NewGenericWriter[RefBlock](f,
		parquet.Compression(&zstd.Codec{}),
		parquet.SortingWriterConfig(parquet.SortingColumns(parquet.Ascending("start"))))
	w.meta.Chroms = append(w.meta.Chroms, FragmentChrom{Name: chrom, File: name})
	w.cur = &w.meta.Chroms[len(w.meta.Chroms)-1]
	w.seen[CanonKey(chrom)] = true
	w.ended = 0
	return nil
}

func (w *FragmentWriter) flush() error {
	if len(w.buf) == 0 {
		return nil
	}
	if _, err := w.pw.Write(w.buf); err != nil {
		return fmt.Errorf("writing %s: %w", w.cur.File, err)
	}
	w.buf = w.buf[:0]
	return nil
}

func (w *FragmentWriter) closeChrom() error {
	if w.pw == nil {
		return nil
	}
	if err := w.flush(); err != nil {
		return err
	}
	if err := w.pw.Close(); err != nil {
		return fmt.Errorf("finishing %s: %w", w.cur.File, err)
	}
	if err := w.out.Close(); err != nil {
		return fmt.Errorf("closing %s: %w", w.cur.File, err)
	}
	w.pw, w.out = nil, nil
	return nil
}

// Finish closes the last file and writes fragment.json, which is what makes the
// fragment complete. A sample with no block above the gate is a legitimate
// fragment of zero chromosomes: it says the gVCF was read and vouched for
// nothing, which is different from never having been read.
func (w *FragmentWriter) Finish() error {
	if w.done {
		return fmt.Errorf("fragment already finished")
	}
	if err := w.closeChrom(); err != nil {
		return err
	}
	w.done = true
	w.meta.Complete = true
	if w.meta.Created.IsZero() {
		w.meta.Created = time.Now().UTC()
	}
	if w.meta.Chroms == nil {
		w.meta.Chroms = []FragmentChrom{}
	}
	f, err := w.sink.Create(FragmentMetaFile)
	if err != nil {
		return fmt.Errorf("creating %s in %s: %w", FragmentMetaFile, w.sink.Describe(), err)
	}
	enc := json.NewEncoder(f)
	enc.SetIndent("", "  ")
	if err := enc.Encode(w.meta); err != nil {
		f.Close()
		return fmt.Errorf("encoding %s: %w", FragmentMetaFile, err)
	}
	if err := f.Close(); err != nil {
		return fmt.Errorf("writing %s to %s: %w", FragmentMetaFile, w.sink.Describe(), err)
	}
	return nil
}

// Meta returns what the fragment records about itself so far.
func (w *FragmentWriter) Meta() FragmentMeta { return w.meta }

// Discard removes everything the writer created. Files are removed WITHOUT
// being closed, for the reason Writer.Discard gives: closing would write valid
// footers into files about to disappear, and a failed unlink would then leave a
// well-formed partial fragment. There is no marker either way, so attach
// refuses what survives.
func (w *FragmentWriter) Discard() error {
	w.done = true
	var errs []error
	aborter, remote := w.sink.(Aborter)
	for _, name := range w.files {
		var err error
		if remote {
			err = aborter.Abort(name)
		} else {
			err = w.sink.Remove(name)
		}
		if err != nil {
			errs = append(errs, fmt.Errorf("removing %s: %w", name, err))
		}
	}
	if err := w.sink.Remove(FragmentMetaFile); err != nil {
		errs = append(errs, err)
	}
	return errors.Join(errs...)
}

// ReadFragmentMeta reads a fragment's sidecar, local or remote, refusing one that
// is not marked complete or comes from a newer layout.
func ReadFragmentMeta(ctx context.Context, dir string) (*FragmentMeta, error) {
	locator := joinStore(trimStoreDir(dir), FragmentMetaFile)
	src, err := iosource.Open(ctx, locator)
	if err != nil {
		return nil, fmt.Errorf("%s has no readable %s, so it is not a finished fragment: %w",
			dir, FragmentMetaFile, err)
	}
	defer src.Close()
	size, err := src.Size()
	if err != nil {
		return nil, err
	}
	raw, err := io.ReadAll(io.NewSectionReader(src, 0, size))
	if err != nil {
		return nil, fmt.Errorf("reading %s: %w", locator, err)
	}
	var m FragmentMeta
	if err := json.Unmarshal(raw, &m); err != nil {
		return nil, fmt.Errorf("parsing %s: %w", locator, err)
	}
	if !m.Complete {
		return nil, fmt.Errorf("%s is marked incomplete; rebuild the fragment", locator)
	}
	if m.FormatVersion > FragmentFormatVersion {
		return nil, fmt.Errorf("%s has fragment format %d; this build reads up to %d",
			locator, m.FormatVersion, FragmentFormatVersion)
	}
	return &m, nil
}

// FragmentBlocks streams one chromosome's blocks from a fragment, in order,
// until fn returns false.
func FragmentBlocks(ctx context.Context, dir string, c FragmentChrom, fn func(RefBlock) bool) error {
	t, err := openTable(ctx, joinStore(trimStoreDir(dir), c.File))
	if err != nil {
		return err
	}
	set := &shardSet{name: c.File, shards: []*shard{{m: t}}}
	defer set.Close()
	return scanParquet(set, fn)
}
