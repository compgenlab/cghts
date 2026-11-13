package varstore

import (
	"container/heap"
	"context"
	"errors"
	"fmt"
	"io"
	"iter"
	"math"
	"os"
	"path/filepath"
	"slices"
	"sort"
	"strings"
	"time"

	"github.com/parquet-go/parquet-go"
	"github.com/parquet-go/parquet-go/compress/zstd"
)

// Attaching reference-block fragments to a volume.
//
// ONE MERGE, TWO USES. A volume's ref_blocks table is sorted by position across
// every sample, so one person's rows are interleaved with everyone's in every
// file. Any change -- building the table, or swapping one person's blocks after
// their gVCF is redone -- therefore rewrites every file that person appears in,
// which is nearly all of them. What differs is only what is read:
//
//   - REBUILD reads fragments alone. The fragments are the source of truth, and
//     the table is replaced wholesale. This is the first attach, and the way to
//     redo one from scratch.
//   - REPLACE reads the existing table with chosen samples filtered out, plus
//     fragments for those samples. The table stands in for every fragment not
//     being changed, so they need not have been kept.
//
// Both are the same k-way merge over sorted inputs; the existing table is one
// more input. That keeps one write path, one set of checks and one manifest
// rewrite, instead of a second implementation to keep correct.
//
// NOTHING IS OVERWRITTEN IN PLACE. Output goes to a new generation of files, and
// the manifest -- rewritten last -- is what switches the volume over. Until it
// lands the volume reads exactly as before, so a run that dies costs nothing but
// orphaned files, and a replace can read the table it is replacing. The previous
// generation is removed once the new manifest is written.
//
// SHARDS FOLLOW THE CATALOG. On a sharded volume the output is cut where the
// catalog's shards begin, and a block crossing a cut is split, so each file stays
// bounded and covers one coordinate range. A ref_blocks shard's First/Last are
// the extent of the blocks in it -- a partition of the chromosome -- unlike the
// catalog's, which are site ranges.

// ErrRefBlocksExist refuses to attach over a volume that already has reference
// blocks without being told whether to rebuild or to replace.
var ErrRefBlocksExist = errors.New("volume already has ref_blocks")

// DefaultAttachBatch bounds how many inputs one merge holds open at once.
const DefaultAttachBatch = 100

// Fragment is a fragment's directory and its sidecar, read once and reused for
// every volume it is attached to.
type Fragment struct {
	Dir  string
	Meta FragmentMeta
}

// LoadFragments reads fragment sidecars and checks they belong together: each
// complete, no sample twice, and all built the same way. Blocks banded, gated or
// gap-merged differently do not mean the same thing by "covered", and a table
// mixing them would answer one gate differently for different people.
func LoadFragments(ctx context.Context, dirs []string) ([]Fragment, error) {
	out := make([]Fragment, 0, len(dirs))
	bySample := map[string]string{}
	for _, d := range dirs {
		m, err := ReadFragmentMeta(ctx, d)
		if err != nil {
			return nil, err
		}
		if prev, dup := bySample[m.SampleID]; dup {
			return nil, fmt.Errorf("%s appears in two fragments, %s and %s; a sample's blocks "+
				"must come from one source, or its coverage would be counted twice", m.SampleID, prev, d)
		}
		bySample[m.SampleID] = d
		if len(out) > 0 {
			if err := sameBuild(out[0], Fragment{Dir: d, Meta: *m}); err != nil {
				return nil, err
			}
		}
		out = append(out, Fragment{Dir: d, Meta: *m})
	}
	return out, nil
}

func sameBuild(a, b Fragment) error {
	switch {
	case a.Meta.MinDP != b.Meta.MinDP:
		return fmt.Errorf("%s was built at --min-dp %d and %s at %d", a.Dir, a.Meta.MinDP, b.Dir, b.Meta.MinDP)
	case !slices.Equal(a.Meta.DepthBands, b.Meta.DepthBands):
		return fmt.Errorf("%s was banded %v and %s %v", a.Dir, a.Meta.DepthBands, b.Dir, b.Meta.DepthBands)
	case a.Meta.MaxGap != b.Meta.MaxGap:
		return fmt.Errorf("%s merged gaps up to %d and %s up to %d; their blocks claim different things",
			a.Dir, a.Meta.MaxGap, b.Dir, b.Meta.MaxGap)
	case a.Meta.Passthrough != b.Meta.Passthrough:
		return fmt.Errorf("%s and %s disagree on --passthrough", a.Dir, b.Dir)
	}
	return nil
}

// AttachOpts says what an attach reads.
type AttachOpts struct {
	// Rebuild replaces any existing ref_blocks with the fragments alone.
	Rebuild bool
	// Replace names samples whose existing rows are dropped; fragments for them,
	// if given, take their place. A replaced sample with no fragment is removed.
	Replace []string
	// Batch bounds how many inputs are merged at once. Zero is DefaultAttachBatch.
	Batch int
	// TempDir holds intermediate merge files; empty is the system default. Always
	// local: intermediates are written and read once, and never belong remotely.
	TempDir string
	// ExtraChroms are chromosomes to take from the fragments beyond those the
	// volume already holds. Each that gains a block is added to the volume's
	// census, which is what lets a store route a query there. Set by
	// AttachRefBlocksToStore; see it for why.
	ExtraChroms []string
}

// AttachResult reports what an attach wrote.
type AttachResult struct {
	Rows    int64
	Shards  int
	Samples int // samples with at least one block in the new table
}

// AttachRefBlocks merges fragments into one volume's ref_blocks.
//
// The fragments are genome-wide; only their blocks on the volume's own
// chromosomes are read. Every volume gets a table even if no fragment has a block
// on it, because the table's presence is the claim "these samples were checked
// here" -- an empty one says nobody was covered, where a missing one would say
// nobody looked, and an archive needs every volume to say the same.
func AttachRefBlocks(ctx context.Context, base string, frags []Fragment, opts AttachOpts) (*AttachResult, error) {
	base = TrimStoreSuffix(base)
	man, err := ReadVolumeManifestContext(ctx, base)
	if err != nil {
		return nil, err
	}
	if !man.Complete {
		return nil, fmt.Errorf("%s is not a finished volume", base)
	}

	roster := make(map[string]bool, len(man.Samples))
	for _, s := range man.Samples {
		roster[s] = true
	}
	for _, f := range frags {
		if !roster[f.Meta.SampleID] {
			return nil, fmt.Errorf("fragment %s is for %s, who is not among %s's samples",
				f.Dir, f.Meta.SampleID, base)
		}
		if f.Meta.MinDP != man.Params.MinDP {
			return nil, fmt.Errorf("fragment %s was gated at --min-dp %d but %s at %d; blocks and "+
				"catalog runs would answer a gate differently", f.Dir, f.Meta.MinDP, base, man.Params.MinDP)
		}
		if !slices.Equal(f.Meta.DepthBands, man.Params.DepthBands) {
			return nil, fmt.Errorf("fragment %s was banded %v but %s %v", f.Dir, f.Meta.DepthBands,
				base, man.Params.DepthBands)
		}
	}
	for i := 1; i < len(frags); i++ {
		if err := sameBuild(frags[0], frags[i]); err != nil {
			return nil, err
		}
	}

	existing := man.Params.RefBlocks
	keepExisting := false
	switch {
	case existing && !opts.Rebuild && len(opts.Replace) == 0:
		return nil, fmt.Errorf("%w: %s; rebuild it from fragments, or name the samples to replace",
			ErrRefBlocksExist, base)
	case !existing && len(opts.Replace) > 0:
		return nil, fmt.Errorf("%s has no ref_blocks, so there is nothing to replace; attach to build it", base)
	case opts.Rebuild && len(opts.Replace) > 0:
		return nil, fmt.Errorf("rebuild and replace are alternatives: a rebuild already replaces everyone")
	case len(opts.Replace) > 0:
		keepExisting = true
	}
	maxGap := man.Params.MaxGap
	if len(frags) > 0 {
		maxGap = frags[0].Meta.MaxGap
	}
	replaced := map[string]bool{}
	if keepExisting {
		for _, s := range opts.Replace {
			replaced[s] = true
		}
		for _, f := range frags {
			if !replaced[f.Meta.SampleID] {
				return nil, fmt.Errorf("fragment %s is for %s, who is not being replaced; their "+
					"existing blocks would be kept beside the new ones", f.Dir, f.Meta.SampleID)
			}
		}
		if len(frags) > 0 && maxGap != man.Params.MaxGap {
			return nil, fmt.Errorf("the fragments merged gaps up to %d but %s's blocks up to %d; "+
				"rebuild to change it", maxGap, base, man.Params.MaxGap)
		}
	}

	var vol *ParquetVolume
	if keepExisting {
		if vol, err = OpenParquetContext(ctx, base); err != nil {
			return nil, err
		}
		defer vol.Close()
	}

	sink, err := OpenSink(base)
	if err != nil {
		return nil, err
	}
	tmp, err := os.MkdirTemp(opts.TempDir, "refblocks-merge-")
	if err != nil {
		return nil, err
	}
	defer os.RemoveAll(tmp)

	batch := opts.Batch
	switch {
	case batch <= 0:
		batch = DefaultAttachBatch
	case batch < 2:
		batch = 2 // a batch of one merges nothing, and the tree would never shrink
	}
	gen := fmt.Sprintf("%x", time.Now().UnixNano())
	out := &segmentWriter{sink: sink, gen: gen}
	samples := map[string]bool{}

	chroms := make([]string, 0, len(man.Chromosomes)+len(opts.ExtraChroms))
	held := map[string]bool{}
	for _, c := range man.Chromosomes {
		chroms = append(chroms, c.Name)
		held[CanonKey(c.Name)] = true
	}
	var added []string
	for _, c := range opts.ExtraChroms {
		if !held[CanonKey(c)] {
			held[CanonKey(c)] = true
			chroms = append(chroms, c)
			added = append(added, c)
		}
	}

	for _, name := range chroms {
		chrom := ChromCensus{Name: name}
		var inputs []blockSource
		for _, f := range frags {
			for _, c := range f.Meta.Chroms {
				if SameChrom(c.Name, chrom.Name) {
					f, c := f, c
					inputs = append(inputs, func() iter.Seq2[RefBlock, error] {
						return fileBlocks(func(fn func(RefBlock) bool) error {
							return FragmentBlocks(ctx, f.Dir, c, fn)
						})
					})
				}
			}
		}
		if vol != nil && vol.hasRefBlocks {
			name := chrom.Name
			inputs = append(inputs, func() iter.Seq2[RefBlock, error] {
				return fileBlocks(func(fn func(RefBlock) bool) error {
					// Every block on the chromosome: a shard filter by name only,
					// since a position bound would prune everything.
					return scanParquetPruned(vol.refBlocks, scanFilter{
						shard:    func(si ShardInfo) bool { return SameChrom(si.Chrom, name) },
						rowGroup: keepAll,
					},
						func(b RefBlock) bool {
							if !SameChrom(b.Chrom, name) || replaced[b.SampleID] {
								return true
							}
							return fn(b)
						})
				})
			})
		}
		merged, err := mergeTree(inputs, batch, tmp)
		if err != nil {
			out.discard()
			return nil, err
		}
		if err := out.chromosome(chrom.Name, catalogCuts(man, chrom.Name), merged, samples); err != nil {
			out.discard()
			return nil, err
		}
	}

	// The switch. Everything before this left the volume reading as it did.
	old := man.Tables[RefBlocksTable]
	man.Tables[RefBlocksTable] = TableInfo{Rows: out.rows, Bytes: out.bytes, Shards: out.shards}
	// A chromosome taken in for its blocks joins the census, which records rows
	// rather than intent: one that gained none stays out, since a census entry
	// is what routes a query to this volume.
	for _, name := range added {
		var c *ChromCensus
		for _, sh := range out.shards {
			if !SameChrom(sh.Chrom, name) {
				continue
			}
			if c == nil {
				man.Chromosomes = append(man.Chromosomes, ChromCensus{Name: name, FirstPos: sh.First, LastPos: sh.Last})
				c = &man.Chromosomes[len(man.Chromosomes)-1]
			}
			c.FirstPos = min(c.FirstPos, sh.First)
			c.LastPos = max(c.LastPos, sh.Last)
		}
	}
	man.Params.RefBlocks = true
	man.Params.MaxGap = maxGap
	if err := WriteManifestTo(sink, *man); err != nil {
		out.discard()
		return nil, err
	}

	// The previous generation, now unreferenced. Best effort: a leftover file is
	// unreachable, where failing here would report an attach that succeeded as
	// one that did not.
	if existing {
		if len(old.Shards) == 0 {
			_ = sink.Remove(TableFile(RefBlocksTable))
		}
		for _, sh := range old.Shards {
			_ = sink.Remove(sh.Name)
		}
	}
	return &AttachResult{Rows: out.rows, Shards: len(out.shards), Samples: len(samples)}, nil
}

// catalogCuts returns where a chromosome's ref_blocks files are cut: the first
// site of every catalog shard after the first. Unsharded, there are none.
func catalogCuts(man *VolumeManifest, chrom string) []int32 {
	var firsts []int32
	for _, sh := range man.Tables[CatalogTable].Shards {
		if SameChrom(sh.Chrom, chrom) {
			firsts = append(firsts, sh.First)
		}
	}
	sort.Slice(firsts, func(i, j int) bool { return firsts[i] < firsts[j] })
	if len(firsts) > 0 {
		firsts = firsts[1:]
	}
	return firsts
}

// blockSource opens one sorted stream of blocks. Opening is deferred so a merge
// holds only its own batch of inputs open.
type blockSource func() iter.Seq2[RefBlock, error]

// fileBlocks adapts a callback scan to an iterator that carries its error.
func fileBlocks(scan func(func(RefBlock) bool) error) iter.Seq2[RefBlock, error] {
	return func(yield func(RefBlock, error) bool) {
		stopped := false
		err := scan(func(b RefBlock) bool {
			if !yield(b, nil) {
				stopped = true
				return false
			}
			return true
		})
		if err != nil && !stopped {
			yield(RefBlock{}, err)
		}
	}
}

// mergeTree merges any number of sorted inputs, at most batch at a time: a larger
// set is merged in batches into local intermediate files, which are merged in
// turn. Ten thousand fragments therefore never hold ten thousand files open.
func mergeTree(inputs []blockSource, batch int, tmp string) (iter.Seq2[RefBlock, error], error) {
	level := 0
	for len(inputs) > batch {
		var next []blockSource
		for i := 0; i < len(inputs); i += batch {
			group := inputs[i:min(i+batch, len(inputs))]
			path := filepath.Join(tmp, fmt.Sprintf("l%d-%05d.parquet", level, i/batch))
			if err := writeIntermediate(path, mergeSorted(group)); err != nil {
				return nil, err
			}
			next = append(next, func() iter.Seq2[RefBlock, error] {
				return fileBlocks(func(fn func(RefBlock) bool) error {
					return scanLocalBlocks(path, fn)
				})
			})
		}
		inputs = next
		level++
	}
	return mergeSorted(inputs), nil
}

func writeIntermediate(path string, blocks iter.Seq2[RefBlock, error]) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	w := newRefBlockWriter(f)
	buf := make([]RefBlock, 0, batchSize)
	for b, err := range blocks {
		if err != nil {
			f.Close()
			return err
		}
		buf = append(buf, b)
		if len(buf) == batchSize {
			if _, err := w.Write(buf); err != nil {
				f.Close()
				return err
			}
			buf = buf[:0]
		}
	}
	if _, err := w.Write(buf); err != nil {
		f.Close()
		return err
	}
	if err := w.Close(); err != nil {
		f.Close()
		return err
	}
	return f.Close()
}

func scanLocalBlocks(path string, fn func(RefBlock) bool) error {
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()
	r := parquet.NewGenericReader[RefBlock](f)
	defer r.Close()
	buf := make([]RefBlock, batchSize)
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
	}
}

func newRefBlockWriter(w io.Writer) *parquet.GenericWriter[RefBlock] {
	return parquet.NewGenericWriter[RefBlock](w,
		parquet.Compression(&zstd.Codec{}),
		parquet.SortingWriterConfig(parquet.SortingColumns(
			parquet.Ascending("chrom"), parquet.Ascending("start"))),
		parquet.BloomFilters(parquet.SplitBlockFilter(10, "sample_id")))
}

// blockLess is the table's order: position, then sample, so the result does not
// depend on which input a block came from.
func blockLess(a, b RefBlock) bool {
	if a.Start != b.Start {
		return a.Start < b.Start
	}
	return a.SampleID < b.SampleID
}

type mergeHead struct {
	b    RefBlock
	next func() (RefBlock, error, bool)
	stop func()
}

type mergeHeap []*mergeHead

func (h mergeHeap) Len() int           { return len(h) }
func (h mergeHeap) Less(i, j int) bool { return blockLess(h[i].b, h[j].b) }
func (h mergeHeap) Swap(i, j int)      { h[i], h[j] = h[j], h[i] }
func (h *mergeHeap) Push(x any)        { *h = append(*h, x.(*mergeHead)) }
func (h *mergeHeap) Pop() any {
	old := *h
	x := old[len(old)-1]
	*h = old[:len(old)-1]
	return x
}

// mergeSorted is a k-way merge of sorted inputs. An input out of order is an
// error rather than something to sort around: every input here was written
// sorted, so disorder means a damaged or foreign file.
func mergeSorted(inputs []blockSource) iter.Seq2[RefBlock, error] {
	return func(yield func(RefBlock, error) bool) {
		h := &mergeHeap{}
		defer func() {
			for _, hd := range *h {
				hd.stop()
			}
		}()
		for _, open := range inputs {
			next, stop := iter.Pull2(open())
			b, err, ok := next()
			if err != nil {
				stop()
				yield(RefBlock{}, err)
				return
			}
			if !ok {
				stop()
				continue
			}
			heap.Push(h, &mergeHead{b: b, next: next, stop: stop})
		}
		for h.Len() > 0 {
			hd := (*h)[0]
			cur := hd.b
			if !yield(cur, nil) {
				return
			}
			b, err, ok := hd.next()
			switch {
			case err != nil:
				yield(RefBlock{}, err)
				return
			case !ok:
				hd.stop()
				heap.Pop(h)
			default:
				if blockLess(b, cur) {
					yield(RefBlock{}, fmt.Errorf("an input is out of order: %s %s:%d after %s:%d",
						b.SampleID, b.Chrom, b.Start, cur.Chrom, cur.Start))
					return
				}
				hd.b = b
				heap.Fix(h, 0)
			}
		}
	}
}

// segmentWriter writes one volume's new ref_blocks generation, a file per
// chromosome segment, and records the shard index for it.
type segmentWriter struct {
	sink   Sink
	gen    string
	shards []ShardInfo
	rows   int64
	bytes  int64
	files  []string

	out io.WriteCloser
	pw  *parquet.GenericWriter[RefBlock]
	cur *ShardInfo
	buf []RefBlock
}

// chromosome writes one chromosome's merged blocks, cut at the given positions.
// A block crossing a cut is split, each piece keeping the block's depth and
// quality, so every file holds only blocks inside its own range.
func (w *segmentWriter) chromosome(chrom string, cuts []int32, blocks iter.Seq2[RefBlock, error], samples map[string]bool) error {
	// Segment i covers [cuts[i-1], hi(i)], from position 1 for the first.
	hi := func(i int) int32 {
		if i == len(cuts) {
			return math.MaxInt32
		}
		return cuts[i] - 1
	}
	seg := 0
	var pending []RefBlock // remainders of split blocks, all starting at hi(seg)+1

	var place func(b RefBlock) error
	place = func(b RefBlock) error {
		if b.End > hi(seg) {
			rest := b
			rest.Start = hi(seg) + 1
			pending = append(pending, rest)
			b.End = hi(seg)
		}
		return w.write(chrom, b)
	}
	advance := func() error {
		if err := w.closeSegment(); err != nil {
			return err
		}
		seg++
		carried := pending
		pending = nil
		sort.Slice(carried, func(i, j int) bool { return carried[i].SampleID < carried[j].SampleID })
		for _, b := range carried {
			if err := place(b); err != nil {
				return err
			}
		}
		return nil
	}

	for b, err := range blocks {
		if err != nil {
			return err
		}
		samples[b.SampleID] = true
		for b.Start > hi(seg) {
			if err := advance(); err != nil {
				return err
			}
		}
		if err := place(b); err != nil {
			return err
		}
	}
	for len(pending) > 0 {
		if err := advance(); err != nil {
			return err
		}
	}
	return w.closeSegment()
}

func (w *segmentWriter) write(chrom string, b RefBlock) error {
	b.Chrom = chrom // the volume's own spelling, whatever the fragment used
	if w.pw == nil {
		name := fmt.Sprintf("%s/%s-%05d.parquet", RefBlocksTable, w.gen, len(w.shards))
		f, err := w.sink.Create(name)
		if err != nil {
			return fmt.Errorf("creating %s: %w", name, err)
		}
		w.files = append(w.files, name)
		w.out, w.pw = f, newRefBlockWriter(f)
		w.shards = append(w.shards, ShardInfo{Name: name, Chrom: chrom, First: b.Start, Last: b.End})
		w.cur = &w.shards[len(w.shards)-1]
	}
	if b.Start < w.cur.First {
		w.cur.First = b.Start
	}
	if b.End > w.cur.Last {
		w.cur.Last = b.End
	}
	w.cur.Rows++
	w.rows++
	w.buf = append(w.buf, b)
	if len(w.buf) == batchSize {
		if _, err := w.pw.Write(w.buf); err != nil {
			return err
		}
		w.buf = w.buf[:0]
	}
	return nil
}

func (w *segmentWriter) closeSegment() error {
	if w.pw == nil {
		return nil
	}
	if _, err := w.pw.Write(w.buf); err != nil {
		return err
	}
	w.buf = w.buf[:0]
	if err := w.pw.Close(); err != nil {
		return err
	}
	if err := w.out.Close(); err != nil {
		return err
	}
	if size, ok, err := w.sink.Stat(w.cur.Name); err == nil && ok {
		w.cur.Bytes = size
		w.bytes += size
	}
	w.pw, w.out, w.cur = nil, nil, nil
	return nil
}

// discard removes a generation that will never be referenced. Best effort: the
// manifest was not switched, so whatever survives is unreachable.
func (w *segmentWriter) discard() {
	aborter, remote := w.sink.(Aborter)
	for _, name := range w.files {
		if remote {
			_ = aborter.Abort(name)
		} else {
			_ = w.sink.Remove(name)
		}
	}
}

// StoreAttachResult reports an attach across a store.
type StoreAttachResult struct {
	Volumes map[string]*AttachResult
	// Created are blocks-only volumes made for chromosomes no volume held.
	Created []string
}

// AttachRefBlocksToStore attaches fragments to a store -- an archive of volumes
// or a single volume -- so that EVERY chromosome a fragment has blocks on is
// answerable, not only those with a catalog site.
//
// WHY THIS IS MORE THAN ATTACHING TO EACH VOLUME. A store routes a query by
// chromosome to the one volume that holds it, and a volume holds the
// chromosomes it has rows on -- until now, catalog sites. A chromosome with no
// catalog site (chrM, unplaced contigs, chrY in a cohort of women) therefore had
// no volume, and its blocks had nowhere to be read from: every query there came
// back not assayed, however well covered it was. Putting every chromosome into
// every volume is not the answer either -- each pVCF header declares every
// contig, so each block would land in every volume and a locus would have as
// many answers as there are volumes.
//
// So each such chromosome gets a home. On an archive it is a new BLOCKS-ONLY
// volume -- the store's samples and settings, an empty catalog, and the
// chromosome's blocks -- registered in the archive manifest, which is rewritten
// last. On a single volume, the volume itself takes the chromosome. Either way
// the chromosome joins a census, so it routes, and a query there is answered
// from blocks like any other off-catalog position. A declared contig nobody has
// a block on gets no volume: "no volume" and "nothing covered" both answer not
// assayed.
func AttachRefBlocksToStore(ctx context.Context, base string, frags []Fragment, opts AttachOpts, program, command string) (*StoreAttachResult, error) {
	base = TrimStoreSuffix(base)
	res := &StoreAttachResult{Volumes: map[string]*AttachResult{}}
	fragChroms := fragmentChroms(frags)

	if !IsStore(ctx, base) {
		man, err := ReadVolumeManifestContext(ctx, base)
		if err != nil {
			return nil, err
		}
		o := opts
		o.ExtraChroms = missingChroms(fragChroms, censusNames(man))
		r, err := AttachRefBlocks(ctx, base, frags, o)
		if err != nil {
			return nil, err
		}
		res.Volumes[base] = r
		return res, nil
	}

	smf, err := readStoreManifest(ctx, base)
	if err != nil {
		return nil, err
	}
	var names, held []string
	for _, v := range smf.Volumes {
		names = append(names, v.Name)
		held = append(held, v.Chroms...)
	}
	orphans := missingChroms(fragChroms, held)

	// Template for new volumes: the first volume's settings, which the archive
	// already requires every volume to share.
	var tmpl *VolumeManifest
	if len(names) > 0 {
		if tmpl, err = ReadVolumeManifestContext(ctx, volumeAt(base, names[0])); err != nil {
			return nil, err
		}
	}

	for _, n := range names {
		r, err := AttachRefBlocks(ctx, volumeAt(base, n), frags, opts)
		if err != nil {
			return nil, fmt.Errorf("volume %s: %w", n, err)
		}
		res.Volumes[n] = r
	}
	for _, chrom := range orphans {
		if tmpl == nil {
			break
		}
		name := volumeDirName(chrom)
		loc := volumeAt(base, name)
		if _, err := ReadVolumeManifestContext(ctx, loc); err == nil {
			return nil, fmt.Errorf("%s already holds a volume named %s that the archive does not list; "+
				"remove it or add it with BuildStore before attaching", base, name)
		}
		if err := writeEmptyVolume(loc, tmpl, program, command); err != nil {
			return nil, fmt.Errorf("creating a blocks-only volume for %s: %w", chrom, err)
		}
		// A volume made now has no blocks to keep or replace: it is built.
		r, err := AttachRefBlocks(ctx, loc, frags, AttachOpts{
			Batch: opts.Batch, TempDir: opts.TempDir, ExtraChroms: []string{chrom},
		})
		if err != nil {
			return nil, fmt.Errorf("volume %s: %w", name, err)
		}
		res.Volumes[name] = r
		res.Created = append(res.Created, name)
		names = append(names, name)
	}

	// Last: the archive manifest, so the archive and its volumes agree that
	// blocks are present, and so the new volumes are routed.
	man, err := BuildStore(ctx, base, names, program, command)
	if err != nil {
		return nil, err
	}
	if err := WriteStoreManifest(base, *man); err != nil {
		return nil, err
	}
	return res, nil
}

func volumeAt(base, name string) string {
	if isLocator(name) || strings.HasPrefix(name, "/") {
		return name
	}
	return joinStore(base, name)
}

// isLocator reports whether a volume name is a locator in its own right.
func isLocator(name string) bool { return strings.Contains(name, "://") }

// volumeDirName turns a chromosome name into a directory name. Contig names may
// hold characters a path should not ("HLA-A*01:01", "chrUn_KI270302v1").
func volumeDirName(chrom string) string {
	b := []byte(chrom)
	for i, c := range b {
		ok := c == '.' || c == '_' || c == '-' ||
			(c >= '0' && c <= '9') || (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')
		if !ok {
			b[i] = '_'
		}
	}
	return string(b)
}

// writeEmptyVolume writes a volume with the template's samples and settings and
// no rows -- what an archive requires of a volume, so that it agrees.
func writeEmptyVolume(loc string, tmpl *VolumeManifest, program, command string) error {
	w, err := NewWriter(loc, WriterOpts{
		Samples:    tmpl.Samples,
		MinDP:      tmpl.Params.MinDP,
		NoCallable: tmpl.Params.NoCallable,
		DepthBands: tmpl.Params.DepthBands,
		DepthRatio: tmpl.Params.DepthRatio,
		Info:       tmpl.Params.Info,
		Format:     tmpl.Params.Format,
		Meta:       tmpl.Meta,
		Program:    program,
		Command:    command,
	})
	if err != nil {
		return err
	}
	if err := w.Finish(); err != nil {
		w.Discard()
		return err
	}
	return nil
}

// fragmentChroms is every chromosome any fragment has a block on, in the order
// first seen, spelled as first seen.
func fragmentChroms(frags []Fragment) []string {
	var out []string
	seen := map[string]bool{}
	for _, f := range frags {
		for _, c := range f.Meta.Chroms {
			if c.Rows > 0 && !seen[CanonKey(c.Name)] {
				seen[CanonKey(c.Name)] = true
				out = append(out, c.Name)
			}
		}
	}
	return out
}

func censusNames(m *VolumeManifest) []string {
	out := make([]string, 0, len(m.Chromosomes))
	for _, c := range m.Chromosomes {
		out = append(out, c.Name)
	}
	return out
}

func missingChroms(want, have []string) []string {
	h := map[string]bool{}
	for _, c := range have {
		h[CanonKey(c)] = true
	}
	var out []string
	for _, c := range want {
		if !h[CanonKey(c)] {
			out = append(out, c)
		}
	}
	return out
}
