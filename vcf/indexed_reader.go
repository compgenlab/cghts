package vcf

import (
	"container/list"
	"iter"
	"sort"

	"github.com/compgenlab/hts/htsio/tabix"
)

// vcfCacheWindows bounds the parsed-VcfRecord LRU, in 16 kb windows (matching the
// tabix linear index). Parsing a VCF record (INFO/FORMAT split) is expensive, so
// caching the decoded records — not just the tabix line — means a position-sorted
// stream parses each source record once.
const vcfCacheWindows = 256

// vcfWindowShift matches tabix's linear-index window (16 kb). Kept in sync with
// tabix.windowShift; a query contained in one window is served from the cache.
const vcfWindowShift = 14

// IndexedVcfReader provides random access to a tabix-indexed VCF file
// (BGZF-compressed with a companion .tbi or .csi index).
//
// It caches decoded records per 16 kb window and disables the underlying tabix
// reader's own record cache to avoid holding each window twice. Not safe for
// concurrent use; build one per annotation pass.
type IndexedVcfReader struct {
	tr       *tabix.Reader
	filename string
	header   *VcfHeader
	cache    *vcfRecordCache // nil = disabled
}

// NewIndexedVcfReader opens a tabix-indexed VCF file for random access. The file
// must be BGZF-compressed and have a companion .tbi or .csi index. The caller
// should [IndexedVcfReader.Close] the reader when done.
func NewIndexedVcfReader(filename string) (*IndexedVcfReader, error) {
	tr, err := tabix.NewReader(filename)
	if err != nil {
		return nil, err
	}
	tr.DisableCache() // this reader caches decoded VcfRecords itself (below)
	return &IndexedVcfReader{tr: tr, filename: filename, cache: newVcfRecordCache(vcfCacheWindows)}, nil
}

// Header returns the VCF header, reading it from the start of the BGZF stream on
// first call.
func (r *IndexedVcfReader) Header() (*VcfHeader, error) {
	if r.header == nil {
		hr, err := NewVcfFile(r.filename)
		if err != nil {
			return nil, err
		}
		defer hr.Close()
		h, err := hr.Header()
		if err != nil {
			return nil, err
		}
		r.header = h
	}
	return r.header, nil
}

// Query returns an iterator over the VCF records overlapping the 0-based
// half-open region [start, end) on the given reference. The iterator yields
// (nil, err) and stops if a record line cannot be parsed.
func (r *IndexedVcfReader) Query(ref string, start, end int) (iter.Seq2[*VcfRecord, error], error) {
	header, err := r.Header()
	if err != nil {
		return nil, err
	}

	// Cache path: a query contained in a single 16 kb window is served from — or
	// hydrates — the decoded-record cache, so newRecord runs once per window. The
	// window's records cover every sub-query within it. Wider/edge queries fall
	// through to the streaming parse and are not cached.
	if r.cache != nil && start >= 0 && start < end && (start>>vcfWindowShift) == ((end-1)>>vcfWindowShift) {
		win := start >> vcfWindowShift
		recs, ok := r.cache.get(vcfWinKey{ref, win})
		if !ok {
			winStart := win << vcfWindowShift
			recs, err = r.hydrateWindow(ref, winStart, winStart+(1<<vcfWindowShift), header)
			if err != nil {
				return nil, err
			}
			r.cache.put(vcfWinKey{ref, win}, recs)
		}
		return vcfSliceSeq(recs, start, end), nil
	}

	recs, err := r.tr.Query(ref, start, end)
	if err != nil {
		return nil, err
	}
	return func(yield func(*VcfRecord, error) bool) {
		for rec, err := range recs {
			if err != nil {
				yield(nil, err)
				return
			}
			vr, perr := newRecord(rec.Line, header)
			if perr != nil {
				yield(nil, perr)
				return
			}
			if !yield(vr, nil) {
				return
			}
		}
	}, nil
}

// hydrateWindow decodes every record overlapping the 16 kb window once, pairing each
// with the tabix Start/End used for overlap tests (so filtering matches the tabix
// query exactly, including symbolic/SV END handling).
func (r *IndexedVcfReader) hydrateWindow(ref string, winStart, winEnd int, header *VcfHeader) ([]cachedVcf, error) {
	seq, err := r.tr.Query(ref, winStart, winEnd)
	if err != nil {
		return nil, err
	}
	var out []cachedVcf
	for rec, err := range seq {
		if err != nil {
			return nil, err
		}
		vr, perr := newRecord(rec.Line, header)
		if perr != nil {
			return nil, perr
		}
		out = append(out, cachedVcf{vr: vr, start: rec.Start, end: rec.End})
	}
	return out, nil
}

// HasRef reports whether the index contains the given reference. Use it to skip
// a query that would otherwise fail for a reference absent from this file.
func (r *IndexedVcfReader) HasRef(ref string) bool {
	return r.tr.HasRef(ref)
}

// RefNames returns the contig names present in the index, in reference order.
// Pass the result to [vcf.NewContigConverter] to match query contigs against
// this file's naming scheme.
func (r *IndexedVcfReader) RefNames() []string {
	return r.tr.RefNames()
}

// Close releases resources held by the reader.
func (r *IndexedVcfReader) Close() error {
	return r.tr.Close()
}

// cachedVcf is a decoded record plus the tabix 0-based [start, end) it was indexed
// at, kept so window sub-queries filter overlap identically to a tabix query.
type cachedVcf struct {
	vr    *VcfRecord
	start int
	end   int
}

type vcfWinKey struct {
	ref string
	win int
}

type vcfWinEntry struct {
	key  vcfWinKey
	recs []cachedVcf
}

// vcfRecordCache is a bounded LRU of decoded 16 kb-window records. Not safe for
// concurrent use; one per reader.
type vcfRecordCache struct {
	cap   int
	ll    *list.List
	items map[vcfWinKey]*list.Element
}

func newVcfRecordCache(cap int) *vcfRecordCache {
	if cap < 1 {
		cap = 1
	}
	return &vcfRecordCache{cap: cap, ll: list.New(), items: make(map[vcfWinKey]*list.Element)}
}

func (c *vcfRecordCache) get(k vcfWinKey) ([]cachedVcf, bool) {
	el, ok := c.items[k]
	if !ok {
		return nil, false
	}
	c.ll.MoveToFront(el)
	return el.Value.(*vcfWinEntry).recs, true
}

func (c *vcfRecordCache) put(k vcfWinKey, recs []cachedVcf) {
	if el, ok := c.items[k]; ok {
		c.ll.MoveToFront(el)
		el.Value.(*vcfWinEntry).recs = recs
		return
	}
	el := c.ll.PushFront(&vcfWinEntry{key: k, recs: recs})
	c.items[k] = el
	for c.ll.Len() > c.cap {
		back := c.ll.Back()
		if back == nil {
			break
		}
		c.ll.Remove(back)
		delete(c.items, back.Value.(*vcfWinEntry).key)
	}
}

// vcfSliceSeq yields cached records overlapping [start, end) — end > start &&
// start < end on the tabix coords — matching the tabix query's overlap test.
func vcfSliceSeq(recs []cachedVcf, start, end int) iter.Seq2[*VcfRecord, error] {
	return func(yield func(*VcfRecord, error) bool) {
		hi := sort.Search(len(recs), func(i int) bool { return recs[i].start >= end })
		for i := 0; i < hi; i++ {
			if recs[i].end <= start {
				continue
			}
			if !yield(recs[i].vr, nil) {
				return
			}
		}
	}
}
