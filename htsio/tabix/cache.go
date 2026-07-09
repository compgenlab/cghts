package tabix

import (
	"container/list"
	"iter"
	"sort"
)

// windowShift keys the parsed-record cache by a 16 kb genomic window — the same
// granularity as the tabix linear index (see index.go). A position-sorted stream
// hits many consecutive variants in the same window, so a small LRU serves almost
// every query from memory without re-scanning/re-parsing the chunk.
const windowShift = 14

// DefaultRecordCacheWindows bounds the parsed-record LRU, in windows. Each entry is
// one 16 kb window's records (tens, occasionally hundreds). A sorted sweep needs
// only a couple live at once; the surplus absorbs dense overlapping features and the
// occasional backward query. Mirrors gtf.defaultGeneCacheCap in spirit.
const DefaultRecordCacheWindows = 256

// windowKey identifies a hydrated window: reference id + (start >> windowShift).
type windowKey struct {
	refID int
	win   int
}

type windowEntry struct {
	key  windowKey
	recs []*Record // every record overlapping the window, in file (Start-sorted) order
}

// recordCache is a bounded LRU of hydrated 16 kb-window record slices. It is NOT
// safe for concurrent use (the list/map are mutated on every query); build one per
// Reader and use it from a single goroutine, exactly as the pipeline does. Mirrors
// gtf.IndexedAnnotationSource.
type recordCache struct {
	cap   int
	ll    *list.List // MRU at front; values are *windowEntry
	items map[windowKey]*list.Element
}

func newRecordCache(cap int) *recordCache {
	if cap < 1 {
		cap = 1
	}
	return &recordCache{cap: cap, ll: list.New(), items: make(map[windowKey]*list.Element)}
}

func (c *recordCache) get(k windowKey) ([]*Record, bool) {
	el, ok := c.items[k]
	if !ok {
		return nil, false
	}
	c.ll.MoveToFront(el)
	return el.Value.(*windowEntry).recs, true
}

func (c *recordCache) put(k windowKey, recs []*Record) {
	if el, ok := c.items[k]; ok {
		c.ll.MoveToFront(el)
		el.Value.(*windowEntry).recs = recs
		return
	}
	el := c.ll.PushFront(&windowEntry{key: k, recs: recs})
	c.items[k] = el
	for c.ll.Len() > c.cap {
		back := c.ll.Back()
		if back == nil {
			break
		}
		c.ll.Remove(back)
		delete(c.items, back.Value.(*windowEntry).key)
	}
}

// sliceSeq yields the cached records overlapping [start, end) — rec.End > start &&
// rec.Start < end — matching iterChunks' overlap test. recs are Start-sorted, so a
// binary search bounds the scan by end; the End filter is a cheap linear pass over a
// single window's (few) records.
func sliceSeq(recs []*Record, start, end int) iter.Seq2[*Record, error] {
	return func(yield func(*Record, error) bool) {
		hi := sort.Search(len(recs), func(i int) bool { return recs[i].Start >= end })
		for i := 0; i < hi; i++ {
			if recs[i].End <= start {
				continue
			}
			if !yield(recs[i], nil) {
				return
			}
		}
	}
}
