package varstore

import "fmt"

// Cutting a store into coordinate shards as it is written.
//
// A shard is closed when the site count reaches ShardSites, and always at a
// chromosome change -- so a shard never spans chromosomes, which is what lets
// its First and Last be compared without a chromosome test on every row.
//
// THE THREE MEMBERS ARE CUT TOGETHER, on the same site boundaries. Shard k of
// calls, sites and regions therefore cover the same interval, and a locus query
// reads the k'th of each rather than reconciling two partitionings. It also
// means the bounds recorded for all three are the SITE range: a run is bounded
// by the sites it covers, and a call sits at one.

// noteSite advances the shard bookkeeping for a site about to be written, and
// rotates first when this site belongs in the next shard.
//
// Called before the site is buffered, so the rotation lands between sites
// rather than inside one.
func (w *Writer) noteSite(chrom string, pos int32) error {
	if w.opts.ShardSites > 0 && w.shardSites > 0 {
		if !SameChrom(chrom, w.shardChrom) || w.shardSites >= w.opts.ShardSites {
			if err := w.rotate(); err != nil {
				return err
			}
		}
	}
	if w.shardSites == 0 {
		w.shardChrom, w.shardFirst = chrom, pos
	}
	w.shardLast = pos
	w.shardSites++
	return nil
}

// Rotate closes the current shard and opens the next, on request.
//
// A CALLER THAT HOLDS STATE SPANNING SITES MUST DRIVE THIS ITSELF, before it
// extends that state into a site belonging to the next shard. The converter
// extends each sample's open run to position P and only then writes site P, so
// a rotation discovered while writing P emits a run that covers P into the
// shard that is closing -- and P belongs to the next one. Every sample then
// reads as never assayed at the first site of every shard.
//
// WouldRotate answers "is the next site a boundary" so the caller can close its
// runs and Rotate before touching the record at all.
func (w *Writer) Rotate() error { return w.rotate() }

// rotate closes the current shard and opens the next.
func (w *Writer) rotate() error {
	// The caller's state first, so anything spanning sites lands in the shard it
	// belongs to. The converter's open callable runs are the case this exists
	// for: a run must lie wholly inside one shard, or a query would have to read
	// every earlier shard in case a run started there and reached in.
	if w.opts.BeforeRotate != nil {
		if err := w.opts.BeforeRotate(); err != nil {
			return fmt.Errorf("flushing before the shard boundary at %s:%d: %w",
				w.shardChrom, w.shardLast, err)
		}
	}
	if err := w.closeShard(); err != nil {
		return err
	}
	w.shardIdx++
	return w.openShard()
}

// closeShard flushes and finalises the current shard, recording what it holds.
func (w *Writer) closeShard() error {
	for _, fn := range []func() error{w.flushCalls, w.flushSites, w.flushRegions} {
		if err := fn(); err != nil {
			return err
		}
	}
	// Exactly one of the two sites writers is live; see openShard.
	closers := []interface{ Close() error }{w.calls, w.regions}
	if w.sitesAny != nil {
		closers = append(closers, w.sitesAny)
	} else {
		closers = append(closers, w.sites)
	}
	for _, c := range closers {
		if err := c.Close(); err != nil {
			return err
		}
	}
	// The sinks for this shard, which are the last three opened. Marked so the
	// final closeFiles does not close them a second time -- which is an error,
	// and one the writer would report as a failed conversion.
	n := len(w.members)
	if n >= 3 {
		for i := n - 3; i < n; i++ {
			if w.members[i].closed {
				continue
			}
			if err := w.members[i].sw.Close(); err != nil {
				return err
			}
			w.members[i].closed = true
		}
	}
	if w.opts.ShardSites <= 0 {
		return nil
	}

	// A shard that received no sites is not recorded. It happens when the input
	// ends exactly on a boundary, and an entry claiming a zero-width range would
	// be a shard no query could ever match and every query would test.
	if w.shardSites == 0 {
		return nil
	}
	for _, member := range []string{CallsMember, SitesMember, RegionsMember} {
		w.shards[member] = append(w.shards[member], ShardInfo{
			Name:  ShardFile(member, w.shardIdx),
			Chrom: w.shardChrom,
			First: w.shardFirst,
			Last:  w.shardLast,
			Rows:  w.shardRows[member],
		})
	}
	return nil
}

// WouldRotate reports whether the next site at this chromosome starts a new
// shard, so a caller can flush state that spans sites BEFORE extending it into
// a site that no longer belongs to the open shard.
//
// THE ORDERING THIS EXISTS FOR, which cost a differential test to find. The
// converter extends each sample's open run to position P and only then writes
// site P. If the rotation is discovered while writing P, the run has already
// been extended to P and is emitted, ending at P, into the shard that is
// closing -- while P itself belongs to the next one. Every sample then reads as
// NEVER ASSAYED at the first site of every shard: 1,796 of 80,000 states on a
// 200-sample fixture, which is nine boundaries times two hundred people.
//
// BeforeRotate cannot fix that on its own, because by the time it fires the run
// already reaches into the new shard. The caller has to ask first.
func (w *Writer) WouldRotate(chrom string) bool {
	if w.opts.ShardSites <= 0 || w.shardSites == 0 {
		return false
	}
	return !SameChrom(chrom, w.shardChrom) || w.shardSites >= w.opts.ShardSites
}
