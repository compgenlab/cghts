package varstore

import (
	"context"
	"fmt"
	"strings"
)

// A member split by coordinate.
//
// WHY SPLIT AT ALL, when Parquet already prunes row groups by their statistics.
// Because the statistics live in the FOOTER, at the end of the file: a locus
// query against a whole-genome member must fetch and parse a footer describing
// hundreds of gigabytes before it can decide to read none of it. Over object
// storage that is a large GET and a large parse to answer "no". A shard index
// in the manifest answers the same question having read only the manifest.
//
// Row-group pruning still happens INSIDE the shard that survives. The two are
// the same idea at two scales, and shard pruning is the one that scales with
// the genome rather than with the file.

// shardSet is how a member is read: one file, or many with an index.
//
// An unsplit store is a set of one whose range admits everything, so nothing
// downstream branches on which it is -- the reader treats every member as split
// and single files simply never get skipped.
type shardSet struct {
	name   string
	shards []*shard
}

type shard struct {
	info   ShardInfo
	m      *member
	bounds bool // false when the shard has no declared range and must be read
}

// openShardSet opens a member, split or not.
//
// The manifest decides. A member with no shard index is opened under its own
// name, which is every store written before splitting existed -- and its single
// shard declares no bounds, so no filter can skip it.
func openShardSet(ctx context.Context, base, member string, info MemberInfo) (*shardSet, error) {
	set := &shardSet{name: member}
	if len(info.Shards) == 0 {
		m, err := openMember(ctx, MemberPath(base, member))
		if err != nil {
			return nil, err
		}
		set.shards = []*shard{{m: m}}
		return set, nil
	}
	for _, si := range info.Shards {
		m, err := openMember(ctx, joinStore(base, si.Name))
		if err != nil {
			set.Close()
			return nil, fmt.Errorf("opening shard %s: %w", si.Name, err)
		}
		set.shards = append(set.shards, &shard{info: si, m: m, bounds: true})
	}
	return set, nil
}

func (s *shardSet) Close() error {
	if s == nil {
		return nil
	}
	var errs []string
	for _, sh := range s.shards {
		if err := sh.m.Close(); err != nil {
			errs = append(errs, err.Error())
		}
	}
	if len(errs) > 0 {
		return fmt.Errorf("closing %s: %s", s.name, strings.Join(errs, "; "))
	}
	return nil
}

// rows totals the set, for the manifest check.
func (s *shardSet) rows() int64 {
	var n int64
	for _, sh := range s.shards {
		if sh.bounds {
			n += sh.info.Rows
		}
	}
	return n
}

// single returns the lone member of an unsplit set, and nil otherwise. For the
// few readers that genuinely want one file -- the key/value metadata on calls,
// which a split member carries on its first shard.
func (s *shardSet) single() *member {
	if s == nil || len(s.shards) == 0 {
		return nil
	}
	return s.shards[0].m
}

// openOptionalShardSet opens a member that a store need not have.
//
// Absence is legitimate only where the manifest recorded nothing in it, which
// verifyAgainstManifest enforces -- so this reports presence rather than
// deciding what it means.
func openOptionalShardSet(ctx context.Context, base, member string, info MemberInfo) (*shardSet, bool, error) {
	if len(info.Shards) > 0 {
		set, err := openShardSet(ctx, base, member, info)
		return set, err == nil, err
	}
	m, present, err := openOptionalMember(ctx, MemberPath(base, member))
	if err != nil || !present {
		return &shardSet{name: member}, present, err
	}
	return &shardSet{name: member, shards: []*shard{{m: m}}}, true, nil
}

// split reports whether this member is stored as several coordinate shards.
func (s *shardSet) split() bool {
	return s != nil && len(s.shards) > 0 && s.shards[0].bounds
}
