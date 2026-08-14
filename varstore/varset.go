package varstore

import (
	"compress/gzip"
	"context"
	"encoding/json"
	"fmt"
	"io"
	"sort"
	"strings"
	"sync"
	"time"

	"github.com/compgenlab/cghts/iosource"
)

// A VARSET is a set of stores read as one.
//
// WHY IT EXISTS. A callset almost never arrives as a single file: whole-genome
// data ships one VCF per chromosome, and converting it gives one store per
// chromosome. Twenty-five stores is the wrong unit to manage, to register, to
// grant access to, or to reason about -- a person has one genome, and the
// system should let somebody hold one object that means it.
//
// AND IT IS WHAT MAKES SMALL SHARDS AFFORDABLE, which is the less obvious half.
// A shard index lives in a manifest read on every open, so its size is a budget:
//
//	chr22 at 500 sites a shard    2,000 shards x3 members    ~100 KB gzipped
//	whole genome at 500            52,000 shards x3          ~2.6 MB gzipped
//
// A single whole-genome store therefore forces a choice between coarse shards
// and an unreadable manifest -- and coarse shards are exactly what makes a locus
// query slow. Split by chromosome, each index covers one chromosome and a query
// opens only the chromosomes it touches, so the choice disappears.
//
// THE MEMBERS ARE DISJOINT BY CHROMOSOME. That is what separates a set from the
// composite PARTS a caller may build above it: parts overlap and need a
// conflict policy, because they are independent observations of the same
// people. Set members observe different regions of the genome, so a locus is
// answered by exactly one of them and there is nothing to resolve. Modelling
// per-chromosome stores as parts would mean asking every member about every
// locus so that all but one could answer "not mine".

// SetManifestFile is the marker that makes a directory a set.
const SetManifestFile = "varset.json.gz"

// SetManifestVersion is the format version this package writes.
const SetManifestVersion = 1

// SetManifest describes a set and the agreement its members are held to.
type SetManifest struct {
	FormatVersion int    `json:"format_version"`
	Complete      bool   `json:"complete"`
	Created       string `json:"created"`
	Program       string `json:"program,omitempty"`
	Command       string `json:"command,omitempty"`

	// Meta is what the caller said the set IS, as for a store.
	Meta map[string]string `json:"meta,omitempty"`

	// Samples is the roster every member must hold, in order.
	//
	// RECORDED ONCE, HERE, and this is the authority rather than a summary. A
	// set whose members disagreed about who is in it would answer a query with
	// a different population per chromosome, and no caller could see that from
	// the outside.
	Samples []string `json:"samples"`

	// Params are the conversion settings every member must share. A member
	// converted under a different depth gate does not mean the same thing by
	// "callable", and a query spanning both would be comparing two claims.
	Params ManifestParams `json:"params"`

	Members []SetMember `json:"members"`
}

// SetMember is one store in a set.
type SetMember struct {
	// Name locates the member. A bare name is relative to the set directory,
	// which keeps a set one thing to copy, move or delete; a locator with a
	// scheme is used as given, so members may live apart when they must.
	Name string `json:"name"`

	// Chroms are the chromosomes this member holds, and they must not overlap
	// another member's -- a locus has to have exactly one member that can
	// answer it.
	Chroms []string `json:"chroms"`

	Sites int64 `json:"sites"`
	Calls int64 `json:"calls"`
}

// VarSet reads a set of stores as one.
//
// Members open LAZILY and at most once: a chr22 query opens one member, not
// twenty-five. Each is verified against the set's declared agreement as it
// opens, so drift is caught when a member is first read rather than requiring
// every manifest to be fetched up front -- which would undo the laziness.
type VarSet struct {
	base string
	man  *SetManifest

	// byChrom maps a canonical chromosome to the member index holding it.
	byChrom map[string]int

	mu     sync.Mutex
	opened []*ParquetStore
	once   []*sync.Once
	errs   []error
}

// OpenSet opens a set, reading only its manifest.
func OpenSet(ctx context.Context, base string) (*VarSet, error) {
	base = TrimStoreSuffix(base)
	man, err := readSetManifest(ctx, base)
	if err != nil {
		return nil, err
	}
	if !man.Complete {
		return nil, fmt.Errorf("%s: the set manifest is not marked complete; the set was not finished", base)
	}

	s := &VarSet{
		base: base, man: man,
		byChrom: map[string]int{},
		opened:  make([]*ParquetStore, len(man.Members)),
		once:    make([]*sync.Once, len(man.Members)),
		errs:    make([]error, len(man.Members)),
	}
	for i := range man.Members {
		s.once[i] = &sync.Once{}
		for _, c := range man.Members[i].Chroms {
			k := CanonKey(c)
			if prev, dup := s.byChrom[k]; dup {
				return nil, fmt.Errorf(
					"%s: %s is held by both %s and %s; a set's members must not overlap, "+
						"or a locus has two answers and no rule for choosing",
					base, c, man.Members[prev].Name, man.Members[i].Name)
			}
			s.byChrom[k] = i
		}
	}
	return s, nil
}

// IsSet reports whether a locator names a set rather than a store.
func IsSet(ctx context.Context, base string) bool {
	_, err := readSetManifest(ctx, TrimStoreSuffix(base))
	return err == nil
}

// Samples is the roster, from the set manifest.
func (s *VarSet) Samples() ([]string, error) { return s.man.Samples, nil }

// Manifest returns the set's own manifest.
func (s *VarSet) SetManifest() *SetManifest { return s.man }

// Close releases every member that was opened.
func (s *VarSet) Close() error {
	s.mu.Lock()
	defer s.mu.Unlock()
	var errs []string
	for _, m := range s.opened {
		if m == nil {
			continue
		}
		if err := m.Close(); err != nil {
			errs = append(errs, err.Error())
		}
	}
	if len(errs) > 0 {
		return fmt.Errorf("closing set %s: %s", s.base, strings.Join(errs, "; "))
	}
	return nil
}

// member opens the store holding a chromosome, or reports that none does.
func (s *VarSet) member(chrom string) (*ParquetStore, bool, error) {
	i, ok := s.byChrom[CanonKey(chrom)]
	if !ok {
		// Not an error: a set covering the autosomes asked about chrM answers
		// "nothing here", which is NotAssayed and not a failure.
		return nil, false, nil
	}
	st, err := s.openMemberAt(i)
	if err != nil {
		return nil, false, err
	}
	return st, true, nil
}

func (s *VarSet) openMemberAt(i int) (*ParquetStore, error) {
	s.once[i].Do(func() {
		m := s.man.Members[i]
		locator := m.Name
		if !strings.Contains(locator, "://") && !strings.HasPrefix(locator, "/") {
			locator = joinStore(s.base, locator)
		}
		st, err := OpenParquetContext(context.Background(), locator)
		if err != nil {
			s.errs[i] = fmt.Errorf("opening set member %s: %w", m.Name, err)
			return
		}
		if err := s.agrees(m, st); err != nil {
			st.Close()
			s.errs[i] = err
			return
		}
		s.mu.Lock()
		s.opened[i] = st
		s.mu.Unlock()
	})
	return s.opened[i], s.errs[i]
}

// agrees holds a member to what the set declares.
//
// CHECKED AS THE MEMBER OPENS, not at OpenSet, so a chr22 query still opens one
// member. The set manifest is the authority and each member is measured against
// it -- which catches a member swapped in from another conversion, the case that
// would otherwise answer with a different population for one chromosome and
// nothing to show for it.
func (s *VarSet) agrees(m SetMember, st *ParquetStore) error {
	got, err := st.Samples()
	if err != nil {
		return err
	}
	if len(got) != len(s.man.Samples) {
		return fmt.Errorf(
			"set member %s holds %d samples but the set declares %d; a query spanning it "+
				"would answer with a different population for those chromosomes",
			m.Name, len(got), len(s.man.Samples))
	}
	for i := range got {
		if got[i] != s.man.Samples[i] {
			return fmt.Errorf(
				"set member %s has %q where the set declares %q at position %d; "+
					"genotype columns are positional and a mismatch attributes calls to the wrong person",
				m.Name, got[i], s.man.Samples[i], i)
		}
	}
	if mm := st.Manifest(); mm != nil {
		if mm.Params.MinDP != s.man.Params.MinDP {
			return fmt.Errorf(
				"set member %s was converted at --min-dp %d but the set declares %d; "+
					"they do not mean the same thing by callable",
				m.Name, mm.Params.MinDP, s.man.Params.MinDP)
		}
		if !sameBands(mm.Params.DepthBands, s.man.Params.DepthBands) {
			return fmt.Errorf(
				"set member %s was banded %v but the set declares %v; a recorded depth "+
					"means a different span in each",
				m.Name, mm.Params.DepthBands, s.man.Params.DepthBands)
		}
	}
	return nil
}

func sameBands(a, b []int32) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// Chroms lists every chromosome the set can answer about, in member order.
func (s *VarSet) Chroms() []string {
	var out []string
	for _, m := range s.man.Members {
		out = append(out, m.Chroms...)
	}
	sort.Strings(out)
	return out
}

// SetManifestPath is the set's marker file.
func SetManifestPath(base string) string { return joinStore(base, SetManifestFile) }

// readSetManifest reads a set's manifest, and fails cleanly when the locator is
// not a set.
func readSetManifest(ctx context.Context, base string) (*SetManifest, error) {
	locator := SetManifestPath(base)
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
	var m SetManifest
	if err := json.NewDecoder(zr).Decode(&m); err != nil {
		return nil, fmt.Errorf("reading %s: %w", locator, err)
	}
	return &m, nil
}

// WriteSetManifest writes a set's manifest through a sink.
//
// LAST, like a store's, and for the same reason: it is the completion marker.
// A set whose manifest exists is a set whose members were all written, and one
// that appeared first would vouch for members that may not be there.
func WriteSetManifest(base string, m SetManifest) error {
	sink, err := OpenSink(base)
	if err != nil {
		return err
	}
	f, err := sink.Create(SetManifestFile)
	if err != nil {
		return fmt.Errorf("creating the set manifest in %s: %w", sink.Describe(), err)
	}
	zw := gzip.NewWriter(f)
	enc := json.NewEncoder(zw)
	enc.SetIndent("", "  ")
	if err := enc.Encode(m); err != nil {
		f.Close()
		return fmt.Errorf("encoding set manifest: %w", err)
	}
	if err := zw.Close(); err != nil {
		f.Close()
		return fmt.Errorf("compressing set manifest: %w", err)
	}
	return f.Close()
}

// BuildSet describes a set over stores that already exist, verifying that they
// agree before writing anything.
//
// SEPARATE FROM CONVERSION on purpose. A whole-genome callset arrives as one
// VCF per chromosome and is converted one store at a time -- often on different
// machines, often days apart -- so the set is assembled from finished stores
// rather than produced by a single run. It also means a set can be rebuilt, or
// a member replaced, without reconverting anything.
func BuildSet(ctx context.Context, base string, members []string, program, command string) (*SetManifest, error) {
	if len(members) == 0 {
		return nil, fmt.Errorf("a set needs at least one member")
	}
	man := SetManifest{
		FormatVersion: SetManifestVersion,
		Complete:      true,
		Created:       time.Now().UTC().Format(time.RFC3339),
		Program:       program,
		Command:       command,
	}
	seen := map[string]string{}
	for i, locator := range members {
		full := locator
		if !strings.Contains(full, "://") && !strings.HasPrefix(full, "/") {
			full = joinStore(base, full)
		}
		st, err := OpenParquetContext(ctx, full)
		if err != nil {
			return nil, fmt.Errorf("reading prospective member %s: %w", locator, err)
		}
		samples, err := st.Samples()
		if err != nil {
			st.Close()
			return nil, err
		}
		sm := st.Manifest()

		if i == 0 {
			man.Samples = append([]string{}, samples...)
			if sm != nil {
				man.Params = sm.Params
			}
		} else {
			// THE AGREEMENT IS ESTABLISHED HERE, once, rather than rediscovered
			// on every query. A set that could not be built is better than one
			// that answers differently per chromosome.
			if err := agreeWith(man, locator, samples, sm); err != nil {
				st.Close()
				return nil, err
			}
		}

		var chroms []string
		for _, c := range sm.Chromosomes {
			if prev, dup := seen[CanonKey(c.Name)]; dup {
				st.Close()
				return nil, fmt.Errorf(
					"%s and %s both hold %s; a set's members must not overlap, or a locus "+
						"has two answers and no rule for choosing", prev, locator, c.Name)
			}
			seen[CanonKey(c.Name)] = locator
			chroms = append(chroms, c.Name)
		}
		sort.Strings(chroms)

		var sites, calls int64
		if sm != nil {
			sites, calls = sm.Counts.Sites, sm.Counts.Calls
		}
		man.Members = append(man.Members, SetMember{
			Name: locator, Chroms: chroms, Sites: sites, Calls: calls,
		})
		st.Close()
	}
	return &man, nil
}

func agreeWith(man SetManifest, locator string, samples []string, sm *Manifest) error {
	if len(samples) != len(man.Samples) {
		return fmt.Errorf(
			"%s holds %d samples but the set's first member holds %d; a set answers one "+
				"population, not a different one per chromosome",
			locator, len(samples), len(man.Samples))
	}
	for i := range samples {
		if samples[i] != man.Samples[i] {
			return fmt.Errorf(
				"%s has %q where the first member has %q at position %d; genotype columns "+
					"are positional and a mismatch attributes calls to the wrong person",
				locator, samples[i], man.Samples[i], i)
		}
	}
	if sm == nil {
		return nil
	}
	if sm.Params.MinDP != man.Params.MinDP {
		return fmt.Errorf(
			"%s was converted at --min-dp %d but the set is %d; they do not mean the same "+
				"thing by callable", locator, sm.Params.MinDP, man.Params.MinDP)
	}
	if !sameBands(sm.Params.DepthBands, man.Params.DepthBands) {
		return fmt.Errorf(
			"%s was banded %v but the set is %v; a recorded depth means a different span in each",
			locator, sm.Params.DepthBands, man.Params.DepthBands)
	}
	return nil
}
