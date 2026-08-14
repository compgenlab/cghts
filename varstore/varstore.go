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

// A VARSTORE is a collection of VOLUMES read as one.
//
// The three words this package uses, since two of them used to be one word:
//
//	varstore   the whole archive -- what a caller holds, names and queries
//	volume     one chromosome's worth of it, a complete store on its own
//	shard      a coordinate range within a volume, x3 tables
//
// Volume and shard split by COORDINATE. The three parquet tables -- calls,
// sites, regions -- split by DATA KIND, which is a different axis: every shard
// has all three. They are tables, not a rung on the ladder.
//
// WHY THE ARCHIVE EXISTS. A callset almost never arrives as a single file:
// whole-genome data ships one VCF per chromosome, and converting it gives one
// volume per chromosome. Twenty-five volumes is the wrong unit to manage, to
// register, to grant access to, or to reason about -- a person has one genome,
// and the system should let somebody hold one object that means it.
//
// AND IT IS WHAT MAKES SMALL SHARDS AFFORDABLE, which is the less obvious half.
// A shard index lives in a manifest read on every open, so its size is a budget:
//
//	chr22 at 500 sites a shard    2,000 shards x3 tables    ~100 KB gzipped
//	whole genome at 500            52,000 shards x3         ~2.6 MB gzipped
//
// A single undivided archive therefore forces a choice between coarse shards
// and an unreadable manifest -- and coarse shards are exactly what makes a locus
// query slow. Split by chromosome, each index covers one chromosome and a query
// opens only the chromosomes it touches, so the choice disappears.
//
// THE VOLUMES ARE DISJOINT BY CHROMOSOME. That is what separates an archive
// from the composite PARTS a caller may build above it: parts overlap and need
// a conflict policy, because they are independent observations of the same
// people. Volumes observe different regions of the genome, so a locus is
// answered by exactly one of them and there is nothing to resolve. Modelling
// per-chromosome volumes as parts would mean asking every volume about every
// locus so that all but one could answer "not mine".

// StoreManifestFile is the marker that makes a directory an archive.
const StoreManifestFile = "varstore.json.gz"

// StoreManifestVersion is the format version this package writes.
const StoreManifestVersion = 1

// StoreManifest describes an archive and the agreement its volumes are held to.
type StoreManifest struct {
	FormatVersion int    `json:"format_version"`
	Complete      bool   `json:"complete"`
	Created       string `json:"created"`
	Program       string `json:"program,omitempty"`
	Command       string `json:"command,omitempty"`

	// Meta is what the caller said the archive IS, as for a volume.
	Meta map[string]string `json:"meta,omitempty"`

	// Samples is the roster every volume must hold, in order.
	//
	// RECORDED ONCE, HERE, and this is the authority rather than a summary. A
	// store whose volumes disagreed about who is in it would answer a query with
	// a different population per chromosome, and no caller could see that from
	// the outside.
	Samples []string `json:"samples"`

	// Params are the conversion settings every volume must share. A volume
	// converted under a different depth gate does not mean the same thing by
	// "callable", and a query spanning both would be comparing two claims.
	Params ManifestParams `json:"params"`

	Volumes []VolumeInfo `json:"volumes"`
}

// VolumeInfo is one volume of an archive.
type VolumeInfo struct {
	// Name locates the volume. A bare name is relative to the store directory,
	// which keeps an archive one thing to copy, move or delete; a locator with a
	// scheme is used as given, so volumes may live apart when they must.
	Name string `json:"name"`

	// Chroms are the chromosomes this volume holds, and they must not overlap
	// another volume's -- a locus has to have exactly one volume that can
	// answer it.
	Chroms []string `json:"chroms"`

	Sites int64 `json:"sites"`
	Calls int64 `json:"calls"`
}

// VarStore reads a collection of volumes as one.
//
// Volumes open LAZILY and at most once: a chr22 query opens one volume, not
// twenty-five. Each is verified against the store's declared agreement as it
// opens, so drift is caught when a volume is first read rather than requiring
// every manifest to be fetched up front -- which would undo the laziness.
type VarStore struct {
	base string
	man  *StoreManifest

	// byChrom maps a canonical chromosome to the volume index holding it.
	byChrom map[string]int

	mu     sync.Mutex
	opened []*ParquetVolume
	once   []*sync.Once
	errs   []error
}

// OpenStore opens an archive, reading only its manifest.
func OpenStore(ctx context.Context, base string) (*VarStore, error) {
	base = TrimStoreSuffix(base)
	man, err := readStoreManifest(ctx, base)
	if err != nil {
		return nil, err
	}
	if !man.Complete {
		return nil, fmt.Errorf("%s: the store manifest is not marked complete; the store was not finished", base)
	}

	s := &VarStore{
		base: base, man: man,
		byChrom: map[string]int{},
		opened:  make([]*ParquetVolume, len(man.Volumes)),
		once:    make([]*sync.Once, len(man.Volumes)),
		errs:    make([]error, len(man.Volumes)),
	}
	for i := range man.Volumes {
		s.once[i] = &sync.Once{}
		for _, c := range man.Volumes[i].Chroms {
			k := CanonKey(c)
			if prev, dup := s.byChrom[k]; dup {
				return nil, fmt.Errorf(
					"%s: %s is held by both %s and %s; a store's volumes must not overlap, "+
						"or a locus has two answers and no rule for choosing",
					base, c, man.Volumes[prev].Name, man.Volumes[i].Name)
			}
			s.byChrom[k] = i
		}
	}
	return s, nil
}

// IsStore reports whether a locator names an archive rather than a bare volume.
func IsStore(ctx context.Context, base string) bool {
	_, err := readStoreManifest(ctx, TrimStoreSuffix(base))
	return err == nil
}

// Samples is the roster, from the store manifest.
func (s *VarStore) Samples() ([]string, error) { return s.man.Samples, nil }

// StoreManifest returns the archive's own manifest.
func (s *VarStore) StoreManifest() *StoreManifest { return s.man }

// Close releases every volume that was opened.
func (s *VarStore) Close() error {
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
		return fmt.Errorf("closing store %s: %s", s.base, strings.Join(errs, "; "))
	}
	return nil
}

// volume opens the store holding a chromosome, or reports that none does.
func (s *VarStore) volume(chrom string) (*ParquetVolume, bool, error) {
	i, ok := s.byChrom[CanonKey(chrom)]
	if !ok {
		// Not an error: an archive covering the autosomes asked about chrM answers
		// "nothing here", which is NotAssayed and not a failure.
		return nil, false, nil
	}
	st, err := s.openVolumeAt(i)
	if err != nil {
		return nil, false, err
	}
	return st, true, nil
}

func (s *VarStore) openVolumeAt(i int) (*ParquetVolume, error) {
	s.once[i].Do(func() {
		m := s.man.Volumes[i]
		locator := m.Name
		if !strings.Contains(locator, "://") && !strings.HasPrefix(locator, "/") {
			locator = joinStore(s.base, locator)
		}
		st, err := OpenParquetContext(context.Background(), locator)
		if err != nil {
			s.errs[i] = fmt.Errorf("opening volume %s: %w", m.Name, err)
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

// agrees holds a volume to what the store declares.
//
// CHECKED AS THE MEMBER OPENS, not at OpenStore, so a chr22 query still opens one
// volume. The store manifest is the authority and each volume is measured against
// it -- which catches a volume swapped in from another conversion, the case that
// would otherwise answer with a different population for one chromosome and
// nothing to show for it.
func (s *VarStore) agrees(m VolumeInfo, st *ParquetVolume) error {
	got, err := st.Samples()
	if err != nil {
		return err
	}
	if len(got) != len(s.man.Samples) {
		return fmt.Errorf(
			"volume %s holds %d samples but the store declares %d; a query spanning it "+
				"would answer with a different population for those chromosomes",
			m.Name, len(got), len(s.man.Samples))
	}
	for i := range got {
		if got[i] != s.man.Samples[i] {
			return fmt.Errorf(
				"volume %s has %q where the store declares %q at position %d; "+
					"genotype columns are positional and a mismatch attributes calls to the wrong person",
				m.Name, got[i], s.man.Samples[i], i)
		}
	}
	if mm := st.VolumeManifest(); mm != nil {
		if mm.Params.MinDP != s.man.Params.MinDP {
			return fmt.Errorf(
				"volume %s was converted at --min-dp %d but the store declares %d; "+
					"they do not mean the same thing by callable",
				m.Name, mm.Params.MinDP, s.man.Params.MinDP)
		}
		if !sameBands(mm.Params.DepthBands, s.man.Params.DepthBands) {
			return fmt.Errorf(
				"volume %s was banded %v but the store declares %v; a recorded depth "+
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

// Chroms lists every chromosome the store can answer about, in volume order.
func (s *VarStore) Chroms() []string {
	var out []string
	for _, m := range s.man.Volumes {
		out = append(out, m.Chroms...)
	}
	sort.Strings(out)
	return out
}

// StoreManifestPath is the store's marker file.
func StoreManifestPath(base string) string { return joinStore(base, StoreManifestFile) }

// readStoreManifest reads a store's manifest, and fails cleanly when the locator is
// not an archive.
func readStoreManifest(ctx context.Context, base string) (*StoreManifest, error) {
	locator := StoreManifestPath(base)
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
	var m StoreManifest
	if err := json.NewDecoder(zr).Decode(&m); err != nil {
		return nil, fmt.Errorf("reading %s: %w", locator, err)
	}
	return &m, nil
}

// WriteStoreManifest writes a store's manifest through a sink.
//
// LAST, like a store's, and for the same reason: it is the completion marker.
// An archive whose manifest exists is one whose volumes were all written, and
// that appeared first would vouch for volumes that may not be there.
func WriteStoreManifest(base string, m StoreManifest) error {
	sink, err := OpenSink(base)
	if err != nil {
		return err
	}
	f, err := sink.Create(StoreManifestFile)
	if err != nil {
		return fmt.Errorf("creating the store manifest in %s: %w", sink.Describe(), err)
	}
	zw := gzip.NewWriter(f)
	enc := json.NewEncoder(zw)
	enc.SetIndent("", "  ")
	if err := enc.Encode(m); err != nil {
		f.Close()
		return fmt.Errorf("encoding store manifest: %w", err)
	}
	if err := zw.Close(); err != nil {
		f.Close()
		return fmt.Errorf("compressing store manifest: %w", err)
	}
	return f.Close()
}

// BuildStore describes an archive over volumes that already exist, verifying that they
// agree before writing anything.
//
// SEPARATE FROM CONVERSION on purpose. A whole-genome callset arrives as one
// VCF per chromosome and is converted one store at a time -- often on different
// machines, often days apart -- so the store is assembled from finished stores
// rather than produced by a single run. It also means an archive can be rebuilt, or
// a volume replaced, without reconverting anything.
func BuildStore(ctx context.Context, base string, volumes []string, program, command string) (*StoreManifest, error) {
	if len(volumes) == 0 {
		return nil, fmt.Errorf("an archive needs at least one volume")
	}
	man := StoreManifest{
		FormatVersion: StoreManifestVersion,
		Complete:      true,
		Created:       time.Now().UTC().Format(time.RFC3339),
		Program:       program,
		Command:       command,
	}
	seen := map[string]string{}
	for i, locator := range volumes {
		full := locator
		if !strings.Contains(full, "://") && !strings.HasPrefix(full, "/") {
			full = joinStore(base, full)
		}
		st, err := OpenParquetContext(ctx, full)
		if err != nil {
			return nil, fmt.Errorf("reading prospective volume %s: %w", locator, err)
		}
		samples, err := st.Samples()
		if err != nil {
			st.Close()
			return nil, err
		}
		sm := st.VolumeManifest()

		if i == 0 {
			man.Samples = append([]string{}, samples...)
			if sm != nil {
				man.Params = sm.Params
			}
		} else {
			// THE AGREEMENT IS ESTABLISHED HERE, once, rather than rediscovered
			// on every query. An archive that could not be built is better than one
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
					"%s and %s both hold %s; a store's volumes must not overlap, or a locus "+
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
		man.Volumes = append(man.Volumes, VolumeInfo{
			Name: locator, Chroms: chroms, Sites: sites, Calls: calls,
		})
		st.Close()
	}
	return &man, nil
}

func agreeWith(man StoreManifest, locator string, samples []string, sm *VolumeManifest) error {
	if len(samples) != len(man.Samples) {
		return fmt.Errorf(
			"%s holds %d samples but the store's first volume holds %d; an archive answers one "+
				"population, not a different one per chromosome",
			locator, len(samples), len(man.Samples))
	}
	for i := range samples {
		if samples[i] != man.Samples[i] {
			return fmt.Errorf(
				"%s has %q where the first volume has %q at position %d; genotype columns "+
					"are positional and a mismatch attributes calls to the wrong person",
				locator, samples[i], man.Samples[i], i)
		}
	}
	if sm == nil {
		return nil
	}
	if sm.Params.MinDP != man.Params.MinDP {
		return fmt.Errorf(
			"%s was converted at --min-dp %d but the store is %d; they do not mean the same "+
				"thing by callable", locator, sm.Params.MinDP, man.Params.MinDP)
	}
	if !sameBands(sm.Params.DepthBands, man.Params.DepthBands) {
		return fmt.Errorf(
			"%s was banded %v but the store is %v; a recorded depth means a different span in each",
			locator, sm.Params.DepthBands, man.Params.DepthBands)
	}
	return nil
}
