package tabix

import (
	"bytes"
	"encoding/binary"
	"math"
	"math/rand"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/htsio/bgzf"
)

// Record counts live in each reference's pseudo-bin (see RefStats). These tests
// check that cghts writes them, that they are right, that they agree field for
// field with what htslib writes for the same file, and that a pseudo-bin can
// never be read as a bin.

// wantCounts counts a fixture's records per contig, in contig order.
func wantCounts(fx tabixFixture) []uint64 {
	byName := map[string]uint64{}
	for _, l := range fx.lines {
		byName[l[:strings.IndexByte(l, '\t')]]++
	}
	var out []uint64
	for _, c := range fx.contigs {
		if n := byName[c.name]; n > 0 {
			out = append(out, n)
		}
	}
	return out
}

// lineAt reads the line starting at virtual offset vo, or "" at end of file.
func lineAt(t *testing.T, r *Reader, vo bgzf.VirtualOffset) string {
	t.Helper()
	if err := r.ir.SeekToVirtualOffset(vo); err != nil {
		t.Fatalf("seeking to %d:%d: %v", vo.BlockOffset(), vo.WithinBlock(), err)
	}
	line, _ := readLine(r.ir, nil)
	return string(line)
}

func checkCounts(t *testing.T, what string, r *Reader, want []uint64) {
	t.Helper()
	got, ok := r.RecordCounts()
	if !ok {
		t.Fatalf("%s: index records no counts", what)
	}
	if len(got) != len(want) {
		t.Fatalf("%s: %d references, want %d", what, len(got), len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("%s: %s has %d records, want %d", what, r.RefNames()[i], got[i], want[i])
		}
	}
	if n, ok := r.idx.(interface{ NoCoor() (uint64, bool) }).NoCoor(); !ok || n != 0 {
		t.Errorf("%s: n_no_coor = %d, %v; want 0, recorded", what, n, ok)
	}
}

// checkExtents checks each reference's [Begin, End) against the data: Begin
// is its first record, and End is where the next reference's first record
// starts, or the end of the file.
func checkExtents(t *testing.T, what string, r *Reader, fx tabixFixture) {
	t.Helper()
	first := map[string]string{}
	for _, l := range fx.lines {
		name := l[:strings.IndexByte(l, '\t')]
		if _, ok := first[name]; !ok {
			first[name] = l
		}
	}
	names := r.RefNames()
	for i, name := range names {
		st, ok := r.idx.RefStats(i)
		if !ok {
			t.Fatalf("%s: %s has no stats", what, name)
		}
		if got := lineAt(t, r, st.Begin); got != first[name] {
			t.Errorf("%s: %s begins at %q, want its first record %q", what, name, got, first[name])
		}
		next := ""
		if i+1 < len(names) {
			next = first[names[i+1]]
		}
		if got := lineAt(t, r, st.End); got != next {
			t.Errorf("%s: %s ends before %q, want %q", what, name, got, next)
		}
	}
}

func TestIndexRecordsCountsAndExtents(t *testing.T) {
	for _, fx := range tabixFixtures() {
		t.Run(fx.name, func(t *testing.T) {
			path := filepath.Join(t.TempDir(), "f.gz")
			writeFixture(t, path, fx)
			if err := NewIndexWriter(fx.opts).WriteIndex(path); err != nil {
				t.Fatal(err)
			}
			r, err := NewReader(path)
			if err != nil {
				t.Fatal(err)
			}
			defer r.Close()
			checkCounts(t, "cghts index", r, wantCounts(fx))
			checkExtents(t, "cghts index", r, fx)
		})
	}
}

// The sorting Writer builds its index on another path from IndexWriter.
func TestSortingWriterRecordsCounts(t *testing.T) {
	fx := tabixFixtures()[0]
	path := filepath.Join(t.TempDir(), "sorted.bed.gz")
	w := NewWriter(path, NewWriterOpts().BED().AutoIndex())
	shuffled := append([]string{}, fx.lines...)
	rand.New(rand.NewSource(7)).Shuffle(len(shuffled), func(i, j int) { shuffled[i], shuffled[j] = shuffled[j], shuffled[i] })
	for _, l := range shuffled {
		if err := w.Write(l); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	r, err := NewReader(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	checkCounts(t, "sorting writer", r, wantCounts(fx))
}

// The pseudo-bin cghts writes must be the one htslib writes: same counts, same
// extent, for every fixture layout. And htslib's own reader must see it --
// bcftools index -n is the one CLI that prints it, for VCF.
func TestRecordCountsAgreeWithHtslib(t *testing.T) {
	bin := htslibTabix(t)
	for _, fx := range tabixFixtures() {
		t.Run(fx.name, func(t *testing.T) {
			dir := t.TempDir()
			// Named for the format: bcftools finds a VCF's index only by its name.
			hts := filepath.Join(dir, "hts."+fx.preset+".gz")
			cg := filepath.Join(dir, "cghts."+fx.preset+".gz")
			writeFixture(t, hts, fx)
			writeFixture(t, cg, fx)
			if out, err := exec.Command(bin, "-p", fx.preset, hts).CombinedOutput(); err != nil {
				t.Fatalf("tabix: %v\n%s", err, out)
			}
			if err := NewIndexWriter(fx.opts).WriteIndex(cg); err != nil {
				t.Fatal(err)
			}
			onHts, err := NewReader(hts)
			if err != nil {
				t.Fatal(err)
			}
			defer onHts.Close()
			onCg, err := NewReader(cg)
			if err != nil {
				t.Fatal(err)
			}
			defer onCg.Close()

			checkCounts(t, "htslib index", onHts, wantCounts(fx))
			for i, name := range onHts.RefNames() {
				want, _ := onHts.idx.RefStats(i)
				got, ok := onCg.idx.RefStats(onCg.idx.RefID(name))
				if !ok || got != want {
					t.Errorf("%s: cghts pseudo-bin %+v, htslib's %+v", name, got, want)
				}
			}

			if fx.preset != "vcf" {
				return
			}
			bcftools, err := exec.LookPath("bcftools")
			if err != nil {
				return
			}
			out, err := exec.Command(bcftools, "index", "-n", cg).CombinedOutput()
			if err != nil {
				t.Fatalf("bcftools index -n: %v\n%s", err, out)
			}
			if got := strings.TrimSpace(string(out)); got != strconv.Itoa(len(fx.lines)) {
				t.Errorf("bcftools index -n on the cghts index: %s, want %d", got, len(fx.lines))
			}
		})
	}
}

// A CSI carries the pseudo-bin too, at a number set by its depth.
func TestCSIRecordCounts(t *testing.T) {
	bin := htslibTabix(t)
	for _, fx := range tabixFixtures() {
		if fx.layout != "plain" {
			continue
		}
		t.Run(fx.name, func(t *testing.T) {
			path := filepath.Join(t.TempDir(), "f.gz")
			writeFixture(t, path, fx)
			if out, err := exec.Command(bin, "-C", "-p", fx.preset, path).CombinedOutput(); err != nil {
				t.Fatalf("tabix -C: %v\n%s", err, out)
			}
			r, err := NewReader(path)
			if err != nil {
				t.Fatal(err)
			}
			defer r.Close()
			if _, isCSI := r.idx.(*CSIIndex); !isCSI {
				t.Fatalf("read %T, want a CSI", r.idx)
			}
			checkCounts(t, "htslib CSI", r, wantCounts(fx))
			checkExtents(t, "htslib CSI", r, fx)
		})
	}
}

// craftedRef is one reference of a hand-built TBI: its bins, written in the
// order given, pseudo-bin included if the caller lists one.
type craftedRef struct {
	bins []uint32
	chs  [][]Chunk
}

// craftTBI builds a minimal TBI (generic preset, one name per ref) so a test
// can say exactly what an index holds, including what cghts would never write.
func craftTBI(t *testing.T, refs []craftedRef, noCoor *uint64) *BinIndex {
	t.Helper()
	var buf bytes.Buffer
	w := bgzf.NewWriter(&buf)
	put := func(v any) {
		if err := binary.Write(w, binary.LittleEndian, v); err != nil {
			t.Fatal(err)
		}
	}
	w.Write([]byte("TBI\x01"))
	put(int32(len(refs)))
	put([6]int32{0, 1, 2, 3, '#', 0})
	var names []byte
	for i := range refs {
		names = append(names, []byte("chr"+strconv.Itoa(i+1)+"\x00")...)
	}
	put(int32(len(names)))
	w.Write(names)
	for _, ref := range refs {
		put(int32(len(ref.bins)))
		for i, b := range ref.bins {
			put(b)
			put(int32(len(ref.chs[i])))
			for _, c := range ref.chs[i] {
				put(uint64(c.Begin))
				put(uint64(c.End))
			}
		}
		put(int32(0)) // no linear index
	}
	if noCoor != nil {
		put(*noCoor)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	idx, err := LoadTBIFrom(&buf)
	if err != nil {
		t.Fatal(err)
	}
	return idx
}

// An index with no pseudo-bin -- cghts before v0.12.9 -- has unknown counts,
// never zero ones; a reference with no bins at all is known to be empty.
func TestIndexWithoutCountsIsUnknown(t *testing.T) {
	real := []Chunk{{Begin: 1 << 16, End: 2 << 16}}
	idx := craftTBI(t, []craftedRef{
		{bins: []uint32{4681}, chs: [][]Chunk{real}},
		{}, // no records
	}, nil)
	if st, ok := idx.RefStats(0); ok {
		t.Errorf("a reference without a pseudo-bin: %+v, known; want unknown", st)
	}
	if st, ok := idx.RefStats(1); !ok || st != (RefStats{}) {
		t.Errorf("a reference without bins: %+v, %v; want zero, known", st, ok)
	}
	if _, ok := idx.NoCoor(); ok {
		t.Error("n_no_coor reported for an index that ends without one")
	}
	r := &Reader{idx: idx}
	if counts, ok := r.RecordCounts(); ok {
		t.Errorf("RecordCounts = %v, complete; want incomplete", counts)
	}
}

// Whole-contig queries pass math.MaxInt32 as the end. Unclamped, the bin loops
// ran on past 2^29 into other levels' bins and the pseudo-bin, whose second
// chunk is (mapped, unmapped) counts: a BAI with unmapped reads placed on a
// contig would have sent the query to offsets that are counts.
func TestWholeContigQueryNeverReadsThePseudoBin(t *testing.T) {
	for _, b := range reg2bins(0, math.MaxInt32) {
		if b >= tbiPseudoBin {
			t.Fatalf("reg2bins(0, MaxInt32) includes bin %d", b)
		}
	}
	if a, b := reg2bins(0, math.MaxInt32), reg2bins(0, tbiMaxPos); len(a) != len(b) {
		t.Errorf("reg2bins(0, MaxInt32) gives %d bins, reg2bins(0, 2^29) %d", len(a), len(b))
	}
	pseudo, _ := csiPseudoBin(5)
	for _, b := range csiReg2bins(0, math.MaxInt32, 14, 5) {
		if b >= pseudo {
			t.Fatalf("csiReg2bins(0, MaxInt32) includes bin %d", b)
		}
	}

	real := []Chunk{{Begin: 1 << 16, End: 2 << 16}}
	stats := []Chunk{real[0], {Begin: 5, End: 1 << 40}} // 5 mapped, 2^40 unmapped
	idx := craftTBI(t, []craftedRef{{
		bins: []uint32{4681, tbiPseudoBin},
		chs:  [][]Chunk{real, stats},
	}}, nil)
	got := idx.Query(0, 0, math.MaxInt32)
	if len(got) != 1 || got[0] != real[0] {
		t.Errorf("whole-contig query = %v, want only the real chunk %v", got, real)
	}
	if st, ok := idx.RefStats(0); !ok || st.Mapped != 5 || st.Unmapped != 1<<40 {
		t.Errorf("RefStats = %+v, %v", st, ok)
	}
}
