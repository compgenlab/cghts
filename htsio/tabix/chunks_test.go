package tabix

import (
	"fmt"
	"math/rand"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

// longFeatureBED is coordinate-sorted BED in the shape that matters for a GTF:
// short features everywhere, and every so often one tens of megabases long,
// which the binning scheme can only put in a top-level bin. Those bins then
// hold records from all along the chromosome, interleaved with everyone
// else's -- the layout that exposed both bugs below. Names are random so the
// file does not compress to nothing.
func longFeatureBED(t *testing.T) (path string, lines []string) {
	t.Helper()
	rng := rand.New(rand.NewSource(7))
	for i := 0; i < 20000; i++ {
		start := i * 10000
		lines = append(lines, fmt.Sprintf("chr1\t%d\t%d\tshort_%d_%x", start, start+100, i, rng.Uint64()))
		if i%500 == 0 {
			lines = append(lines, fmt.Sprintf("chr1\t%d\t%d\tlong_%d_%x", start, start+70_000_000, i, rng.Uint64()))
		}
	}
	sort.SliceStable(lines, func(i, j int) bool { return bedStart(lines[i]) < bedStart(lines[j]) })
	path = filepath.Join(t.TempDir(), "long.bed.gz")
	writeBGZF(t, path, lines...)
	if err := NewIndexWriter(NewWriterOpts().BED()).WriteIndex(path); err != nil {
		t.Fatal(err)
	}
	return path, lines
}

func bedStart(line string) int {
	var s int
	fmt.Sscanf(strings.Split(line, "\t")[1], "%d", &s)
	return s
}

func bedEnd(line string) int {
	var e int
	fmt.Sscanf(strings.Split(line, "\t")[2], "%d", &e)
	return e
}

// Queries return exactly the overlapping records, once each. Reads now follow
// the index chunk by chunk, jumping the gaps between chunks, where they used to
// read on from the first chunk until the region was passed; the jumps must not
// lose or repeat a line.
func TestQueryAcrossChunkGapsMatchesAScan(t *testing.T) {
	path, lines := longFeatureBED(t)
	r, err := NewReader(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()

	starts := make([]int, len(lines))
	ends := make([]int, len(lines))
	for i, l := range lines {
		starts[i], ends[i] = bedStart(l), bedEnd(l)
	}
	rng := rand.New(rand.NewSource(11))
	for q := 0; q < 300; q++ {
		start := rng.Intn(200_000_000)
		width := 1 + rng.Intn(16000) // within one 16 kb window: the cached path
		if q%3 == 0 {
			width = 1 + rng.Intn(3_000_000) // wider: the direct chunk scan
		}
		end := start + width
		var want []string
		for i, l := range lines {
			if starts[i] < end && ends[i] > start {
				want = append(want, l)
			}
		}
		got := queryLines(t, r, "chr1", start, end)
		sort.Strings(want)
		sort.Strings(got)
		if strings.Join(got, "\n") != strings.Join(want, "\n") {
			t.Fatalf("chr1:%d-%d: got %d records, want %d", start, end, len(got), len(want))
		}
	}
}

// A small window's chunks span a small part of the file. The writer used to
// extend each bin's one chunk over every record of that bin, so a top-level
// bin's chunk ran from its first long feature to its last -- most of the file
// -- and every query, which always reads the top-level bins, read it.
func TestIndexChunksAreLocal(t *testing.T) {
	path, _ := longFeatureBED(t)
	st, err := os.Stat(path)
	if err != nil {
		t.Fatal(err)
	}
	idx, err := LoadTBI(path + ".tbi")
	if err != nil {
		t.Fatal(err)
	}
	ref := idx.RefID("chr1")
	worst := int64(0)
	for start := 0; start < 200_000_000; start += 1_000_000 {
		var span int64
		for _, c := range idx.Query(ref, start, start+16384) {
			span += c.End.BlockOffset() - c.Begin.BlockOffset()
		}
		worst = max(worst, span)
	}
	if limit := st.Size() / 10; worst > limit {
		t.Errorf("a 16 kb window's chunks span %d compressed bytes of a %d-byte file; want under %d", worst, st.Size(), limit)
	}
}
