package tabix

import (
	"fmt"
	"math/rand"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

// htslibTabix returns the path of htslib's tabix, skipping the test when it is
// not installed. CI installs it; locally these tests run when it is present.
func htslibTabix(t *testing.T) string {
	t.Helper()
	p, err := exec.LookPath("tabix")
	if err != nil {
		t.Skip("htslib tabix not on PATH")
	}
	return p
}

// tabixLines runs `tabix FILE chr1:START-END` (1-based inclusive) for the
// 0-based half-open [start, end) and returns its lines, sorted.
func tabixLines(t *testing.T, bin, path string, start, end int) []string {
	t.Helper()
	out, err := exec.Command(bin, path, fmt.Sprintf("chr1:%d-%d", start+1, end)).Output()
	if err != nil {
		t.Fatalf("tabix %s chr1:%d-%d: %v", path, start+1, end, err)
	}
	var lines []string
	for _, l := range strings.Split(strings.TrimRight(string(out), "\n"), "\n") {
		if l != "" {
			lines = append(lines, l)
		}
	}
	sort.Strings(lines)
	return lines
}

// randomRegions yields n regions over the long-feature layout: most within one
// 16 kb window (the reader's cached path), a third wider (its direct scan).
func randomRegions(n int) [][2]int {
	rng := rand.New(rand.NewSource(23))
	var out [][2]int
	for i := 0; i < n; i++ {
		start := rng.Intn(200_000_000)
		width := 1 + rng.Intn(16000)
		if i%3 == 0 {
			width = 1 + rng.Intn(3_000_000)
		}
		out = append(out, [2]int{start, start + width})
	}
	return out
}

// cghts reading an index htslib wrote returns what htslib returns. The reader
// now follows chunk ends; htslib's chunks are its own, so this checks the
// reader against chunk lists cghts did not produce.
func TestReaderMatchesHtslibOnAnHtslibIndex(t *testing.T) {
	bin := htslibTabix(t)
	src, _ := longFeatureBED(t)
	path := filepath.Join(t.TempDir(), "htslib.bed.gz")
	b, err := os.ReadFile(src)
	if err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(path, b, 0o644); err != nil {
		t.Fatal(err)
	}
	if out, err := exec.Command(bin, "-p", "bed", path).CombinedOutput(); err != nil {
		t.Fatalf("tabix -p bed: %v\n%s", err, out)
	}

	r, err := NewReader(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	for _, q := range randomRegions(100) {
		want := tabixLines(t, bin, path, q[0], q[1])
		got := queryLines(t, r, "chr1", q[0], q[1])
		sort.Strings(got)
		if strings.Join(got, "\n") != strings.Join(want, "\n") {
			t.Fatalf("chr1:%d-%d on an htslib index: cghts %d records, htslib %d", q[0], q[1], len(got), len(want))
		}
	}
}

// htslib reading an index cghts wrote returns exactly the overlapping records.
// This is the writer checked by htslib's own reader: a chunk list that missed
// a record, or pointed somewhere wrong, would show here.
func TestHtslibReadsACghtsIndex(t *testing.T) {
	bin := htslibTabix(t)
	path, lines := longFeatureBED(t) // indexed by cghts
	starts := make([]int, len(lines))
	ends := make([]int, len(lines))
	for i, l := range lines {
		starts[i], ends[i] = bedStart(l), bedEnd(l)
	}
	for _, q := range randomRegions(100) {
		var want []string
		for i, l := range lines {
			if starts[i] < q[1] && ends[i] > q[0] {
				want = append(want, l)
			}
		}
		sort.Strings(want)
		got := tabixLines(t, bin, path, q[0], q[1])
		if strings.Join(got, "\n") != strings.Join(want, "\n") {
			t.Fatalf("chr1:%d-%d on a cghts index: htslib %d records, want %d", q[0], q[1], len(got), len(want))
		}
	}
}

// cghts's index reads about as little as htslib's. For the same file, the
// compressed span a 16 kb window's chunks cover in each index is compared;
// the writer used to give each bin one chunk from its first record to its
// last, so a window covered most of the file where htslib's covered a block.
func TestIndexChunksAreAsLocalAsHtslib(t *testing.T) {
	bin := htslibTabix(t)
	src, _ := longFeatureBED(t) // cghts index beside it
	htsPath := filepath.Join(t.TempDir(), "htslib.bed.gz")
	b, err := os.ReadFile(src)
	if err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(htsPath, b, 0o644); err != nil {
		t.Fatal(err)
	}
	if out, err := exec.Command(bin, "-p", "bed", htsPath).CombinedOutput(); err != nil {
		t.Fatalf("tabix -p bed: %v\n%s", err, out)
	}
	worst := func(path string) int64 {
		idx, err := LoadTBI(path + ".tbi")
		if err != nil {
			t.Fatal(err)
		}
		ref := idx.RefID("chr1")
		var w int64
		for start := 0; start < 200_000_000; start += 1_000_000 {
			var span int64
			for _, c := range idx.Query(ref, start, start+16384) {
				span += c.End.BlockOffset() - c.Begin.BlockOffset()
			}
			w = max(w, span)
		}
		return w
	}
	cg, hts := worst(src), worst(htsPath)
	// A block's slack either way: the two writers may cut a chunk on either
	// side of a block boundary.
	if limit := 2*hts + 65536; cg > limit {
		t.Errorf("worst 16 kb window spans %d compressed bytes in cghts's index, %d in htslib's; want at most %d", cg, hts, limit)
	}
}
