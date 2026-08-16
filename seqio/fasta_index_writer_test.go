package seqio

import (
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"math/rand"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/htsio/bgzf"
)

const idxFasta = ">chr1 with a description\n" +
	"ACGTACGTACGTACGTACGT\nACGTACGTACGTACGTACGT\nACGTACGT\n" +
	">chr2\n" +
	"TTTTAAAACCCCGGGG\nTTTTAAAACCCCGGGG\nTTTT\n" +
	">chr3\n" +
	"NNNNNNNNNN\n"

// bigFasta spans several BGZF blocks, so a .gzi carries real entries. A small
// fixture only ever proves the count is zero.
func bigFasta() string {
	var b strings.Builder
	line := strings.Repeat("ACGT", 15) // 60 bases
	for c := 1; c <= 3; c++ {
		fmt.Fprintf(&b, ">chr%d some description\n", c)
		for i := 0; i < 1500; i++ {
			b.WriteString(line)
			b.WriteByte('\n')
		}
		b.WriteString("ACGTACGT\n") // short final line
	}
	return b.String()
}

// bgzipTo compresses body into path using this package's own BGZF writer.
func bgzipTo(t *testing.T, path, body string) {
	t.Helper()
	w, err := bgzf.NewBGZipFile(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := w.Write([]byte(body)); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
}

// samtools is the oracle. Our indexes are only correct if htslib's readers
// accept them, and the cheapest proof is producing byte-identical output to
// htslib's own writer.
func samtools(t *testing.T) string {
	t.Helper()
	p, err := exec.LookPath("samtools")
	if err != nil {
		t.Skip("samtools not installed; skipping the htslib byte comparison")
	}
	return p
}

// The index has to be usable, which means the reader in this package must read
// correct bases through it — including across a line boundary, where a wrong
// LineWidth shows up.
func TestBuildFastaIndexRoundTrips(t *testing.T) {
	dir := t.TempDir()
	p := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(p, []byte(idxFasta), 0o644); err != nil {
		t.Fatal(err)
	}
	entries, err := BuildFastaIndex(p)
	if err != nil {
		t.Fatal(err)
	}
	if len(entries) != 3 {
		t.Fatalf("indexed %d sequences, want 3", len(entries))
	}
	if entries[0].Name != "chr1" || entries[0].Length != 48 {
		t.Errorf("chr1 = %+v, want name chr1 length 48", entries[0])
	}

	r, err := NewIndexedFastaReader(p)
	if err != nil {
		t.Fatalf("our own reader rejected the index: %v", err)
	}
	defer r.Close()
	got, err := r.GetSequenceRange("chr2", 14, 18)
	if err != nil {
		t.Fatal(err)
	}
	if strings.ToUpper(string(got)) != "GGTT" {
		t.Errorf("chr2:14-18 = %q, want GGTT (spans a line boundary)", got)
	}
}

func TestFaiMatchesHtslib(t *testing.T) {
	sam := samtools(t)
	// Separate directories: both writers put the index beside the file.
	ours, theirs := t.TempDir(), t.TempDir()
	for _, d := range []string{ours, theirs} {
		if err := os.WriteFile(filepath.Join(d, "ref.fa"), []byte(idxFasta), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	if out, err := exec.Command(sam, "faidx", filepath.Join(theirs, "ref.fa")).CombinedOutput(); err != nil {
		t.Fatalf("samtools faidx: %v\n%s", err, out)
	}
	if _, err := BuildFastaIndex(filepath.Join(ours, "ref.fa")); err != nil {
		t.Fatal(err)
	}
	want, _ := os.ReadFile(filepath.Join(theirs, "ref.fa.fai"))
	got, _ := os.ReadFile(filepath.Join(ours, "ref.fa.fai"))
	if string(got) != string(want) {
		t.Errorf(".fai differs from htslib.\n--- ours ---\n%s--- samtools ---\n%s", got, want)
	}
}

// Both index the same compressed bytes, so the block boundaries are fixed and
// the .gzi must agree exactly — this catches an off-by-one in the block walk
// that a "does it load" check would not. The EOF marker is the specific trap:
// counting it yields one entry too many.
func TestFaiAndGziMatchHtslibForBGZF(t *testing.T) {
	sam := samtools(t)

	src := t.TempDir()
	gz := filepath.Join(src, "ref.fa.gz")
	bgzipTo(t, gz, bigFasta())
	body, err := os.ReadFile(gz)
	if err != nil {
		t.Fatal(err)
	}

	ours, theirs := t.TempDir(), t.TempDir()
	for _, d := range []string{ours, theirs} {
		if err := os.WriteFile(filepath.Join(d, "ref.fa.gz"), body, 0o644); err != nil {
			t.Fatal(err)
		}
	}
	if o, err := exec.Command(sam, "faidx", filepath.Join(theirs, "ref.fa.gz")).CombinedOutput(); err != nil {
		t.Fatalf("samtools faidx (bgzf): %v\n%s", err, o)
	}
	if _, err := BuildFastaIndex(filepath.Join(ours, "ref.fa.gz")); err != nil {
		t.Fatal(err)
	}
	for _, ext := range []string{".fai", ".gzi"} {
		want, err := os.ReadFile(filepath.Join(theirs, "ref.fa.gz"+ext))
		if err != nil {
			t.Fatalf("samtools wrote no %s: %v", ext, err)
		}
		got, err := os.ReadFile(filepath.Join(ours, "ref.fa.gz"+ext))
		if err != nil {
			t.Fatalf("we wrote no %s: %v", ext, err)
		}
		if !bytes.Equal(got, want) {
			t.Errorf("%s differs from htslib (ours %d bytes, samtools %d)", ext, len(got), len(want))
		}
	}
}

// A .fai addresses uncompressed offsets, so it is only usable with a container
// that can seek to one. Plain gzip cannot, and an unusable index is worse than
// a refusal.
func TestBuildFastaIndexRefusesPlainGzip(t *testing.T) {
	dir := t.TempDir()
	p := filepath.Join(dir, "ref.fa.gz")
	f, err := os.Create(p)
	if err != nil {
		t.Fatal(err)
	}
	zw := gzip.NewWriter(f) // no BGZF extra field
	if _, err := zw.Write([]byte(idxFasta)); err != nil {
		t.Fatal(err)
	}
	zw.Close()
	f.Close()

	if _, err := BuildFastaIndex(p); err == nil {
		t.Fatal("a plain-gzip FASTA was accepted; its index would be unusable")
	} else if !strings.Contains(err.Error(), "BGZF") {
		t.Errorf("error does not name the problem: %v", err)
	}
}

// A ragged record cannot be described by a .fai, and producing one anyway would
// give wrong coordinates rather than an error.
func TestBuildFastaIndexRejectsRaggedLines(t *testing.T) {
	dir := t.TempDir()
	p := filepath.Join(dir, "ragged.fa")
	if err := os.WriteFile(p, []byte(">chr1\nACGTACGT\nACG\nACGTACGT\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if _, err := BuildFastaIndex(p); err == nil {
		t.Error("uneven line lengths were accepted")
	}
}

// Indexing must not scale with the file.
//
// The first version of this decompressed the whole stream into a bytes.Buffer to
// find block boundaries — which passes every correctness test and then dies on
// the only input anyone cares about: a human genome is ~3 GB uncompressed, and
// the worker was OOM-killed indexing GRCh38.
//
// HOW THIS IS MEASURED IS THE WHOLE DIFFICULTY, and three earlier versions of it
// failed correct code. Reading the heap after the call measured how recently a
// GC had run. Sampling from a goroutine during the call measured how much
// garbage had piled up between cycles, because HeapAlloc counts what the
// collector has not reached yet. Collecting aggressively with SetGCPercent
// narrowed that and did not close it: an aggressive GC target still needs CPU to
// meet, so on a loaded machine the same correct code read 4 MB idle and 15 MB
// busy, and tripped a 15 MB threshold. The instrument was at fault every time.
//
// So the probe now runs INSIDE Read, synchronously, on the implementation's own
// goroutine. Nothing else is allocating at that instant, so a forced collection
// followed by a reading gives the true live set. A busy machine makes this
// slower and cannot make it wrong, which is what the previous versions lacked.
//
// Making the reader the instrument also means these test indexFasta and the bgzf
// block indexer directly, rather than through the path BuildFastaIndex opens for
// itself. Both take an io.Reader, and both are where a buffering regression
// would actually land.

// probeReader samples the live heap on every Read of the stream beneath it.
type probeReader struct {
	r     io.Reader
	base  uint64
	peak  int64
	reads int
}

func (p *probeReader) Read(b []byte) (int, error) {
	// Synchronous, on the caller's goroutine: the implementation is between
	// allocations here, so a forced collection settles the heap to exactly what
	// is live. This is the difference between measuring the program and
	// measuring the scheduler.
	runtime.GC()
	var m runtime.MemStats
	runtime.ReadMemStats(&m)
	if grew := int64(m.HeapAlloc) - int64(p.base); grew > p.peak {
		p.peak = grew
	}
	p.reads++
	return p.r.Read(b)
}

func newProbe(r io.Reader) *probeReader {
	runtime.GC()
	var m runtime.MemStats
	runtime.ReadMemStats(&m)
	return &probeReader{r: r, base: m.HeapAlloc}
}

// hugeFasta returns a body large enough that buffering it is unmistakable.
// Distinct from bigFasta, which is a small three-sequence fixture.
//
// THE BASES ARE PSEUDO-RANDOM, not a repeated motif, and that is not
// decoration. A repeated motif compresses 61 MB down to 0.18 MB, which the
// block indexer's 1 MB buffer then swallows in a single Read -- so the probe
// samples once, before anything has been decoded, and reports a heap of nearly
// nothing however the implementation behaves. Incompressible input keeps the
// compressed stream large enough to be read in many pieces, which is what makes
// the measurement mean something. Seeded, so the fixture is identical run to
// run.
func hugeFasta(lines int) string {
	rng := rand.New(rand.NewSource(1))
	const bases = "ACGT"
	var b strings.Builder
	b.Grow(lines*61 + 16)
	b.WriteString(">chr1 big\n")
	// Bulk random, then two bits per base: rng.Intn per base costs more than
	// everything else in the test put together.
	line := make([]byte, 60)
	for i := 0; i < lines; i++ {
		rng.Read(line)
		for j := range line {
			line[j] = bases[line[j]&3]
		}
		b.Write(line)
		b.WriteByte('\n')
	}
	return b.String()
}

// liveCeiling is what a streaming implementation may hold.
//
// indexFasta reads through a fixed 1 MB bufio and the block indexer through
// another, so the working set is a couple of megabytes plus one decompressed
// block, whatever the input size. What is retained beyond that is O(sequences)
// and O(blocks), never O(bytes). 8 MB leaves room for allocator slack and still
// sits an order of magnitude under the 61 MB a buffering implementation needs.
const liveCeiling = 8 << 20

func TestIndexFastaDoesNotBufferTheStream(t *testing.T) {
	if testing.Short() {
		t.Skip("allocates a multi-megabyte fixture")
	}
	body := hugeFasta(1_000_000)
	p := newProbe(strings.NewReader(body))

	entries, err := indexFasta(p)
	if err != nil {
		t.Fatal(err)
	}
	if len(entries) != 1 {
		t.Fatalf("indexed %d sequences, want 1", len(entries))
	}
	// A fixture read in one gulp would prove nothing about streaming.
	if p.reads < 8 {
		t.Fatalf("only %d reads, so the probe never sampled mid-stream", p.reads)
	}
	if p.peak > liveCeiling {
		t.Errorf("live heap reached %d bytes indexing a %d-byte stream -- this is buffering it, "+
			"which OOMs on a real genome", p.peak, len(body))
	}
	t.Logf("peak live %.1f MB over %d reads of a %.1f MB stream",
		float64(p.peak)/1e6, p.reads, float64(len(body))/1e6)
}

func TestBlockIndexingDoesNotBufferTheStream(t *testing.T) {
	if testing.Short() {
		t.Skip("allocates a multi-megabyte fixture")
	}
	dir := t.TempDir()
	gz := filepath.Join(dir, "big.fa.gz")
	body := hugeFasta(300_000)
	bgzipTo(t, gz, body)

	f, err := os.Open(gz)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()

	// The probe sits UNDER the block indexer, so it samples while blocks are
	// decoded and their boundaries accumulated -- which is where the original
	// bug lived.
	p := newProbe(f)
	ir := bgzf.NewIndexingReader(p)
	entries, err := indexFasta(ir)
	if err != nil {
		t.Fatal(err)
	}
	if len(entries) != 1 {
		t.Fatalf("indexed %d sequences, want 1", len(entries))
	}
	if p.reads < 4 {
		t.Fatalf("only %d reads, so the probe never sampled mid-stream", p.reads)
	}
	if p.peak > liveCeiling {
		t.Errorf("live heap reached %d bytes indexing a %d-byte stream -- this is buffering it, "+
			"which OOMs on a real genome", p.peak, len(body))
	}
	if n := len(ir.Blocks()); n == 0 {
		t.Fatal("no blocks recorded, so the fixture is not bgzf")
	}
	t.Logf("peak live %.1f MB over %d reads, %d blocks recorded",
		float64(p.peak)/1e6, p.reads, len(ir.Blocks()))
}

// And nothing is retained once BuildFastaIndex returns: an index holding onto
// the blocks would keep a genome resident for as long as the caller kept it.
//
// This reading IS safe from outside, because it is taken after the call with
// nothing else running -- two forced collections settle it exactly.
func TestBuildFastaIndexRetainsNothing(t *testing.T) {
	if testing.Short() {
		t.Skip("allocates a multi-megabyte fixture")
	}
	dir := t.TempDir()
	gz := filepath.Join(dir, "big.fa.gz")
	body := hugeFasta(300_000)
	bgzipTo(t, gz, body)

	runtime.GC()
	var before runtime.MemStats
	runtime.ReadMemStats(&before)

	entries, err := BuildFastaIndex(gz)
	if err != nil {
		t.Fatal(err)
	}
	if len(entries) != 1 {
		t.Fatalf("indexed %d sequences, want 1", len(entries))
	}

	runtime.GC()
	runtime.GC()
	var after runtime.MemStats
	runtime.ReadMemStats(&after)
	uncompressed := int64(len(body))
	if held := int64(after.HeapAlloc) - int64(before.HeapAlloc); held > uncompressed/100 {
		t.Errorf("the index still holds %d bytes of a %d-byte file after it was built",
			held, uncompressed)
	}
	runtime.KeepAlive(entries)
}
