package tabix

import (
	"bytes"
	"fmt"
	"math/rand"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/htsio/bgzf"
)

// The tabix reader and index writer are checked here against htslib's tabix,
// the reference implementation. For every fixture -- BED, GTF and VCF, laid out
// plainly, as concatenated BGZF streams and as one block per line -- each
// region is answered four ways and three must equal the first:
//
//	A  htslib tabix on an htslib-written index   (the reference)
//	B  htslib tabix on a cghts-written index     (the writer, checked by htslib)
//	C  cghts reader on an htslib-written index   (the reader, on chunks it did not write)
//	D  cghts reader on a cghts-written index     (the two together)
//
// Regions are random -- within one 16 kb window (the reader's cached path) and
// wider (its direct scan) -- plus the edges: whole chromosomes, position 0,
// past a chromosome's end, and each side of record boundaries. And a cghts
// index must be about as local as htslib's: a 16 kb window's chunks spanning
// most of the file is what made every query slow (compgenlab/cghts#62).
//
// The tests skip without tabix on PATH; CI installs it. Slow by design: this
// code changes rarely, and when it does it should be checked hard.

func htslibTabix(t *testing.T) string {
	t.Helper()
	p, err := exec.LookPath("tabix")
	if err != nil {
		t.Skip("htslib tabix not on PATH")
	}
	return p
}

type contig struct {
	name   string
	length int
}

type tabixFixture struct {
	name    string
	preset  string      // tabix -p
	opts    *WriterOpts // the same configuration for cghts's writer
	header  []string
	contigs []contig
	lines   []string // records, sorted by contig then start
	layout  string   // "plain", "concat" (header and body as separate BGZF streams) or "blocks" (one block per line)
}

// Long features land in the top-level bins, interleaved with everything else:
// the layout that matters for chunk handling.
func bedRecords(rng *rand.Rand, contigs []contig, per int) []string {
	var lines []string
	for _, c := range contigs {
		step := c.length / per
		var cl []string
		for i := 0; i < per; i++ {
			start := i * step
			cl = append(cl, fmt.Sprintf("%s\t%d\t%d\tshort_%d_%x", c.name, start, start+50+rng.Intn(500), i, rng.Uint32()))
			if i%97 == 0 {
				end := min(c.length, start+1_000_000+rng.Intn(c.length/2))
				cl = append(cl, fmt.Sprintf("%s\t%d\t%d\tlong_%d_%x", c.name, start, end, i, rng.Uint32()))
			}
		}
		sort.SliceStable(cl, func(i, j int) bool { return field(cl[i], 1) < field(cl[j], 1) })
		lines = append(lines, cl...)
	}
	return lines
}

func gtfRecords(rng *rand.Rand, contigs []contig, genes int) []string {
	var lines []string
	for _, c := range contigs {
		step := c.length / genes
		var cl []string
		for g := 0; g < genes; g++ {
			start := g*step + 1 // GTF is 1-based
			span := 2000 + rng.Intn(50_000)
			if g%40 == 0 {
				span = 1_000_000 + rng.Intn(c.length/3)
			}
			end := min(c.length, start+span)
			strand := "+-"[rng.Intn(2):][:1]
			attr := fmt.Sprintf(`gene_id "G%s_%d"; transcript_id "T%s_%d";`, c.name, g, c.name, g)
			cl = append(cl, fmt.Sprintf("%s\tt\tgene\t%d\t%d\t.\t%s\t.\t%s", c.name, start, end, strand, attr))
			cl = append(cl, fmt.Sprintf("%s\tt\ttranscript\t%d\t%d\t.\t%s\t.\t%s", c.name, start, end, strand, attr))
			for e := start; e < end; e += 1 + rng.Intn(max(1, (end-start)/3)) {
				cl = append(cl, fmt.Sprintf("%s\tt\texon\t%d\t%d\t.\t%s\t.\t%s", c.name, e, min(end, e+100+rng.Intn(300)), strand, attr))
			}
		}
		sort.SliceStable(cl, func(i, j int) bool { return field(cl[i], 3) < field(cl[j], 3) })
		lines = append(lines, cl...)
	}
	return lines
}

// VCF ends come from the record: len(REF), or INFO/END for a gVCF block. (Not
// SVLEN: cghts reads it for symbolic alleles and htslib's tabix does not, a
// deliberate difference this would report.)
func vcfRecords(rng *rand.Rand, contigs []contig, per int) []string {
	var lines []string
	for _, c := range contigs {
		step := c.length / per
		for i := 0; i < per; i++ {
			pos := i*step + 1
			switch {
			case i%50 == 0:
				end := min(c.length, pos+500_000+rng.Intn(c.length/3))
				lines = append(lines, fmt.Sprintf("%s\t%d\t.\tA\t<NON_REF>\t.\t.\tEND=%d\tGT\t0/0", c.name, pos, end))
			case i%7 == 0:
				lines = append(lines, fmt.Sprintf("%s\t%d\t.\tA%s\tA\t50\tPASS\t.\tGT\t0/1", c.name, pos, strings.Repeat("C", 1+rng.Intn(400))))
			default:
				lines = append(lines, fmt.Sprintf("%s\t%d\t.\tA\tG\t50\tPASS\t.\tGT\t0/1", c.name, pos))
			}
		}
	}
	return lines
}

// field returns tab-separated column i as an int.
func field(line string, i int) int {
	var n int
	fmt.Sscanf(strings.Split(line, "\t")[i], "%d", &n)
	return n
}

func tabixFixtures() []tabixFixture {
	rng := rand.New(rand.NewSource(29))
	three := []contig{{"chr1", 200_000_000}, {"chr2", 50_000_000}, {"chrX", 8_000_000}}
	// Sized so a chunk spanning a chromosome is megabytes, far past a block's
	// slack: a small fixture lets a chromosome-wide chunk pass for local.
	bed := bedRecords(rng, three, 40000)
	gtf := gtfRecords(rng, three, 3000)
	vcfHeader := []string{"##fileformat=VCFv4.2",
		"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">",
		"##contig=<ID=chr1,length=200000000>", "##contig=<ID=chr2,length=50000000>", "##contig=<ID=chrX,length=8000000>",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"}
	vcf := vcfRecords(rng, three, 30000)
	gtfHeader := []string{"#description: test", "#format: gtf"}
	return []tabixFixture{
		{"bed", "bed", NewWriterOpts().BED(), nil, three, bed, "plain"},
		{"bed-blocks", "bed", NewWriterOpts().BED(), nil, three[1:], bedRecords(rng, three[1:], 600), "blocks"},
		{"gtf", "gff", NewWriterOpts().GFF().Meta('#'), gtfHeader, three, gtf, "plain"},
		{"gtf-concat", "gff", NewWriterOpts().GFF().Meta('#'), gtfHeader, three, gtf, "concat"},
		{"vcf", "vcf", NewWriterOpts().VCF(), vcfHeader, three, vcf, "plain"},
		{"vcf-concat", "vcf", NewWriterOpts().VCF(), vcfHeader, three, vcf, "concat"},
	}
}

func writeFixture(t *testing.T, path string, fx tabixFixture) {
	t.Helper()
	var out bytes.Buffer
	stream := func(lines []string, perLine bool) {
		w := bgzf.NewWriter(&out)
		for _, l := range lines {
			if _, err := w.Write([]byte(l + "\n")); err != nil {
				t.Fatal(err)
			}
			if perLine {
				if err := w.Flush(); err != nil {
					t.Fatal(err)
				}
			}
		}
		if err := w.Close(); err != nil {
			t.Fatal(err)
		}
	}
	switch fx.layout {
	case "concat":
		stream(fx.header, false)
		stream(fx.lines, false)
	case "blocks":
		stream(append(append([]string{}, fx.header...), fx.lines...), true)
	default:
		stream(append(append([]string{}, fx.header...), fx.lines...), false)
	}
	if err := os.WriteFile(path, out.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}
}

type region struct {
	ref        string
	start, end int // 0-based half-open; end < 0 is the whole contig
}

func (r region) String() string {
	if r.end < 0 {
		return r.ref
	}
	return fmt.Sprintf("%s:%d-%d", r.ref, r.start+1, r.end)
}

func fixtureRegions(rng *rand.Rand, fx tabixFixture) []region {
	var out []region
	for _, c := range fx.contigs {
		out = append(out,
			region{c.name, 0, -1},
			region{c.name, 0, 1},
			region{c.name, c.length - 10, c.length + 1_000_000},
			region{c.name, c.length + 1000, c.length + 2000})
		for i := 0; i < 150; i++ {
			start := rng.Intn(c.length)
			width := 1 + rng.Intn(16000)
			if i%3 == 0 {
				width = 1 + rng.Intn(c.length/10)
			}
			out = append(out, region{c.name, start, start + width})
		}
	}
	// Each side of real record boundaries.
	for i := 0; i < 60; i++ {
		l := fx.lines[rng.Intn(len(fx.lines))]
		ref := strings.SplitN(l, "\t", 2)[0]
		var s int
		switch fx.preset {
		case "bed":
			s = field(l, 1)
		case "gff":
			s = field(l, 3) - 1
		default:
			s = field(l, 1) - 1
		}
		out = append(out, region{ref, s, s + 1}, region{ref, max(0, s-1), s}, region{ref, s + 1, s + 2})
	}
	return out
}

func tabixCLI(t *testing.T, bin, path string, r region) []string {
	t.Helper()
	out, err := exec.Command(bin, path, r.String()).Output()
	if err != nil {
		t.Fatalf("tabix %s %s: %v", filepath.Base(path), r, err)
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

func cghtsQuery(t *testing.T, rd *Reader, r region) []string {
	t.Helper()
	end := r.end
	if end < 0 {
		end = 1 << 29
	}
	got := queryLines(t, rd, r.ref, r.start, end)
	sort.Strings(got)
	return got
}

func TestTabixAgreesWithHtslib(t *testing.T) {
	bin := htslibTabix(t)
	for _, fx := range tabixFixtures() {
		t.Run(fx.name, func(t *testing.T) {
			dir := t.TempDir()
			hts := filepath.Join(dir, "hts.gz")
			cg := filepath.Join(dir, "cghts.gz")
			writeFixture(t, hts, fx)
			writeFixture(t, cg, fx)
			if out, err := exec.Command(bin, "-p", fx.preset, hts).CombinedOutput(); err != nil {
				t.Fatalf("tabix -p %s: %v\n%s", fx.preset, err, out)
			}
			if err := NewIndexWriter(fx.opts).WriteIndex(cg); err != nil {
				t.Fatalf("cghts index: %v", err)
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

			regions := fixtureRegions(rand.New(rand.NewSource(31)), fx)
			nonEmpty := 0
			for _, r := range regions {
				want := tabixCLI(t, bin, hts, r)
				if len(want) > 0 {
					nonEmpty++
				}
				check := func(what string, got []string) {
					if strings.Join(got, "\n") != strings.Join(want, "\n") {
						t.Fatalf("%s %s: %d records, htslib on its own index %d", what, r, len(got), len(want))
					}
				}
				check("htslib on the cghts index", tabixCLI(t, bin, cg, r))
				check("cghts on the htslib index", cghtsQuery(t, onHts, r))
				check("cghts on its own index", cghtsQuery(t, onCg, r))
			}
			if nonEmpty < len(regions)/4 {
				t.Fatalf("only %d of %d regions returned anything; the fixture is not exercising queries", nonEmpty, len(regions))
			}

			// Locality: the worst 16 kb window's chunks, cghts's index against htslib's.
			worst := func(path string) int64 {
				idx, err := LoadTBI(path + ".tbi")
				if err != nil {
					t.Fatal(err)
				}
				var w int64
				for _, c := range fx.contigs {
					ref := idx.RefID(c.name)
					for start := 0; start < c.length; start += c.length / 200 {
						var span int64
						for _, ch := range idx.Query(ref, start, start+16384) {
							span += ch.End.BlockOffset() - ch.Begin.BlockOffset()
						}
						w = max(w, span)
					}
				}
				return w
			}
			cw, hw := worst(cg), worst(hts)
			st, _ := os.Stat(cg)
			t.Logf("worst 16 kb window: cghts index %d, htslib index %d compressed bytes; file %d", cw, hw, st.Size())
			// A block's slack either way: the writers may cut a chunk on either
			// side of a block boundary.
			if limit := 2*hw + 65536; cw > limit {
				t.Errorf("worst 16 kb window spans %d compressed bytes in cghts's index, %d in htslib's; want at most %d", cw, hw, limit)
			}
		})
	}
}

// A line ending exactly at a full block's end: the next line starts the next
// block, and the index must say so. It used to record that line's start as
// byte 0 of the block before -- the middle of a line -- because the offset
// 65536 wrapped; the record could not be found, and htslib, reading the cghts
// index, failed outright.
func TestALineEndingAtAFullBlockIsFound(t *testing.T) {
	// Lines of exactly 64 bytes: 1024 of them fill the first block exactly.
	// The rest are on chr2, so the index must start a chunk -- not extend one
	// -- at the boundary, which is where a wrapped offset does its damage.
	var lines []string
	for i := 0; i < 1100; i++ {
		chrom := "chr1"
		if i >= 1024 {
			chrom = "chr2"
		}
		l := fmt.Sprintf("%s\t%d\t%d\t", chrom, i*100, i*100+50)
		l += strings.Repeat("n", 63-len(l))
		lines = append(lines, l)
	}
	path := filepath.Join(t.TempDir(), "edge.bed.gz")
	writeBGZF(t, path, lines...)
	if err := NewIndexWriter(NewWriterOpts().BED()).WriteIndex(path); err != nil {
		t.Fatal(err)
	}
	r, err := NewReader(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	for _, i := range []int{1023, 1024, 1025} { // last of block 1, first and second of block 2
		got := queryLines(t, r, strings.Fields(lines[i])[0], i*100, i*100+1)
		if len(got) != 1 || got[0] != lines[i] {
			t.Errorf("record %d (around the block boundary): got %q", i, got)
		}
	}
	if bin, err := exec.LookPath("tabix"); err == nil {
		for _, i := range []int{1023, 1024, 1025} {
			got := tabixCLI(t, bin, path, region{strings.Fields(lines[i])[0], i * 100, i*100 + 1})
			if len(got) != 1 || got[0] != lines[i] {
				t.Errorf("htslib on the cghts index, record %d: got %q", i, got)
			}
		}
	}
}
