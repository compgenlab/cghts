package gtf

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/htsio/tabix"
)

// indexGTF writes lines as a bgzipped, tabix-indexed GTF.
func indexGTF(t *testing.T, lines ...string) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), "g.gtf.gz")
	w := tabix.NewWriter(path, tabix.NewWriterOpts().Columns(1, 4, 5).Meta('#').AutoIndex())
	for _, l := range lines {
		if err := w.Write(l); err != nil {
			t.Fatalf("write %q: %v", l, err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	return path
}

func gene(chrom, id, name, extra string) string {
	attrs := `gene_id "` + id + `"; gene_name "` + name + `";`
	if extra != "" {
		attrs += " " + extra
	}
	return chrom + "\tsrc\tgene\t100\t200\t.\t+\t.\t" + attrs
}

// An indexed GTF answers Provides by sampling the file.
//
// The whole-file model knows because it has read everything; this one has read
// nothing until queried, which is the point of it. Before this it could not
// answer at all and the annotator assumed yes — declaring a biotype field for
// every GTF, including the ones that have none.
func TestAnIndexedSourceReportsWhatItCarries(t *testing.T) {
	withBio := indexGTF(t,
		gene("chr1", "g1", "A", `gene_biotype "protein_coding";`),
		gene("chr1", "g2", "B", `gene_biotype "lncRNA";`))
	s, err := NewIndexedAnnotationSource(withBio, nil)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	if !s.Provides("biotype") {
		t.Error("a GTF carrying biotypes reported none")
	}
	if s.Provides("status") {
		t.Error("a GTF carrying no gene_status reported one")
	}

	// And a GTF with neither says so, rather than assuming.
	without := indexGTF(t, gene("chr1", "g1", "A", ""), gene("chr1", "g2", "B", ""))
	s2, err := NewIndexedAnnotationSource(without, nil)
	if err != nil {
		t.Fatal(err)
	}
	defer s2.Close()
	if s2.Provides("biotype") {
		t.Error("a GTF with no biotypes reported some; the header would declare " +
			"a field no record ever carries")
	}
}

// The answer matches the whole-file model's, which is the contract.
//
// Two implementations of "what does this GTF carry" that disagree would put a
// different header on the same file depending on how it happened to be opened.
func TestBothSourcesAgreeOnWhatIsProvided(t *testing.T) {
	lines := []string{
		gene("chr1", "g1", "A", `gene_biotype "protein_coding"; gene_status "KNOWN";`),
		gene("chr2", "g2", "B", `gene_biotype "lncRNA";`),
	}
	path := indexGTF(t, lines...)

	plain := filepath.Join(t.TempDir(), "g.gtf")
	if err := os.WriteFile(plain, []byte(strings.Join(lines, "\n")+"\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	mem, err := NewAnnotationSource(plain, nil)
	if err != nil {
		t.Fatal(err)
	}
	idx, err := NewIndexedAnnotationSource(path, nil)
	if err != nil {
		t.Fatal(err)
	}
	defer idx.Close()

	for _, key := range []string{"biotype", "status", "gene_id", "gene_name", "strand"} {
		if got, want := idx.Provides(key), mem.Provides(key); got != want {
			t.Errorf("Provides(%q): indexed %v, in-memory %v — the two models "+
				"disagree about the same file", key, got, want)
		}
	}
	// Not vacuous: this fixture carries both optional fields.
	if !mem.Provides("biotype") || !mem.Provides("status") {
		t.Fatal("the fixture carries neither optional field, so the comparison " +
			"proves nothing")
	}
}

// A GTF whose first contig has no usable records is not evidence of absence.
//
// The sample says "cannot tell", and unknown resolves to declaring the field: a
// declared field never written is a header line nobody reads, while a written
// one never declared is a record a strict parser may reject.
func TestAnEmptySampleReportsUnknownRatherThanAbsent(t *testing.T) {
	// Only a malformed row, so nothing parses.
	path := indexGTF(t, "chr1\tsrc\tgene\t100\t200\t.\t+\t.\t")
	s, err := NewIndexedAnnotationSource(path, nil)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	if !s.Provides("biotype") {
		t.Error("an unreadable sample reported the field absent; unknown has to " +
			"resolve to declaring it")
	}
}

// The sample is taken once, not per call.
func TestTheSampleIsTakenOnce(t *testing.T) {
	s, err := NewIndexedAnnotationSource(
		indexGTF(t, gene("chr1", "g1", "A", `gene_biotype "protein_coding";`)), nil)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	first := s.Provides("biotype")
	if !s.probed {
		t.Fatal("the sample was not recorded as taken")
	}
	if second := s.Provides("biotype"); second != first {
		t.Errorf("two calls disagreed: %v then %v", first, second)
	}
}

// The sample stops after probeLines lines, and what lies past them is not seen.
//
// This is the cost of sampling rather than scanning, stated as a test rather
// than left to be discovered: a GTF whose first thousand lines carry no biotype
// but whose later rows do will be reported as carrying none, and its BIOTYPE
// field will not be declared.
//
// The trade is deliberate. The alternative is reading the whole file to answer
// one boolean, which costs exactly what the index exists to avoid — and a GTF
// mixing schemas partway through is not a thing any pipeline emits.
func TestTheSampleStopsAtTheLineBudget(t *testing.T) {
	// probeLines rows with no biotype, then one with. A probe that read past
	// its budget would find the last one.
	lines := make([]string, 0, probeLines+1)
	for i := 0; i < probeLines; i++ {
		lines = append(lines, geneAt("chr1", 100+i*10, "g", "A", ""))
	}
	lines = append(lines, geneAt("chr1", 100+probeLines*10, "gz", "Z", `gene_biotype "protein_coding";`))

	s, err := NewIndexedAnnotationSource(indexGTF(t, lines...), nil)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	if s.Provides("biotype") {
		t.Error("the probe read past its line budget; the bound is what keeps " +
			"this from scanning the file it was indexed to avoid scanning")
	}

	// And the same file, read wholly, does find it — so the fixture really does
	// carry a biotype and this is a difference in how far each looked, not a
	// difference in the data.
	plain := filepath.Join(t.TempDir(), "g.gtf")
	if err := os.WriteFile(plain, []byte(strings.Join(lines, "\n")+"\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	mem, err := NewAnnotationSource(plain, nil)
	if err != nil {
		t.Fatal(err)
	}
	if !mem.Provides("biotype") {
		t.Fatal("the whole-file model found no biotype either; the fixture is wrong")
	}
}

// geneAt is gene() at a chosen coordinate, so a test can write many rows.
func geneAt(chrom string, pos int, id, name, extra string) string {
	attrs := `gene_id "` + id + `"; gene_name "` + name + `";`
	if extra != "" {
		attrs += " " + extra
	}
	return chrom + "\tsrc\tgene\t" + itoa(pos) + "\t" + itoa(pos+50) + "\t.\t+\t.\t" + attrs
}

func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	var b []byte
	for n > 0 {
		b = append([]byte{byte('0' + n%10)}, b...)
		n /= 10
	}
	return string(b)
}
