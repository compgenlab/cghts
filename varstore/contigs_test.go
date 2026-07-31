package varstore

import (
	"path/filepath"
	"strings"
	"testing"
)

// TestContigsRoundTrip pins that the source's ##contig lines survive conversion.
//
// They are kept because a store is expected to be exported back to VCF, and those
// lines are how a VCF says which reference it was called against. Derived from the
// calls alone the best available is a bare ID for whichever contigs a query
// happened to name -- no lengths, and nothing at all for a query that named none.
func TestContigsRoundTrip(t *testing.T) {
	want := []string{
		"##contig=<ID=chr1,length=248956422>",
		"##contig=<ID=chr2,length=242193529>",
	}
	base := filepath.Join(t.TempDir(), "s")
	w, err := NewWriter(base, WriterOpts{
		Samples: []string{"S1"},
		MinDP:   10,
		Contigs: want,
	})
	if err != nil {
		t.Fatal(err)
	}
	if err := w.WriteCall(Call{
		SampleID: "S1", Chrom: "chr1", Pos: 100, Ref: "A", Alt: "G",
		GT: "0/1", DP: 30, ADRef: Missing, ADAlt: Missing, GQ: 99,
	}); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	s, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	got := s.Contigs()
	if strings.Join(got, "|") != strings.Join(want, "|") {
		t.Errorf("Contigs() = %v\nwant %v", got, want)
	}
	// Lengths are the part a query cannot reconstruct, so check one explicitly.
	if len(got) > 0 && !strings.Contains(got[0], "length=248956422") {
		t.Errorf("contig length did not survive: %q", got[0])
	}
}

// TestContigsAbsentIsNotAnError pins that a store written without contigs -- every
// store made before this was recorded, and any source that declares none -- reads
// back as nothing rather than failing or yielding an empty string.
func TestContigsAbsentIsNotAnError(t *testing.T) {
	base := filepath.Join(t.TempDir(), "s")
	w, err := NewWriter(base, WriterOpts{Samples: []string{"S1"}, MinDP: 10})
	if err != nil {
		t.Fatal(err)
	}
	if err := w.WriteCall(Call{
		SampleID: "S1", Chrom: "chr1", Pos: 100, Ref: "A", Alt: "G",
		GT: "0/1", DP: 30, ADRef: Missing, ADAlt: Missing, GQ: 99,
	}); err != nil {
		t.Fatal(err)
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	s, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	if got := s.Contigs(); len(got) != 0 {
		t.Errorf("Contigs() = %v, want nothing -- an unset key must not become [\"\"]", got)
	}
}
