package gtstore

import (
	"strings"
	"testing"
)

// IgnoreRefBlocks answers as though no blocks had been attached, so one store
// can be compared with and without its gVCF blocks.
func TestIgnoreRefBlocksOnAStore(t *testing.T) {
	path := writeAnnotVCF(t)
	blocks := []RefBlock{
		{SampleID: "S1", Chrom: "chr1", Start: 150, End: 250, MinDP: 20},
		{SampleID: "S2", Chrom: "chr1", Start: 150, End: 250, MinDP: 20},
	}
	ps, err := OpenParquet(convertForTest(t, path, blocks))
	if err != nil {
		t.Fatal(err)
	}
	defer ps.Close()

	off := Locus{Chrom: "chr1", Pos: 180, Ref: "A", Alt: "T"}      // between catalog sites
	excluded := Locus{Chrom: "chr1", Pos: 200, Ref: "C", Alt: "T"} // imputed, filtered out
	catalog := Locus{Chrom: "chr1", Pos: 100, Ref: "A", Alt: "G"}  // answered by runs
	base := Query{Loci: []Locus{off, excluded, catalog}, IncludeRef: true, Gate: annotGate,
		ExcludeInfo: []string{"IMPUTED"}}

	with := queryRows(t, ps, base)
	ignoring := base
	ignoring.IgnoreRefBlocks = true
	without := queryRows(t, ps, ignoring)

	count := func(rows []string, src string) int {
		n := 0
		for _, r := range rows {
			if strings.Contains(r, src) {
				n++
			}
		}
		return n
	}
	if count(with, "[ref_block]") == 0 {
		t.Fatalf("fixture answers nothing from blocks:\n  %s", strings.Join(with, "\n  "))
	}
	if n := count(without, "[ref_block]"); n != 0 {
		t.Errorf("IgnoreRefBlocks still answered %d row(s) from blocks:\n  %s", n, strings.Join(without, "\n  "))
	}
	// Everything not from a block is unchanged.
	var keep []string
	for _, r := range with {
		if !strings.Contains(r, "[ref_block]") {
			keep = append(keep, r)
		}
	}
	if strings.Join(keep, "\n") != strings.Join(without, "\n") {
		t.Errorf("catalog answers changed:\nwith blocks (minus block rows):\n  %s\nignoring:\n  %s",
			strings.Join(keep, "\n  "), strings.Join(without, "\n  "))
	}
}

// The same on a gVCF: its block records are skipped, its variant records kept.
func TestIgnoreRefBlocksOnAGvcf(t *testing.T) {
	_, indexed := writeGvcf(t)
	q := Query{Spans: []Span{{Chrom: "chr1", Start: 0, End: 30000}}, IncludeRef: true}
	with := rowsOf(t, indexed, q)
	q.IgnoreRefBlocks = true
	without := rowsOf(t, indexed, q)
	blocks := 0
	for _, c := range with {
		if c.Alt == "." {
			blocks++
		}
	}
	if blocks == 0 {
		t.Fatal("fixture has no block rows")
	}
	for _, c := range without {
		if c.Alt == "." || c.RefSource == RefFromBlock {
			t.Errorf("IgnoreRefBlocks still reported a block row: %+v", c)
		}
	}
	if len(without) != len(with)-blocks {
		t.Errorf("%d rows ignoring blocks, want %d (all but the %d block rows)", len(without), len(with)-blocks, blocks)
	}
}
