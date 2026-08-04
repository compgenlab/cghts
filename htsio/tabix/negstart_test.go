package tabix

import "testing"

// A negative start is a legitimate input: an annotator extending a window past
// the start of a contig produces one, as does ParseRegion on "chr1:0-1000".
// It used to index the linear index at -1 and panic -- and only for TBI, since
// the CSI path converts to uint32 and merely misses, so the same query crashed
// with one index kind and succeeded with the other.
func TestQueryClampsNegativeStart(t *testing.T) {
	tr, err := NewReader("../../vcf/testdata/sample.vcf.gz")
	if err != nil {
		t.Fatal(err)
	}
	defer tr.Close()
	names := tr.RefNames()
	if len(names) == 0 {
		t.Fatal("fixture has no references")
	}

	count := func(start, end int) int {
		seq, err := tr.Query(names[0], start, end)
		if err != nil {
			t.Fatalf("Query(%d, %d): %v", start, end, err)
		}
		n := 0
		for _, err := range seq {
			if err != nil {
				t.Fatal(err)
			}
			n++
		}
		return n
	}

	// Clamping to 0 must not change the answer: [-100, end) and [0, end) select
	// the same records, since no record can start before zero.
	from0 := count(0, 1_000_000)
	if from0 == 0 {
		t.Fatal("fixture yielded no records, so this proves nothing")
	}
	if got := count(-100, 1_000_000); got != from0 {
		t.Errorf("negative start gave %d records, [0,end) gave %d", got, from0)
	}
}
