package tabix

import (
	"path/filepath"
	"testing"
)

// htslib records each reference's record count in a pseudo-bin; these fixtures
// were indexed by htslib (4 records on chr1, 1 on chr2).
func TestRecordCountFromHtslibIndexes(t *testing.T) {
	for _, path := range []string{"../testdata/test.vcf.gz", "../testdata/test_csi.vcf.gz"} {
		r, err := NewReader(path)
		if err != nil {
			t.Fatal(err)
		}
		for ref, want := range map[string]uint64{"chr1": 4, "chr2": 1} {
			if n, ok := r.RecordCount(ref); !ok || n != want {
				t.Errorf("%s: RecordCount(%s) = %d, %v; want %d", path, ref, n, ok, want)
			}
		}
		if _, ok := r.RecordCount("chrNope"); ok {
			t.Errorf("%s: an unknown reference reported a count", path)
		}
		r.Close()
	}
}

// An index without the pseudo-bin -- one from cghts before it wrote one --
// reports no count rather than zero, so a caller can fall back.
func TestRecordCountAbsentIsNotZero(t *testing.T) {
	idx := &BinIndex{refs: []binIndex{{bins: map[uint32][]Chunk{4681: {{Begin: 1, End: 2}}}}}}
	if n, ok := idx.RecordCount(0); ok {
		t.Errorf("RecordCount = %d, true; want not ok", n)
	}
	if _, ok := idx.RecordCount(5); ok {
		t.Error("out-of-range reference reported a count")
	}
}

// cghts's own writer records counts too, as htslib does.
func TestIndexWriterRecordsCounts(t *testing.T) {
	path := filepath.Join(t.TempDir(), "v.vcf.gz")
	writeBGZF(t, path,
		"##fileformat=VCFv4.2",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
		"chr1\t100\t.\tA\tG\t.\tPASS\t.",
		"chr1\t200\t.\tA\tG\t.\tPASS\t.",
		"chr1\t70000\t.\tA\tG\t.\tPASS\t.",
		"chr2\t500\t.\tC\tT\t.\tPASS\t.")
	if err := NewIndexWriter(NewWriterOpts().VCF()).WriteIndex(path); err != nil {
		t.Fatal(err)
	}
	r, err := NewReader(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	for ref, want := range map[string]uint64{"chr1": 3, "chr2": 1} {
		if n, ok := r.RecordCount(ref); !ok || n != want {
			t.Errorf("RecordCount(%s) = %d, %v; want %d", ref, n, ok, want)
		}
	}
	// The pseudo-bin must not disturb queries.
	if got := queryLines(t, r, "chr1", 0, 100000); len(got) != 3 {
		t.Errorf("chr1 query = %d rows, want 3", len(got))
	}
}
