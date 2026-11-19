package gtstore

import (
	"os"
	"path/filepath"
	"sync/atomic"
	"testing"

	"github.com/parquet-go/parquet-go"
)

type countingReaderAt struct {
	f     *os.File
	reads atomic.Int64
}

func (c *countingReaderAt) ReadAt(p []byte, off int64) (int, error) {
	c.reads.Add(1)
	return c.f.ReadAt(p, off)
}

// Opening a store must not read every row group's bloom filter. parquet-go's
// default open does -- one random read per row group per filtered column --
// which on cluster storage left gtstore create waiting for minutes on volumes
// it only needed the footers of.
func TestOpeningATableReadsOnlyItsFooter(t *testing.T) {
	base := filepath.Join(t.TempDir(), "store")
	w, err := NewWriter(base, WriterOpts{Samples: []string{"S1", "S2"}, MinDP: 10, RowGroupSize: 10})
	if err != nil {
		t.Fatal(err)
	}
	for i := int32(0); i < 2000; i++ {
		for _, s := range []string{"S1", "S2"} {
			if err := w.WriteCall(Call{SampleID: s, Chrom: "chr1", Pos: 100 + i, Ref: "A", Alt: "G", GT: "0/1",
				DP: 30, ADRef: Missing, ADAlt: Missing, GQ: Missing}); err != nil {
				t.Fatal(err)
			}
		}
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
	f, err := os.Open(CallsPath(base))
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	st, _ := f.Stat()

	eager := &countingReaderAt{f: f}
	pf, err := parquet.OpenFile(eager, st.Size())
	if err != nil {
		t.Fatal(err)
	}
	groups := len(pf.RowGroups())
	if groups < 100 {
		t.Fatalf("fixture has %d row groups; the test needs many", groups)
	}
	footer := &countingReaderAt{f: f}
	if _, err := openFooter(footer, st.Size()); err != nil {
		t.Fatal(err)
	}
	scan := &countingReaderAt{f: f}
	if _, err := openForScan(scan, st.Size()); err != nil {
		t.Fatal(err)
	}
	if n := footer.reads.Load(); n > 4 {
		t.Errorf("openFooter made %d reads over %d row groups; want a handful (the default made %d)",
			n, groups, eager.reads.Load())
	}
	if n := scan.reads.Load(); n > 6 {
		t.Errorf("openForScan made %d reads over %d row groups; bloom filters should load only on use (default %d)",
			n, groups, eager.reads.Load())
	}
	// And the store still opens and answers.
	ps, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer ps.Close()
	calls, err := CollectCalls(ps, Query{Loci: []Locus{{Chrom: "chr1", Pos: 1500, Ref: "A", Alt: "G"}}})
	if err != nil || len(calls) != 2 {
		t.Errorf("locus query after a lazy open: %d calls, %v", len(calls), err)
	}
}
