package vcf

import (
	"io"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

// writeChunk writes a VCF chunk (header + the given records) to filename.
func writeChunk(t *testing.T, filename string, h *VcfHeader, recs []*VcfRecord) {
	t.Helper()
	w, err := OpenVcfWriter(filename)
	if err != nil {
		t.Fatalf("OpenVcfWriter(%s): %v", filename, err)
	}
	if err := w.WriteHeader(h); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for _, rec := range recs {
		if err := w.WriteRecord(rec); err != nil {
			t.Fatalf("WriteRecord: %v", err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
}

func TestChunkedVcfReader(t *testing.T) {
	r := openTestFile(t)
	h, err := r.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	all := readAll(t, r)
	r.Close()

	dir := t.TempDir()
	// Split the 5 sample records across three numbered chunks (2 + 2 + 1).
	writeChunk(t, filepath.Join(dir, "split.1.vcf.gz"), h, all[0:2])
	writeChunk(t, filepath.Join(dir, "split.2.vcf.gz"), h, all[2:4])
	writeChunk(t, filepath.Join(dir, "split.3.vcf.gz"), h, all[4:5])

	c, err := NewChunkedVcfReader(filepath.Join(dir, "split.1.vcf.gz"))
	if err != nil {
		t.Fatalf("NewChunkedVcfReader: %v", err)
	}
	defer c.Close()

	var got []string
	for {
		rec, err := c.NextRecord()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("NextRecord: %v", err)
		}
		got = append(got, rec.Chrom+":"+strconv.Itoa(rec.Pos))
	}
	want := []string{"chr1:100", "chr1:200", "chr1:300", "chr2:500", "chr2:1000"}
	if strings.Join(got, ",") != strings.Join(want, ",") {
		t.Errorf("chunked sequence = %v, want %v", got, want)
	}
}

func TestChunkedVcfReaderPadding(t *testing.T) {
	r := openTestFile(t)
	h, _ := r.Header()
	all := readAll(t, r)
	r.Close()

	dir := t.TempDir()
	// Zero-padded names must round-trip (003 -> 004 ...).
	writeChunk(t, filepath.Join(dir, "z.003.vcf.gz"), h, all[0:1])
	writeChunk(t, filepath.Join(dir, "z.004.vcf.gz"), h, all[1:2])

	c, err := NewChunkedVcfReader(filepath.Join(dir, "z.003.vcf.gz"))
	if err != nil {
		t.Fatalf("NewChunkedVcfReader: %v", err)
	}
	defer c.Close()
	n := 0
	for {
		_, err := c.NextRecord()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("NextRecord: %v", err)
		}
		n++
	}
	if n != 2 {
		t.Errorf("padded chunk read %d records, want 2", n)
	}

	if _, err := NewChunkedVcfReader(filepath.Join(dir, "unnumbered.vcf.gz")); err == nil {
		t.Errorf("expected error for an unnumbered filename")
	}
}

func TestHeaderUnionAccessors(t *testing.T) {
	const src = `##fileformat=VCFv4.2
##reference=hg38
##ALT=<ID=DEL,Description="Deletion">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1	.	A	<DEL>	.	PASS	DP=5
`
	rd, err := NewVcfReader(strings.NewReader(src))
	if err != nil {
		t.Fatalf("NewVcfReader: %v", err)
	}
	h, err := rd.Header()
	if err != nil {
		t.Fatalf("Header: %v", err)
	}
	if ids := h.AltIDs(); len(ids) != 1 || ids[0] != "DEL" {
		t.Errorf("AltIDs = %v, want [DEL]", ids)
	}
	if d, ok := h.AltDef("DEL"); !ok || d.Description != "Deletion" {
		t.Errorf("AltDef(DEL) = %+v, ok=%v", d, ok)
	}
	found := false
	for _, line := range h.OtherLines() {
		if line == "##reference=hg38" {
			found = true
		}
	}
	if !found {
		t.Errorf("OtherLines missing ##reference: %v", h.OtherLines())
	}
}
