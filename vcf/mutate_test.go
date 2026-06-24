package vcf

import (
	"bytes"
	"testing"
)

// serializeLine builds a record from a raw line, applies fn, and returns the
// writer output (which uses the verbatim path when clean, reconstruct when dirty).
func writeRec(t *testing.T, rec *VcfRecord) string {
	t.Helper()
	var buf bytes.Buffer
	w := NewVcfWriter(&buf)
	if err := w.WriteRecord(rec); err != nil {
		t.Fatalf("WriteRecord: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	return buf.String()
}

func TestWriteVerbatimWhenClean(t *testing.T) {
	const line = "chr1\t100\trs1\tA\tG\t50.0\tPASS\tDP=30;AF=0.5;DB\tGT:AD\t0/0:28,2\t0/1:15,15"
	rec, err := newRecord(line, nil)
	if err != nil {
		t.Fatal(err)
	}
	if got := writeRec(t, rec); got != line+"\n" {
		t.Errorf("clean record not verbatim.\n got: %q\nwant: %q", got, line+"\n")
	}
}

func TestAddInfoFlagSerialize(t *testing.T) {
	const line = "chr1\t100\trs1\tA\tG\t50.0\tPASS\tDP=30\tGT\t0/0\t0/1"
	rec, _ := newRecord(line, nil)
	rec.AddInfoFlag("CG_INSERT")
	rec.AddInfo("CG_INSLEN", "3")
	if !rec.Dirty() {
		t.Fatal("record should be dirty after AddInfo")
	}
	want := "chr1\t100\trs1\tA\tG\t50\tPASS\tDP=30;CG_INSERT;CG_INSLEN=3\tGT\t0/0\t0/1\n"
	if got := writeRec(t, rec); got != want {
		t.Errorf("serialize mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestSetIDSerialize(t *testing.T) {
	const line = "chr1\t100\t.\tA\tG\t.\t.\t.\t"
	// Trailing tab makes 8 columns (no samples); INFO=".".
	rec, err := newRecord("chr1\t100\t.\tA\tG\t.\t.\t.", nil)
	if err != nil {
		t.Fatal(err)
	}
	rec.SetID("chr1_100_A_G")
	// FILTER was "." on input, so it reconstructs as "." (not PASS).
	want := "chr1\t100\tchr1_100_A_G\tA\tG\t.\t.\t.\n"
	if got := writeRec(t, rec); got != want {
		t.Errorf("SetID serialize mismatch.\n got: %q\nwant: %q", got, want)
	}
	_ = line
}

func TestAddFormatAllSamples(t *testing.T) {
	const line = "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/0:28,2\t0/1:15,15"
	rec, _ := newRecord(line, nil)
	// Add CG_DS to every sample (like the Dosage annotator).
	for i := 0; i < rec.NumSamples(); i++ {
		if err := rec.AddFormat(i, "CG_DS", "0"); err != nil {
			t.Fatal(err)
		}
	}
	// FORMAT keys come from sample 0: GT:AD:CG_DS.
	want := "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT:AD:CG_DS\t0/0:28,2:0\t0/1:15,15:0\n"
	if got := writeRec(t, rec); got != want {
		t.Errorf("AddFormat serialize mismatch.\n got: %q\nwant: %q", got, want)
	}
}

func TestSerializeFilterMissingVsPass(t *testing.T) {
	// FILTER "." should reconstruct as "." (not PASS); a real code stays.
	dot, _ := newRecord("chr1\t100\t.\tA\tG\t.\t.\t.", nil)
	dot.AddInfoFlag("X")
	if got := writeRec(t, dot); got != "chr1\t100\t.\tA\tG\t.\t.\tX\n" {
		t.Errorf("FILTER '.' mismatch: %q", got)
	}
	filt, _ := newRecord("chr1\t100\t.\tA\tG\t.\tlowqual\t.", nil)
	filt.AddInfoFlag("X")
	if got := writeRec(t, filt); got != "chr1\t100\t.\tA\tG\t.\tlowqual\tX\n" {
		t.Errorf("FILTER code mismatch: %q", got)
	}
}

func TestFormatTrailingMissingTrim(t *testing.T) {
	// A short sample (missing trailing AD) reconstructs with the trailing
	// missing trimmed; GT is always kept.
	const line = "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/0:28,2\t0/0"
	rec, _ := newRecord(line, nil)
	rec.AddInfoFlag("X") // make dirty so samples reconstruct
	want := "chr1\t100\t.\tA\tG\t.\tPASS\tX\tGT:AD\t0/0:28,2\t0/0\n"
	if got := writeRec(t, rec); got != want {
		t.Errorf("trailing-trim mismatch.\n got: %q\nwant: %q", got, want)
	}
}
