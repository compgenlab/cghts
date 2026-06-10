package seqio

import (
	"os"
	"path/filepath"
	"testing"
)

func writeTempFai(t *testing.T, content string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), "test.fai")
	if err := os.WriteFile(p, []byte(content), 0644); err != nil {
		t.Fatal(err)
	}
	return p
}

// TestParseFaiIndexRejectsBadGeometry asserts that malformed .fai line geometry
// is rejected at parse time, before it can drive a division-by-zero (lineBases)
// or nonsensical byte offsets (lineWidth) in loadChunk.
func TestParseFaiIndexRejectsBadGeometry(t *testing.T) {
	cases := map[string]string{
		"zero lineBases":     "chr1\t100\t6\t0\t61\n",
		"zero lineWidth":     "chr1\t100\t6\t60\t0\n",
		"negative length":    "chr1\t-5\t6\t60\t61\n",
		"negative lineBases": "chr1\t100\t6\t-1\t61\n",
		"negative lineWidth": "chr1\t100\t6\t60\t-1\n",
	}
	for name, content := range cases {
		t.Run(name, func(t *testing.T) {
			if _, _, err := parseFaiIndex(writeTempFai(t, content)); err == nil {
				t.Fatalf("expected error for %q, got nil", name)
			}
		})
	}
}

// TestParseFaiIndexValid confirms a well-formed .fai still parses.
func TestParseFaiIndexValid(t *testing.T) {
	entries, names, err := parseFaiIndex(writeTempFai(t, "chr1\t100\t6\t60\t61\n"))
	if err != nil {
		t.Fatalf("valid .fai rejected: %v", err)
	}
	if len(names) != 1 || entries["chr1"] == nil || entries["chr1"].LineBases != 60 {
		t.Fatalf("unexpected parse result: names=%v entries=%v", names, entries)
	}
}
