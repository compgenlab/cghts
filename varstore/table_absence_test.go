package varstore

import (
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"testing"
)

// An optional table that is absent is an answer; one that exists but cannot be
// read is a failure. openOptionalTable used to discard every error alike, so a
// table with bad permissions, an I/O error or a symlink loop was reported as
// one that had never been written -- and the store opened with a quietly
// missing half rather than saying what was wrong.
func TestAnUnreadableTableIsNotReportedAsAbsent(t *testing.T) {
	if runtime.GOOS == "windows" {
		t.Skip("permission bits do not work this way on windows")
	}
	if os.Geteuid() == 0 {
		t.Skip("running as root; permission bits do not deny access")
	}

	base := filepath.Join(t.TempDir(), "cohort")
	buildCensusStore(t, base, WriterOpts{MinDP: 10})

	// The regions table is optional -- --no-callable stores legitimately have
	// none -- which is exactly why an unreadable one must not pass for missing.
	regions := RegionsPath(base)
	if err := os.Chmod(regions, 0o000); err != nil {
		t.Fatal(err)
	}
	defer os.Chmod(regions, 0o644)

	_, err := OpenParquet(base)
	if err == nil {
		t.Fatal("OpenParquet succeeded with an unreadable regions table; it was " +
			"treated as absent, so --hom-ref would answer from a store missing half its evidence")
	}
	if !strings.Contains(err.Error(), "regions") {
		t.Errorf("error %q should name the table that could not be read", err)
	}
}

// The genuinely-absent case still opens, since that is what --no-callable
// produces and the manifest records the zero row count.
func TestAnAbsentOptionalTableStillOpens(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	w, err := NewWriter(base, WriterOpts{Samples: []string{"S1"}, MinDP: 10, NoCallable: true})
	if err != nil {
		t.Fatal(err)
	}
	if err := w.WriteSite(Site{Chrom: "chr1", Pos: 100, Ref: "A", Alt: "T", AC: 1, AN: 2}); err != nil {
		t.Fatal(err)
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
	s, err := OpenParquet(base)
	if err != nil {
		t.Fatalf("a --no-callable store should still open: %v", err)
	}
	s.Close()
}
