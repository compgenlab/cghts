package varstore

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/parquet-go/parquet-go"
)

// writeSomeRows builds a writer at base and puts a few rows in every member,
// without closing it. The caller decides how the conversion ends.
func writeSomeRows(t *testing.T, base string) *Writer {
	t.Helper()
	w, err := NewWriter(base, WriterOpts{Samples: []string{"S1"}, MinDP: 10, RowGroupSize: 1})
	if err != nil {
		t.Fatal(err)
	}
	for pos := int32(100); pos < 110; pos++ {
		if err := w.WriteCall(Call{
			SampleID: "S1", Chrom: "chr1", Pos: pos, Ref: "A", Alt: "T", GT: "0/1", DP: 30,
		}); err != nil {
			t.Fatal(err)
		}
		if err := w.WriteSite(Site{
			Chrom: "chr1", Pos: pos, Ref: "A", Alt: "T", AC: 1, AN: 2, NCarriers: 1, NCalled: 1,
		}); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.WriteRegion(CalledSiteRun{
		SampleID: "S1", Chrom: "chr1", Start: 100, End: 109, NSites: 10,
	}); err != nil {
		t.Fatal(err)
	}
	return w
}

func TestDiscardRemovesEveryMember(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	w := writeSomeRows(t, base)
	if err := w.Discard(); err != nil {
		t.Fatalf("Discard: %v", err)
	}
	for _, p := range []string{CallsPath(base), SitesPath(base), RegionsPath(base)} {
		if _, err := os.Stat(p); !os.IsNotExist(err) {
			t.Errorf("%s still present after Discard (err=%v)", p, err)
		}
	}
}

// A discarded conversion must not leave readable parquet behind, even when the
// unlink cannot happen. Closing the parquet writers first would write a complete
// footer into each member, so a process killed partway through Discard -- or an
// unlink that simply fails -- would leave three well-formed files holding only
// part of the input. That is the one partial store a reader cannot detect.
func TestDiscardDoesNotFinalizeMembers(t *testing.T) {
	if os.Geteuid() == 0 {
		t.Skip("running as root: directory permissions would not prevent the unlink")
	}
	dir := t.TempDir()
	base := filepath.Join(dir, "cohort")
	w := writeSomeRows(t, base)

	// Make the unlink fail, so the files survive Discard and can be inspected.
	if err := os.Chmod(dir, 0o555); err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { os.Chmod(dir, 0o755) })

	err := w.Discard()
	if err == nil {
		t.Fatal("Discard reported success though the members could not be removed")
	}

	for _, p := range []string{CallsPath(base), SitesPath(base), RegionsPath(base)} {
		f, oerr := os.Open(p)
		if oerr != nil {
			continue // removed after all; nothing to check
		}
		st, serr := f.Stat()
		if serr != nil {
			f.Close()
			t.Fatal(serr)
		}
		_, perr := parquet.OpenFile(f, st.Size())
		f.Close()
		if perr == nil {
			t.Errorf("%s has a valid parquet footer after Discard; a discarded "+
				"member must never look like a finished one", filepath.Base(p))
		}
	}
}

// Close must not finish the other members once one has failed. Three valid
// footers where one file is short is undetectable; a failure that stops is not.
func TestCloseStopsAtTheFirstFailure(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	w := writeSomeRows(t, base)

	// Close the calls handle out from under the writer so its flush fails.
	if err := w.files[0].Close(); err != nil {
		t.Fatal(err)
	}

	if err := w.Close(); err == nil {
		t.Fatal("Close succeeded though the calls member could not be written")
	}

	// The members that were never finalized must not carry a footer.
	for _, p := range []string{SitesPath(base), RegionsPath(base)} {
		f, oerr := os.Open(p)
		if oerr != nil {
			t.Fatal(oerr)
		}
		st, _ := f.Stat()
		_, perr := parquet.OpenFile(f, st.Size())
		f.Close()
		if perr == nil {
			t.Errorf("%s was finalized after an earlier member failed; the set "+
				"would look complete while calls is short", filepath.Base(p))
		}
	}
}

// A member that is present but unreadable must not be mistaken for one that was
// deliberately omitted. Before the footer check, a truncated sites file and a
// --no-callable store looked identical at open, and the store answered queries
// as though the runs had never been requested.
func TestTruncatedMemberIsNotReadAsAbsent(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	w := writeSomeRows(t, base)
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}

	// Sanity: it opens cleanly before the damage.
	s, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	s.Close()

	if err := os.Truncate(SitesPath(base), 12); err != nil {
		t.Fatal(err)
	}
	if _, err := OpenParquet(base); err == nil {
		t.Fatal("a truncated sites member opened as though it were simply absent")
	}
}

// The absent case must still be absent: a store with no regions member is what
// --no-callable produces, and it has to keep opening.
func TestAbsentMemberStillOpens(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	w := writeSomeRows(t, base)
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	if err := os.Remove(RegionsPath(base)); err != nil {
		t.Fatal(err)
	}
	s, err := OpenParquet(base)
	if err != nil {
		t.Fatalf("a store without a regions member should open: %v", err)
	}
	defer s.Close()
	if s.hasRegions {
		t.Error("hasRegions is true though the member was removed")
	}
}
