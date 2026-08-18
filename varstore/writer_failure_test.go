package varstore

import (
	"errors"
	"os"
	"path/filepath"
	"testing"

	"github.com/parquet-go/parquet-go"
)

// writeSomeRows builds a writer at base and puts a few rows in every table,
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

func TestDiscardRemovesEveryTable(t *testing.T) {
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
// footer into each table, so a process killed partway through Discard -- or an
// unlink that simply fails -- would leave three well-formed files holding only
// part of the input. That is the one partial store a reader cannot detect.
func TestDiscardDoesNotFinalizeTables(t *testing.T) {
	if os.Geteuid() == 0 {
		t.Skip("running as root: directory permissions would not prevent the unlink")
	}
	base := filepath.Join(t.TempDir(), "cohort")
	w := writeSomeRows(t, base)

	// Make the unlink fail, so the files survive Discard and can be inspected.
	// The store directory itself has to be the read-only one: the tables live
	// inside it, and removing a file needs write permission on its directory.
	if err := os.Chmod(base, 0o555); err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { os.Chmod(base, 0o755) })

	err := w.Discard()
	if err == nil {
		t.Fatal("Discard reported success though the tables could not be removed")
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
				"table must never look like a finished one", filepath.Base(p))
		}
	}
}

// Close must not finish the other tables once one has failed. Three valid
// footers where one file is short is undetectable; a failure that stops is not.
func TestCloseStopsAtTheFirstFailure(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	w := writeSomeRows(t, base)

	// Close the calls table out from under the writer so its flush fails.
	if err := w.tables[0].sw.Close(); err != nil {
		t.Fatal(err)
	}

	if err := w.Close(); err == nil {
		t.Fatal("Close succeeded though the calls table could not be written")
	}

	// The tables that were never finalized must not carry a footer.
	for _, p := range []string{SitesPath(base), RegionsPath(base)} {
		f, oerr := os.Open(p)
		if oerr != nil {
			t.Fatal(oerr)
		}
		st, _ := f.Stat()
		_, perr := parquet.OpenFile(f, st.Size())
		f.Close()
		if perr == nil {
			t.Errorf("%s was finalized after an earlier table failed; the set "+
				"would look complete while calls is short", filepath.Base(p))
		}
	}
}

// A table that is present but unreadable must not be mistaken for one that was
// deliberately omitted. Before the footer check, a truncated sites file and a
// --no-callable store looked identical at open, and the store answered queries
// as though the runs had never been requested.
func TestATruncatedTableIsNotReadAsAbsent(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	w := writeSomeRows(t, base)
	if err := w.Finish(); err != nil {
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
		t.Fatal("a truncated sites table opened as though it were simply absent")
	}
}

// A table that was never populated may go missing without complaint -- that is
// the --no-callable shape -- but one the manifest vouched for may not. Before
// the manifest the two were the same event, so deleting a regions file full of
// runs read as "this store never tracked callable regions" and --hom-ref
// silently changed its answer.
func TestAbsentMemberIsOnlyAllowedWhenItHeldNothing(t *testing.T) {
	t.Run("empty table may vanish", func(t *testing.T) {
		base := filepath.Join(t.TempDir(), "cohort")
		w, err := NewWriter(base, WriterOpts{
			Samples: []string{"S1"}, MinDP: 10, RowGroupSize: 1, NoCallable: true,
		})
		if err != nil {
			t.Fatal(err)
		}
		if err := w.WriteCall(Call{
			SampleID: "S1", Chrom: "chr1", Pos: 100, Ref: "A", Alt: "T", GT: "0/1", DP: 30,
		}); err != nil {
			t.Fatal(err)
		}
		if err := w.WriteSite(Site{
			Chrom: "chr1", Pos: 100, Ref: "A", Alt: "T", AC: 1, AN: 2, NCarriers: 1,
		}); err != nil {
			t.Fatal(err)
		}
		if err := w.Finish(); err != nil {
			t.Fatal(err)
		}
		if err := os.Remove(RegionsPath(base)); err != nil {
			t.Fatal(err)
		}
		s, err := OpenParquet(base)
		if err != nil {
			t.Fatalf("a store whose regions table held nothing should open: %v", err)
		}
		defer s.Close()
		if s.hasRegions {
			t.Error("hasRegions is true though the table was removed")
		}
	})

	t.Run("populated table may not", func(t *testing.T) {
		base := filepath.Join(t.TempDir(), "cohort")
		w := writeSomeRows(t, base)
		if err := w.Finish(); err != nil {
			t.Fatal(err)
		}
		if err := os.Remove(RegionsPath(base)); err != nil {
			t.Fatal(err)
		}
		if _, err := OpenParquet(base); err == nil {
			t.Fatal("a store opened with its regions table deleted, though the " +
				"manifest recorded runs in it")
		}
	})
}

// Site must refuse like its siblings rather than dereference a nil table. An
// absent sites file is a reachable state -- the manifest permits one that
// recorded no rows -- and Site was the only accessor that did not check.
func TestSiteRefusesWithoutTheSitesMember(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	w, err := NewWriter(base, WriterOpts{Samples: []string{"S1"}, MinDP: 10, NoCallable: true})
	if err != nil {
		t.Fatal(err)
	}
	if err := w.WriteCall(Call{
		SampleID: "S1", Chrom: "chr1", Pos: 100, Ref: "A", Alt: "T", GT: "0/1", DP: 30,
	}); err != nil {
		t.Fatal(err)
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
	if err := os.Remove(SitesPath(base)); err != nil {
		t.Fatal(err)
	}

	s, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	if _, _, err := s.Site(Locus{Chrom: "chr1", Pos: 100, Ref: "A", Alt: "T"}); err == nil {
		t.Error("Site succeeded with no sites table")
	} else if !errors.Is(err, ErrNotClassifiable) {
		t.Errorf("Site error is not ErrNotClassifiable: %v", err)
	}
}
