package varstore

import (
	"os"
	"path/filepath"
	"testing"
)

// Every spelling of a store has to land on the same directory, and none of it
// may depend on the filesystem: a remote locator cannot be stat'd, and the
// answer must not quietly change based on what happens to exist.
func TestTrimStoreSuffixSpellings(t *testing.T) {
	for _, tc := range []struct{ in, want string }{
		{"cohort", "cohort"},
		{"cohort/", "cohort"},
		{"cohort//", "cohort"},
		{"cohort/calls.parquet", "cohort"},
		{"cohort/sites.parquet", "cohort"},
		{"cohort/regions.parquet", "cohort"},
		{"/data/stores/cohort", "/data/stores/cohort"},
		{"/data/stores/cohort/", "/data/stores/cohort"},
		{"/data/stores/cohort/calls.parquet", "/data/stores/cohort"},

		// A table named on its own is a store in the current directory.
		{"calls.parquet", ""},

		// Remote locators resolve identically, and the scheme's "//" survives.
		{"s3://bucket/cohort", "s3://bucket/cohort"},
		{"s3://bucket/cohort/", "s3://bucket/cohort"},
		{"s3://bucket/cohort/calls.parquet", "s3://bucket/cohort"},
		{"https://host/d/cohort/regions.parquet", "https://host/d/cohort"},

		// The filesystem root is a directory name, not a decorated empty one.
		{"/", "/"},
	} {
		if got := TrimStoreSuffix(tc.in); got != tc.want {
			t.Errorf("TrimStoreSuffix(%q) = %q, want %q", tc.in, got, tc.want)
		}
	}
}

func TestTablePathKeepsTheScheme(t *testing.T) {
	for _, tc := range []struct{ base, want string }{
		{"cohort", "cohort/calls.parquet"},
		{"cohort/", "cohort/calls.parquet"},
		{"/data/cohort", "/data/cohort/calls.parquet"},
		{"", "calls.parquet"},
		// filepath.Join would clean this to "s3:/bucket/cohort/calls.parquet".
		{"s3://bucket/cohort", "s3://bucket/cohort/calls.parquet"},
		{"https://host/cohort/", "https://host/cohort/calls.parquet"},
	} {
		if got := CallsPath(tc.base); got != tc.want {
			t.Errorf("CallsPath(%q) = %q, want %q", tc.base, got, tc.want)
		}
	}
}

// Reading a store back must work whichever way it is named, since a shell
// completion lands on a table and a script is likely to carry the slash.
func TestStoreOpensUnderEverySpelling(t *testing.T) {
	dir := t.TempDir()
	base := filepath.Join(dir, "cohort")
	w := writeSomeRows(t, base)
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}

	want := -1
	for _, spelling := range []string{base, base + "/", CallsPath(base), SitesPath(base)} {
		s, err := OpenParquet(spelling)
		if err != nil {
			t.Fatalf("OpenParquet(%q): %v", spelling, err)
		}
		got, err := CollectCalls(s, Query{})
		s.Close()
		if err != nil {
			t.Fatalf("CollectCalls via %q: %v", spelling, err)
		}
		if want < 0 {
			want = len(got)
			if want == 0 {
				t.Fatal("no calls read back; the fixture is not exercising anything")
			}
			continue
		}
		if len(got) != want {
			t.Errorf("via %q got %d calls, want %d", spelling, len(got), want)
		}
	}
}

// The overwrite guard keys on the tables, so a directory holding unrelated
// files is a fine target and its contents are left alone.
func TestCheckStoreTargetIgnoresUnrelatedFiles(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	if err := EnsureStoreDir(base); err != nil {
		t.Fatal(err)
	}
	if err := writeFile(filepath.Join(base, "README.txt"), "not a table"); err != nil {
		t.Fatal(err)
	}
	if err := CheckStoreTarget(base, false); err != nil {
		t.Errorf("refused a directory holding no tables: %v", err)
	}

	w := writeSomeRows(t, base)
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
	if err := CheckStoreTarget(base, false); err == nil {
		t.Error("allowed overwriting a populated store without --force")
	}
	if err := CheckStoreTarget(base, true); err != nil {
		t.Errorf("--force still refused: %v", err)
	}
}

func writeFile(path, content string) error {
	return os.WriteFile(path, []byte(content), 0o644)
}
