package varstore

import (
	"bytes"
	"fmt"
	"io"
	"sort"
	"strings"
	"sync"
	"testing"
)

// A store written through a Sink that is not a filesystem.
//
// This is what the interface buys beyond remote output: the writer can be
// exercised with no directory, no disk and no cleanup, which the previous
// version could not be at all. Every writer test had to be a filesystem test.

// memSink keeps tables in memory and records what happened to them, so a test
// can assert on the sequence rather than on leftover files.
type memSink struct {
	mu       sync.Mutex
	tables   map[string]*bytes.Buffer
	removed  []string
	aborted  []string
	openedIn []string
}

func newMemSink() *memSink {
	return &memSink{tables: map[string]*bytes.Buffer{}}
}

func (m *memSink) Create(name string) (io.WriteCloser, error) {
	m.mu.Lock()
	defer m.mu.Unlock()
	buf := &bytes.Buffer{}
	m.tables[name] = buf
	m.openedIn = append(m.openedIn, name)
	return nopCloser{buf}, nil
}

func (m *memSink) Remove(name string) error {
	m.mu.Lock()
	defer m.mu.Unlock()
	m.removed = append(m.removed, name)
	delete(m.tables, name)
	return nil
}

func (m *memSink) Stat(name string) (int64, bool, error) {
	m.mu.Lock()
	defer m.mu.Unlock()
	b, ok := m.tables[name]
	if !ok {
		return 0, false, nil
	}
	return int64(b.Len()), true, nil
}

func (m *memSink) Describe() string { return "memory" }

func (m *memSink) names() []string {
	m.mu.Lock()
	defer m.mu.Unlock()
	out := make([]string, 0, len(m.tables))
	for n := range m.tables {
		out = append(out, n)
	}
	sort.Strings(out)
	return out
}

type nopCloser struct{ io.Writer }

func (nopCloser) Close() error { return nil }

// abortingSink is a memSink that abandons rather than deletes, which is how an
// object store behaves: a table is not there until it is finished.
type abortingSink struct{ *memSink }

func (a abortingSink) Abort(name string) error {
	a.mu.Lock()
	defer a.mu.Unlock()
	a.aborted = append(a.aborted, name)
	delete(a.tables, name)
	return nil
}

func TestAStoreCanBeWrittenWithoutAFilesystem(t *testing.T) {
	sink := newMemSink()
	w, err := NewWriter("memory://store", WriterOpts{
		Sink: sink, Samples: []string{"S1", "S2"}, MinDP: 10, Program: "test",
	})
	if err != nil {
		t.Fatal(err)
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}

	want := []string{"calls.parquet", "manifest.json.gz", "regions.parquet", "sites.parquet"}
	if got := sink.names(); !equal(got, want) {
		t.Errorf("wrote %v, want %v", got, want)
	}

	// The manifest is written LAST and is what makes a store readable. A
	// conversion that stopped earlier must not leave one, so its position in
	// the sequence is part of the contract rather than an accident of ordering.
	if last := sink.openedIn[len(sink.openedIn)-1]; last != VolumeManifestFile {
		t.Errorf("the last table created was %s; the manifest must be written last", last)
	}
}

// A failed conversion ABANDONS rather than deletes when the sink says its
// tables only appear once finished.
//
// The distinction is invisible locally and load-bearing remotely: calling
// Remove on an upload that never completed deletes nothing and leaves its parts
// behind, invisible to a listing and billed.
func TestAbortAbandonsOnASinkThatCannotDelete(t *testing.T) {
	sink := abortingSink{newMemSink()}
	w, err := NewWriter("memory://store", WriterOpts{
		Sink: sink, Samples: []string{"S1"}, Program: "test",
	})
	if err != nil {
		t.Fatal(err)
	}
	if err := w.Discard(); err != nil {
		t.Fatal(err)
	}

	if len(sink.aborted) != 3 {
		t.Errorf("abandoned %v; all three tables should have been abandoned", sink.aborted)
	}
	for _, name := range sink.removed {
		if strings.HasSuffix(name, ".parquet") {
			t.Errorf("%s was deleted rather than abandoned", name)
		}
	}
	if left := sink.names(); len(left) != 0 {
		t.Errorf("%v survived a discarded conversion", left)
	}
}

// An unregistered scheme must be reported as a missing transport, not treated
// as a very strangely named directory.
func TestAnUnknownSchemeIsRefusedByName(t *testing.T) {
	_, err := OpenSink("gs://bucket/store")
	if err == nil {
		t.Fatal("gs:// was accepted with no transport registered for it")
	}
	if !strings.Contains(err.Error(), "gs") {
		t.Errorf("the error does not name the scheme: %v", err)
	}
	if CanWrite("gs://bucket/store") {
		t.Error("CanWrite says gs:// is writable")
	}
	if !CanWrite("/tmp/store") {
		t.Error("CanWrite says a local path is not writable")
	}
}

func equal(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

var _ = fmt.Sprintf

// The overwrite guard has to mean the same thing wherever the store is going.
//
// It used to stat local paths, so against a remote base it found nothing and
// waved through exactly the overwrite it exists to prevent -- silently, because
// the tables it looks for cannot be seen with os.Stat. This is that case: a
// sink holding a store, asked whether it may be written over.
func TestTheOverwriteGuardSeesTablesInAnySink(t *testing.T) {
	sink := newMemSink()

	// Nothing there yet.
	if err := CheckStoreTargetIn(sink, false); err != nil {
		t.Fatalf("refused an empty destination: %v", err)
	}

	// Write a store into it, then ask again.
	w, err := NewWriter("memory://store", WriterOpts{Sink: sink, Samples: []string{"S1"}})
	if err != nil {
		t.Fatal(err)
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}

	err = CheckStoreTargetIn(sink, false)
	if err == nil {
		t.Fatal("a populated destination was accepted; this is the silent-overwrite case")
	}
	// The reader has to be able to go and look at what stopped them.
	if !strings.Contains(err.Error(), "calls.parquet") {
		t.Errorf("the refusal does not name a table: %v", err)
	}
	if CheckStoreTargetIn(sink, true) != nil {
		t.Error("--force did not override the guard")
	}
}

// A partial set stops it too: the tables are only meaningful together, so a
// half-replaced store is worse than either keeping or replacing the old one.
func TestASingleSurvivingTableIsEnoughToStop(t *testing.T) {
	sink := newMemSink()
	w, err := NewWriter("memory://store", WriterOpts{Sink: sink, Samples: []string{"S1"}})
	if err != nil {
		t.Fatal(err)
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
	for _, n := range []string{"calls.parquet", "regions.parquet", VolumeManifestFile} {
		if err := sink.Remove(n); err != nil {
			t.Fatal(err)
		}
	}

	err = CheckStoreTargetIn(sink, false)
	if err == nil {
		t.Fatal("a lone surviving table was not enough to stop an overwrite")
	}
	if !strings.Contains(err.Error(), "sites.parquet") {
		t.Errorf("the refusal does not name the survivor: %v", err)
	}
}

// EnsureStoreDir is for local paths, and says so rather than making a mess.
//
// On an s3:// locator it used to reach MkdirAll and create a local directory
// literally named "s3:", which is a confusing thing to find and says nothing
// about why the store was not written.
func TestEnsureStoreDirRefusesARemoteBase(t *testing.T) {
	err := EnsureStoreDir("s3://bucket/cohort")
	if err == nil {
		t.Fatal("EnsureStoreDir accepted an s3:// base; it would create a local `s3:` directory")
	}
	if !strings.Contains(err.Error(), "s3") {
		t.Errorf("the error does not name the scheme: %v", err)
	}
	if err := EnsureStoreDir(t.TempDir() + "/store"); err != nil {
		t.Errorf("EnsureStoreDir refused a local path: %v", err)
	}
}
