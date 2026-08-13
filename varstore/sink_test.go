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

// memSink keeps members in memory and records what happened to them, so a test
// can assert on the sequence rather than on leftover files.
type memSink struct {
	mu       sync.Mutex
	members  map[string]*bytes.Buffer
	removed  []string
	aborted  []string
	openedIn []string
}

func newMemSink() *memSink {
	return &memSink{members: map[string]*bytes.Buffer{}}
}

func (m *memSink) Create(name string) (io.WriteCloser, error) {
	m.mu.Lock()
	defer m.mu.Unlock()
	buf := &bytes.Buffer{}
	m.members[name] = buf
	m.openedIn = append(m.openedIn, name)
	return nopCloser{buf}, nil
}

func (m *memSink) Remove(name string) error {
	m.mu.Lock()
	defer m.mu.Unlock()
	m.removed = append(m.removed, name)
	delete(m.members, name)
	return nil
}

func (m *memSink) Stat(name string) (int64, bool, error) {
	m.mu.Lock()
	defer m.mu.Unlock()
	b, ok := m.members[name]
	if !ok {
		return 0, false, nil
	}
	return int64(b.Len()), true, nil
}

func (m *memSink) Describe() string { return "memory" }

func (m *memSink) names() []string {
	m.mu.Lock()
	defer m.mu.Unlock()
	out := make([]string, 0, len(m.members))
	for n := range m.members {
		out = append(out, n)
	}
	sort.Strings(out)
	return out
}

type nopCloser struct{ io.Writer }

func (nopCloser) Close() error { return nil }

// abortingSink is a memSink that abandons rather than deletes, which is how an
// object store behaves: a member is not there until it is finished.
type abortingSink struct{ *memSink }

func (a abortingSink) Abort(name string) error {
	a.mu.Lock()
	defer a.mu.Unlock()
	a.aborted = append(a.aborted, name)
	delete(a.members, name)
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
	if last := sink.openedIn[len(sink.openedIn)-1]; last != ManifestFile {
		t.Errorf("the last member created was %s; the manifest must be written last", last)
	}
}

// A failed conversion ABANDONS rather than deletes when the sink says its
// members only appear once finished.
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
		t.Errorf("abandoned %v; all three members should have been abandoned", sink.aborted)
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
