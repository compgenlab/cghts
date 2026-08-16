package varstore

import (
	"context"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// captureStore writes one chromosome capturing an INFO field.
func captureStore(t *testing.T, dir, chrom string, info []InfoField) string {
	t.Helper()
	base := filepath.Join(dir, chrom)
	w, err := NewWriter(base, WriterOpts{
		Samples: []string{"S1", "S2"}, MinDP: 10, RowGroupSize: 8, Info: info,
	})
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i < 4; i++ {
		pos := int32(100 + i*10)
		vals := map[string]any{}
		if len(info) > 0 {
			vals["R2"] = 0.9
		}
		if err := w.WriteSiteInfo(Site{
			Chrom: chrom, Pos: pos, Ref: "G", Alt: "A", AN: 4, NCalled: 2,
		}, vals); err != nil {
			t.Fatal(err)
		}
	}
	for _, s := range []string{"S1", "S2"} {
		if err := w.WriteRegion(CalledSiteRun{
			SampleID: s, Chrom: chrom, Start: 100, End: 130, NSites: 4, MinDP: 30,
		}); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
	return base
}

// An archive answers capability questions exactly as one of its volumes does.
//
// THIS IS A SILENT FAILURE IF IT REGRESSES, which is why it is a test rather
// than an observation. Callers ask a store what it captured through interface
// assertions, because the Store interface deliberately does not carry them. An
// archive missing the methods does not error -- it fails the assertion and is
// reported as a store that captured nothing, which is indistinguishable from a
// conversion that really captured nothing. A registration reading it would
// record no INFO fields for a store holding R2, and every later question about
// how a site was arrived at would be refused for the wrong reason.
func TestAnArchiveReportsWhatItsVolumesCaptured(t *testing.T) {
	dir := t.TempDir()
	info := []InfoField{{Name: "R2", Column: "info_r2", Type: InfoFloat, Number: "1"}}
	captureStore(t, dir, "chr1", info)
	captureStore(t, dir, "chr2", info)

	man, err := BuildStore(context.Background(), dir, []string{"chr1", "chr2"}, "test", "test")
	if err != nil {
		t.Fatal(err)
	}
	if err := WriteStoreManifest(dir, *man); err != nil {
		t.Fatal(err)
	}

	vol, err := OpenParquet(filepath.Join(dir, "chr1"))
	if err != nil {
		t.Fatal(err)
	}
	defer vol.Close()

	arc, err := OpenStore(context.Background(), dir)
	if err != nil {
		t.Fatal(err)
	}
	defer arc.Close()

	// The assertions the adapter actually makes, made here against both shapes.
	vi, ok := any(vol).(interface{ InfoFields() []InfoField })
	if !ok {
		t.Fatal("a volume does not satisfy the InfoFields assertion")
	}
	ai, ok := any(arc).(interface{ InfoFields() []InfoField })
	if !ok {
		t.Fatal("an archive does not satisfy the InfoFields assertion, so it reports capturing nothing")
	}
	if len(ai.InfoFields()) == 0 {
		t.Fatal("the archive reports no captured INFO fields for a store that captured R2")
	}
	if !sameInfo(ai.InfoFields(), vi.InfoFields()) {
		t.Errorf("archive captured %v, volume captured %v", ai.InfoFields(), vi.InfoFields())
	}

	if _, ok := any(arc).(interface {
		SitesInfo(func(Site, InfoRow) bool) error
	}); !ok {
		t.Fatal("an archive does not satisfy the SitesInfo assertion, so no site evidence would land")
	}

	// And the walk covers every volume, not just the first.
	n := 0
	if err := arc.SitesInfo(func(Site, InfoRow) bool { n++; return true }); err != nil {
		t.Fatal(err)
	}
	if n != 8 {
		t.Errorf("the archive-wide provenance walk saw %d sites, want 8 across two volumes", n)
	}

	if arc.SpanSemantics() != vol.SpanSemantics() {
		t.Errorf("archive spans %q, volume spans %q", arc.SpanSemantics(), vol.SpanSemantics())
	}
	if arc.HasCoverage() != vol.HasCoverage() {
		t.Errorf("archive coverage %v, volume coverage %v", arc.HasCoverage(), vol.HasCoverage())
	}
}

// A volume whose capture disagrees with the archive is refused when it opens.
//
// The archive answers for ONE set of captured fields. A volume holding another
// would report a site as unmeasured only because nobody asked it that question,
// which reads exactly like a site that was measured and failed.
func TestAnArchiveRefusesAVolumeThatCapturedSomethingElse(t *testing.T) {
	dir := t.TempDir()
	withR2 := []InfoField{{Name: "R2", Column: "info_r2", Type: InfoFloat, Number: "1"}}
	captureStore(t, dir, "chr1", withR2)
	captureStore(t, dir, "chr2", withR2)

	man, err := BuildStore(context.Background(), dir, []string{"chr1", "chr2"}, "test", "test")
	if err != nil {
		t.Fatal(err)
	}
	if err := WriteStoreManifest(dir, *man); err != nil {
		t.Fatal(err)
	}

	// Replace chr2 with one that captured nothing, leaving the archive's claim
	// in place -- which is what copying a volume in from another conversion
	// looks like from the outside.
	other := t.TempDir()
	plain := captureStore(t, other, "chr2", nil)
	victim := filepath.Join(dir, "chr2")
	if err := os.RemoveAll(victim); err != nil {
		t.Fatal(err)
	}
	if err := os.Rename(plain, victim); err != nil {
		t.Fatal(err)
	}

	arc, err := OpenStore(context.Background(), dir)
	if err != nil {
		t.Fatal(err)
	}
	defer arc.Close()

	// Lazy: the disagreement surfaces when chr2 is first read.
	_, err = arc.Classify(Locus{Chrom: "chr2", Pos: 100, Ref: "G", Alt: "A"}, Gate{MinDP: 10})
	if err == nil {
		t.Fatal("a volume capturing different INFO fields was accepted")
	}
	if !strings.Contains(err.Error(), "captured INFO") {
		t.Errorf("the error does not name the disagreement: %v", err)
	}
}
