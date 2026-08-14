package varstore

import (
	"bytes"
	"compress/gzip"
	"os"
	"path/filepath"
	"reflect"
	"testing"
)

func TestValidMetaKey(t *testing.T) {
	cases := []struct {
		key  string
		want bool
	}{
		{"dataset", true},
		{"reference", true},
		{"a", true},
		{"cohort_1", true},
		{"cohort-1", true},
		{"x9", true},

		{"", false},
		{"Dataset", false},  // uppercase
		{"data set", false}, // space
		{"data.set", false}, // dot -- MetaPrefix already uses it as a separator
		{"data/set", false}, // path separator
		{"café", false},     // non-ASCII
		{"cgkit.source", false},
	}
	for _, c := range cases {
		if got := ValidMetaKey(c.key); got != c.want {
			t.Errorf("ValidMetaKey(%q) = %v, want %v", c.key, got, c.want)
		}
	}
}

// Every reserved key must itself be a valid key, or a caller ranging over
// ReservedMetaKeys to build flags would generate one the writer then rejects.
func TestReservedMetaKeysAreValid(t *testing.T) {
	if len(ReservedMetaKeys) == 0 {
		t.Fatal("ReservedMetaKeys is empty")
	}
	seen := map[string]bool{}
	for _, k := range ReservedMetaKeys {
		if !ValidMetaKey(k) {
			t.Errorf("reserved key %q is not a valid metadata key", k)
		}
		if seen[k] {
			t.Errorf("reserved key %q listed twice", k)
		}
		seen[k] = true
	}
	// The constants are what callers reference; the slice is what they range
	// over. Pin that the slice actually contains them.
	for _, k := range []string{
		MetaKeyDataset, MetaKeyReference, MetaKeyCaller,
		MetaKeyAccession, MetaKeyURL, MetaKeyVersion, MetaKeyDescription,
	} {
		if !seen[k] {
			t.Errorf("reserved constant %q missing from ReservedMetaKeys", k)
		}
	}
}

// An invalid key must stop the conversion before anything on disk is touched.
// NewWriter clears the previous run's manifest as one of its first acts, so a
// key rejected later would have unmade a readable store in order to refuse the
// run meant to replace it.
func TestNewWriterRejectsBadMetaKeyBeforeTouchingDisk(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	buildCensusStore(t, base, WriterOpts{MinDP: 10})

	before, err := ReadVolumeManifest(base)
	if err != nil {
		t.Fatal(err)
	}

	w, err := NewWriter(base, WriterOpts{
		Samples: []string{"S1"},
		Meta:    map[string]string{"Dataset": "x", "bad key": "y", "ok": "z"},
	})
	if err == nil {
		w.Discard()
		t.Fatal("NewWriter accepted invalid metadata keys")
	}
	// Both offenders named, not just the first.
	for _, want := range []string{`"Dataset"`, `"bad key"`} {
		if !bytes.Contains([]byte(err.Error()), []byte(want)) {
			t.Errorf("error %q does not name %s", err, want)
		}
	}

	after, err := ReadVolumeManifest(base)
	if err != nil {
		t.Fatalf("the existing store was damaged by a rejected conversion: %v", err)
	}
	if !reflect.DeepEqual(before, after) {
		t.Error("manifest changed despite the conversion being rejected")
	}
}

func TestMetaRoundTripsThroughManifestAndStore(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	meta := map[string]string{
		MetaKeyDataset:   "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05",
		MetaKeyReference: "GRCh38",
		MetaKeyCaller:    "GATK 4.2.6.1",
		"cohort":         "1kg-phase3", // an unreserved key must survive too
	}
	buildCensusStore(t, base, WriterOpts{MinDP: 10, Meta: meta})

	m, err := ReadVolumeManifest(base)
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(m.Meta, meta) {
		t.Errorf("manifest meta = %v, want %v", m.Meta, meta)
	}

	// And independently off the calls file, which is the copy a query-time
	// caller reads. The two must agree; they are written from one source but
	// travel by different routes.
	s, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	if got := s.Provenance().Meta; !reflect.DeepEqual(got, meta) {
		t.Errorf("Provenance().Meta = %v, want %v", got, meta)
	}

	// Mutating the caller's map after the fact must not reach the manifest.
	meta["dataset"] = "mutated"
	if m.Meta[MetaKeyDataset] == "mutated" {
		t.Error("manifest aliases the caller's map")
	}
}

// A store converted without metadata must read back identically to one written
// before the field existed: absent, not an empty object. Otherwise "not stated"
// and "stated as nothing" become indistinguishable.
func TestNoMetaIsOmittedEntirely(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	buildCensusStore(t, base, WriterOpts{MinDP: 10})

	m, err := ReadVolumeManifest(base)
	if err != nil {
		t.Fatal(err)
	}
	if m.Meta != nil {
		t.Errorf("Meta = %v, want nil", m.Meta)
	}

	f, err := os.Open(VolumeManifestPath(base))
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	zr, err := gzip.NewReader(f)
	if err != nil {
		t.Fatal(err)
	}
	defer zr.Close()
	var buf bytes.Buffer
	if _, err := buf.ReadFrom(zr); err != nil {
		t.Fatal(err)
	}
	if bytes.Contains(buf.Bytes(), []byte(`"meta"`)) {
		t.Errorf("manifest carries a meta key when none was recorded:\n%s", buf.String())
	}

	s, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	if got := s.Provenance().Meta; got != nil {
		t.Errorf("Provenance().Meta = %v, want nil", got)
	}
}

// Caller metadata is namespaced, so it cannot displace provenance the writer
// vouches for. A key spelled "source" must not become cgkit.source.
func TestMetaCannotShadowProvenanceKeys(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	buildCensusStore(t, base, WriterOpts{
		MinDP:   10,
		Program: "cgkit test",
		Sources: []string{"real.vcf"},
		Meta:    map[string]string{"source": "spoofed", "program": "spoofed"},
	})

	s, err := OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()

	p := s.Provenance()
	if p.Source != "real.vcf" {
		t.Errorf("Source = %q, want %q", p.Source, "real.vcf")
	}
	if p.Program != "cgkit test" {
		t.Errorf("Program = %q, want %q", p.Program, "cgkit test")
	}
	if p.Meta["source"] != "spoofed" {
		t.Errorf("caller metadata lost: %v", p.Meta)
	}
}
