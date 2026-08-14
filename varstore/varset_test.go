package varstore

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"testing"
)

// A varset: several stores, disjoint by chromosome, read as one.
//
// The tests that matter are that it answers exactly as its members do, and that
// it REFUSES to be built from members that disagree -- because a set whose
// members disagree answers with a different population per chromosome, and
// nothing outside it could see that.

// chromStore writes a one-chromosome store with three samples.
func chromStore(t *testing.T, dir, chrom string, samples []string, minDP int32) string {
	t.Helper()
	base := filepath.Join(dir, chrom)
	w, err := NewWriter(base, WriterOpts{Samples: samples, MinDP: minDP})
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i < 10; i++ {
		pos := int32(100 + i*10)
		if err := w.WriteSite(Site{
			Chrom: chrom, Pos: pos, Ref: "G", Alt: "A", AN: int32(2 * len(samples)), NCalled: int32(len(samples)),
		}); err != nil {
			t.Fatal(err)
		}
		// The first sample carries at every third site.
		if i%3 == 0 {
			if err := w.WriteCall(Call{
				SampleID: samples[0], Chrom: chrom, Pos: pos, Ref: "G", Alt: "A",
				GT: "0/1", DP: 40, ADRef: 20, ADAlt: 20, GQ: Missing,
			}); err != nil {
				t.Fatal(err)
			}
		}
	}
	for _, name := range samples {
		if err := w.WriteRegion(CalledSiteRun{
			SampleID: name, Chrom: chrom, Start: 100, End: 190, NSites: 10, MinDP: 30,
		}); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Finish(); err != nil {
		t.Fatal(err)
	}
	return base
}

func buildSet(t *testing.T, chroms []string, samples []string) string {
	t.Helper()
	dir := t.TempDir()
	var members []string
	for _, c := range chroms {
		chromStore(t, dir, c, samples, 10)
		members = append(members, c)
	}
	man, err := BuildSet(context.Background(), dir, members, "test", "test")
	if err != nil {
		t.Fatalf("BuildSet: %v", err)
	}
	if err := WriteSetManifest(dir, *man); err != nil {
		t.Fatal(err)
	}
	return dir
}

func TestSetAnswersAsItsMembersDo(t *testing.T) {
	samples := []string{"S1", "S2", "S3"}
	chroms := []string{"chr1", "chr2", "chr7"}
	dir := buildSet(t, chroms, samples)

	set, err := OpenSet(context.Background(), dir)
	if err != nil {
		t.Fatal(err)
	}
	defer set.Close()

	gate := Gate{MinDP: 10}
	for _, chrom := range chroms {
		one, err := OpenParquet(filepath.Join(dir, chrom))
		if err != nil {
			t.Fatal(err)
		}
		for i := 0; i < 10; i++ {
			l := Locus{Chrom: chrom, Pos: int32(100 + i*10), Ref: "G", Alt: "A"}
			a, err := one.Classify(l, gate)
			if err != nil {
				t.Fatal(err)
			}
			b, err := set.Classify(l, gate)
			if err != nil {
				t.Fatal(err)
			}
			by := map[string]State{}
			for _, st := range b {
				by[st.SampleID] = st.State
			}
			for _, st := range a {
				if by[st.SampleID] != st.State {
					t.Errorf("%s %s: member %q, set %q", l, st.SampleID, st.State, by[st.SampleID])
				}
			}
		}
		one.Close()
	}
}

// The shape a set exists for: one question spanning several chromosomes,
// answered by the members that hold them.
func TestSetClassifyManySpansMembers(t *testing.T) {
	samples := []string{"S1", "S2", "S3"}
	dir := buildSet(t, []string{"chr1", "chr2", "chr7"}, samples)
	set, err := OpenSet(context.Background(), dir)
	if err != nil {
		t.Fatal(err)
	}
	defer set.Close()

	var loci []Locus
	for _, chrom := range []string{"chr1", "chr2", "chr7"} {
		loci = append(loci, Locus{Chrom: chrom, Pos: 100, Ref: "G", Alt: "A"})
	}
	// And one on a chromosome no member holds.
	loci = append(loci, Locus{Chrom: "chrX", Pos: 100, Ref: "G", Alt: "A"})

	got, err := set.ClassifyMany(loci, Gate{MinDP: 10})
	if err != nil {
		t.Fatal(err)
	}
	if len(got) != len(loci) {
		t.Fatalf("got %d loci back, asked about %d", len(got), len(loci))
	}
	for _, l := range loci[:3] {
		if len(got[l]) != len(samples) {
			t.Errorf("%s returned %d states, want %d", l, len(got[l]), len(samples))
		}
		if got[l][0].State != StateCarrier {
			t.Errorf("%s %s = %q, want carrier", l, got[l][0].SampleID, got[l][0].State)
		}
	}
	// A chromosome nobody holds is NOT ASSAYED for everyone, which is an answer
	// rather than a failure -- a set covering the autosomes asked about chrX
	// says "nothing here", not "error".
	for _, st := range got[loci[3]] {
		if st.State != StateNotAssayed {
			t.Errorf("chrX %s = %q, want not assayed", st.SampleID, st.State)
		}
	}
}

// A set whose members disagree would answer with a different population per
// chromosome, and nothing outside it could see that. So it cannot be built.
func TestSetRefusesDisagreeingMembers(t *testing.T) {
	t.Run("different rosters", func(t *testing.T) {
		dir := t.TempDir()
		chromStore(t, dir, "chr1", []string{"S1", "S2", "S3"}, 10)
		chromStore(t, dir, "chr2", []string{"S1", "S2"}, 10)
		_, err := BuildSet(context.Background(), dir, []string{"chr1", "chr2"}, "t", "t")
		if err == nil {
			t.Fatal("a set was built from members holding different samples")
		}
	})
	t.Run("different depth gates", func(t *testing.T) {
		dir := t.TempDir()
		samples := []string{"S1", "S2"}
		chromStore(t, dir, "chr1", samples, 10)
		chromStore(t, dir, "chr2", samples, 20)
		_, err := BuildSet(context.Background(), dir, []string{"chr1", "chr2"}, "t", "t")
		if err == nil {
			t.Fatal("a set was built from members converted at different depth gates; " +
				"they do not mean the same thing by callable")
		}
	})
	t.Run("overlapping chromosomes", func(t *testing.T) {
		dir := t.TempDir()
		samples := []string{"S1", "S2"}
		chromStore(t, dir, "chr1", samples, 10)
		// A second store holding the same chromosome, under another name.
		base := filepath.Join(dir, "dup")
		w, _ := NewWriter(base, WriterOpts{Samples: samples, MinDP: 10})
		_ = w.WriteSite(Site{Chrom: "chr1", Pos: 500, Ref: "T", Alt: "C"})
		_ = w.Finish()
		_, err := BuildSet(context.Background(), dir, []string{"chr1", "dup"}, "t", "t")
		if err == nil {
			t.Fatal("a set was built with two members holding chr1; a locus would have " +
				"two answers and no rule for choosing")
		}
	})
}

// A member swapped in from another conversion is caught when it opens, rather
// than requiring every manifest to be read at OpenSet -- which would undo the
// laziness the set exists for.
func TestSetCatchesASwappedMemberOnOpen(t *testing.T) {
	samples := []string{"S1", "S2", "S3"}
	dir := buildSet(t, []string{"chr1", "chr2"}, samples)

	// Replace chr2 with a store holding a different roster.
	other := t.TempDir()
	chromStore(t, other, "chr2", []string{"X1", "X2", "X3"}, 10)
	if err := copyStore(filepath.Join(other, "chr2"), filepath.Join(dir, "chr2")); err != nil {
		t.Skipf("could not stage the swap: %v", err)
	}

	set, err := OpenSet(context.Background(), dir)
	if err != nil {
		t.Fatal(err)
	}
	defer set.Close()

	// chr1 still answers: laziness means the bad member has not been touched.
	if _, err := set.Classify(Locus{Chrom: "chr1", Pos: 100, Ref: "G", Alt: "A"}, Gate{MinDP: 10}); err != nil {
		t.Errorf("chr1 failed because chr2 is wrong: %v", err)
	}
	// chr2 refuses.
	_, err = set.Classify(Locus{Chrom: "chr2", Pos: 100, Ref: "G", Alt: "A"}, Gate{MinDP: 10})
	if err == nil {
		t.Fatal("a member with a different roster answered; its calls would be attributed " +
			"to the wrong people")
	}
}

// A set is opened by OpenStore like anything else, so a caller never has to ask
// which it has.
func TestOpenStoreReturnsASet(t *testing.T) {
	dir := buildSet(t, []string{"chr1", "chr2"}, []string{"S1", "S2"})
	st, err := OpenStore(context.Background(), dir, "")
	if err != nil {
		t.Fatal(err)
	}
	defer st.Close()
	if _, ok := st.(*VarSet); !ok {
		t.Fatalf("OpenStore returned %T, want *VarSet", st)
	}
	samples, err := st.Samples()
	if err != nil || len(samples) != 2 {
		t.Fatalf("samples = %v, %v", samples, err)
	}
}

func copyStore(from, to string) error {
	return execCopy(from, to)
}

func execCopy(from, to string) error {
	// Deliberately simple: the fixture stores are a handful of small files.
	return copyDir(from, to)
}

func copyDir(from, to string) error {
	entries, err := osReadDir(from)
	if err != nil {
		return err
	}
	if err := osMkdirAll(to); err != nil {
		return err
	}
	for _, e := range entries {
		src := filepath.Join(from, e)
		dst := filepath.Join(to, e)
		b, err := osReadFile(src)
		if err != nil {
			// A directory (a shard folder); recurse.
			if err2 := copyDir(src, dst); err2 != nil {
				return fmt.Errorf("%v / %v", err, err2)
			}
			continue
		}
		if err := osWriteFile(dst, b); err != nil {
			return err
		}
	}
	return nil
}

func osReadDir(p string) ([]string, error) {
	es, err := os.ReadDir(p)
	if err != nil {
		return nil, err
	}
	var out []string
	for _, e := range es {
		out = append(out, e.Name())
	}
	return out, nil
}
func osMkdirAll(p string) error { return os.MkdirAll(p, 0o755) }
func osReadFile(p string) ([]byte, error) {
	st, err := os.Stat(p)
	if err != nil {
		return nil, err
	}
	if st.IsDir() {
		return nil, fmt.Errorf("dir")
	}
	return os.ReadFile(p)
}
func osWriteFile(p string, b []byte) error { return os.WriteFile(p, b, 0o644) }
