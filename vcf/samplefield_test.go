package vcf

import (
	"math/rand"
	"strings"
	"testing"
)

// SampleField must report exactly what Sample(i).Get(key).String() does, since
// readers switch between the two for speed alone.
func TestSampleFieldMatchesSample(t *testing.T) {
	lines := []string{
		"chr1\t100\t.\tA\tG\t.\t.\t.\tGT:AD:DP:GQ:PL\t0/1:10,12:22:99:200,0,180\t0/0:30,0:30\t./.\t1|1:0,15:15:45:500,45,0\t0/1::18::",
		"chr1\t200\t.\tA\tG\t.\t.\tIMPUTED\tGT:DS:GP\t0|0:0.001:0.99,0.01,0\t0|1:1.0:0,1,0\t.\t\t1|1",
		"chr1\t300\t.\tA\tG,T\t.\t.\t.\tGT:DP:DP\t1/2:5:9\t0/0:7",
		"chr1\t400\t.\tA\tG\t.\t.\t.\tGT\t0/1\t.",
	}
	for _, line := range lines {
		for _, edit := range []bool{false, true} {
			want, err := newRecord(line, nil)
			if err != nil {
				t.Fatal(err)
			}
			got, _ := newRecord(line, nil)
			keys := append(want.FormatKeys(), "NOPE")
			for i := 0; i < want.NumSamples(); i++ {
				if edit {
					// A sample parsed and edited through Sample must read back
					// the edit, not the raw column.
					a, _ := got.Sample(i)
					a.Set("GT", "9/9")
					w, _ := want.Sample(i)
					w.Set("GT", "9/9")
				}
				attrs, _ := want.Sample(i)
				for _, key := range keys {
					exp := missing
					if v, ok := attrs.Get(key); ok {
						exp = v.String()
					}
					if g := got.SampleField(i, got.FormatIndex(key)); g != exp {
						t.Errorf("%q sample %d %s (edited=%v): SampleField %q, Sample %q",
							line[:12], i, key, edit, g, exp)
					}
				}
			}
		}
	}
}

func TestSampleFieldRandom(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	vals := []string{"", ".", "0/1", "12", "3,4", "a:b"}
	for n := 0; n < 500; n++ {
		nkeys := 1 + rng.Intn(5)
		keys := make([]string, nkeys)
		for k := range keys {
			keys[k] = string(rune('A' + k))
		}
		cols := []string{"chr1", "1", ".", "A", "G", ".", ".", ".", strings.Join(keys, ":")}
		for s := 0; s < 4; s++ {
			nv := rng.Intn(nkeys + 2)
			f := make([]string, nv)
			for k := range f {
				f[k] = vals[rng.Intn(len(vals))]
			}
			cols = append(cols, strings.Join(f, ":"))
		}
		line := strings.Join(cols, "\t")
		want, _ := newRecord(line, nil)
		got, _ := newRecord(line, nil)
		for i := 0; i < want.NumSamples(); i++ {
			attrs, _ := want.Sample(i)
			for _, key := range keys {
				v, _ := attrs.Get(key)
				if g := got.SampleField(i, got.FormatIndex(key)); g != v.String() {
					t.Fatalf("%q sample %d key %s: SampleField %q, Sample %q", line, i, key, g, v.String())
				}
			}
		}
	}
}
