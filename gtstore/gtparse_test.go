package gtstore

import (
	"strconv"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/vcf"
)

// The allocation-free genotype readers must agree with the Split-based ones
// they replaced, on every shape a GT takes in the wild and some it should not.
func oldAlleles(gt string) []string { return strings.Split(strings.ReplaceAll(gt, "|", "/"), "/") }

func oldSplitGT(gt string, altIdx int) (string, bool) {
	if gt == "" {
		return ".", false
	}
	sep := "/"
	if strings.Contains(gt, "|") {
		sep = "|"
	}
	carrier := false
	var out []string
	for _, a := range oldAlleles(gt) {
		switch {
		case a == "." || a == "":
			out = append(out, ".")
		case a == "0":
			out = append(out, "0")
		default:
			n, err := strconv.Atoi(a)
			if err == nil && n == altIdx {
				out = append(out, "1")
				carrier = true
			} else {
				out = append(out, ".")
			}
		}
	}
	return strings.Join(out, sep), carrier
}

func oldAddAlleleCounts(gt string, acc []int32) (an int32) {
	for _, a := range oldAlleles(gt) {
		if a == "" || a == "." {
			continue
		}
		an++
		n, err := strconv.Atoi(a)
		if err != nil || n <= 0 {
			continue
		}
		if n-1 < len(acc) {
			acc[n-1]++
		}
	}
	return an
}

var gtCases = []string{
	"", ".", "0", "1", "2", "0/0", "0/1", "1/0", "1/1", "1|0", "0|1", "1|1", "./.", ".|.",
	"./1", "1/.", "1/2", "2/1", "2|.", "0/2", "10/3", "3/10", "0/1/2", "1|2|3", "0/", "/1", "/",
	"//", "a/1", "1/b", "-1/1", "+1/0", "01/1", "00/0", "999999999999/1", " 1/0", "0 /1",
	"1/1/1/1", "0|0|0", "1/|2",
}

func TestGenotypeReadersMatchSplitVersions(t *testing.T) {
	for _, gt := range gtCases {
		old := oldAlleles(gt)
		if got := gtAlleles(gt); strings.Join(got, "\x00") != strings.Join(old, "\x00") || len(got) != len(old) {
			t.Errorf("gtAlleles(%q) = %q, want %q", gt, got, old)
		}
		var oHas, oAlt bool
		oRef, oSeen := true, false
		for _, a := range old {
			oHas = oHas || (a != "" && a != ".")
			oAlt = oAlt || (a != "" && a != "." && a != "0")
			if a != "" {
				if a != "0" {
					oRef = false
				}
				oSeen = true
			}
		}
		if HasCall(gt) != oHas {
			t.Errorf("HasCall(%q) = %v", gt, !oHas)
		}
		if IsAltCarrier(gt) != oAlt {
			t.Errorf("IsAltCarrier(%q) = %v", gt, !oAlt)
		}
		if IsHomRef(gt) != (oRef && oSeen) {
			t.Errorf("IsHomRef(%q) = %v", gt, !(oRef && oSeen))
		}
		for alt := 1; alt <= 3; alt++ {
			wr, wc := oldSplitGT(gt, alt)
			gr, gc := SplitGT(gt, alt)
			if gr != wr || gc != wc {
				t.Errorf("SplitGT(%q, %d) = %q,%v; want %q,%v", gt, alt, gr, gc, wr, wc)
			}
			if carriesAllele(gt, alt) != wc {
				t.Errorf("carriesAllele(%q, %d) = %v, want %v", gt, alt, !wc, wc)
			}
		}
		wa, ga := make([]int32, 3), make([]int32, 3)
		wn, gn := oldAddAlleleCounts(gt, wa), AddAlleleCounts(gt, ga)
		if wn != gn || strconvAll(wa) != strconvAll(ga) {
			t.Errorf("AddAlleleCounts(%q) = %d %v, want %d %v", gt, gn, ga, wn, wa)
		}
	}
}

func strconvAll(v []int32) string {
	var b strings.Builder
	for _, x := range v {
		b.WriteString(strconv.Itoa(int(x)))
		b.WriteByte(' ')
	}
	return b.String()
}

// A non-carrier is the common case in a wide cohort, and must not allocate.
func TestNonCarrierChecksDoNotAllocate(t *testing.T) {
	acc := make([]int32, 2)
	n := testing.AllocsPerRun(100, func() {
		_ = HasCall("0/0")
		_ = carriesAllele("0|0", 1)
		_ = AddAlleleCounts("0/1", acc)
		_ = IsHomRef("0/0")
	})
	if n != 0 {
		t.Errorf("genotype checks allocated %.0f times per run, want 0", n)
	}
}

// SampleReader must return what the map-based ReadSample did: Sample(i), then
// Get on each of the four keys.
func TestSampleReaderMatchesSampleMap(t *testing.T) {
	const vcfText = "##fileformat=VCFv4.2\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD\n" +
		"chr1\t100\t.\tA\tG\t.\t.\t.\tGT:AD:DP:GQ:PL\t0/1:10,12:22:99:1,2,3\t0/0:30,0:30\t./.\t1|1::.:x\n" +
		"chr1\t200\t.\tA\tG\t.\t.\tIMPUTED\tGT:DS:GP\t0|0:0.1:1,0,0\t0|1\t.\t\n" +
		"chr1\t300\t.\tA\tG\t.\t.\t.\tDP:GT\t7:0/1\t:1/1\t3\t.:.\n"
	r, err := vcf.NewVcfReader(strings.NewReader(vcfText))
	if err != nil {
		t.Fatal(err)
	}
	nrec := 0
	for {
		rec, err := r.NextRecord()
		if err != nil || rec == nil {
			break
		}
		nrec++
		sr := NewSampleReader(rec)
		for i := 0; i < rec.NumSamples(); i++ {
			want := SampleFields{DP: Missing, GQ: Missing}
			attrs, _ := rec.Sample(i)
			if v, ok := attrs.Get("GT"); ok {
				want.GT = v.String()
			}
			if v, ok := attrs.Get("DP"); ok {
				want.DP = atoi32(v.String())
			}
			if v, ok := attrs.Get("GQ"); ok {
				want.GQ = atoi32(v.String())
			}
			if v, ok := attrs.Get("AD"); ok {
				want.AD = v.String()
			}
			if got := sr.Read(i); got != want {
				t.Errorf("%s:%d sample %d: Read %+v, map %+v", rec.Chrom, rec.Pos, i, got, want)
			}
		}
	}
	if nrec != 3 {
		t.Fatalf("read %d records, want 3", nrec)
	}
}
