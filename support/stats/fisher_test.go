package stats

import (
	"math"
	"testing"
)

func approx(a, b, eps float64) bool { return math.Abs(a-b) < eps }

func TestFisherExactPvalue(t *testing.T) {
	f := NewFisherExact()
	// Pvalue(1,1,1,1) = C(2,1)C(2,1)/C(4,2) = 4/6.
	if got := f.Pvalue(1, 1, 1, 1); !approx(got, 4.0/6.0, 1e-12) {
		t.Errorf("Pvalue(1,1,1,1) = %v, want %v", got, 4.0/6.0)
	}
	// Pvalue(3,1,1,3) = C(4,3)C(4,1)/C(8,4) = 16/70.
	if got := f.Pvalue(3, 1, 1, 3); !approx(got, 16.0/70.0, 1e-12) {
		t.Errorf("Pvalue(3,1,1,3) = %v, want %v", got, 16.0/70.0)
	}
}

func TestFisherTwoTailed(t *testing.T) {
	f := NewFisherExact()
	// For [[3,1],[1,3]] the algorithm (strict < threshold, one walk per tail)
	// sums the exact table (16/70) plus the a=0 and a=4 tables (1/70 each):
	// 18/70 = 0.2571428571...
	if got := f.TwoTailedPvalue(3, 1, 1, 3); !approx(got, 18.0/70.0, 1e-12) {
		t.Errorf("TwoTailedPvalue(3,1,1,3) = %v, want %v", got, 18.0/70.0)
	}
	// A symmetric balanced table is maximally non-significant (p == 1).
	if got := f.TwoTailedPvalue(10, 10, 10, 10); !approx(got, 1.0, 1e-9) {
		t.Errorf("TwoTailedPvalue(10,10,10,10) = %v, want 1.0", got)
	}
}

func TestPhred(t *testing.T) {
	cases := []struct {
		in, want float64
	}{
		{0, 255}, {-1, 255}, {1, 0}, {2, 0}, {0.1, 10}, {0.01, 20}, {0.001, 30},
	}
	for _, c := range cases {
		if got := Phred(c.in); !approx(got, c.want, 1e-9) {
			t.Errorf("Phred(%v) = %v, want %v", c.in, got, c.want)
		}
	}
}

func TestLog2(t *testing.T) {
	if got := Log2(8); !approx(got, 3, 1e-12) {
		t.Errorf("Log2(8) = %v, want 3", got)
	}
	if got := Log2(1024); !approx(got, 10, 1e-12) {
		t.Errorf("Log2(1024) = %v, want 10", got)
	}
}
