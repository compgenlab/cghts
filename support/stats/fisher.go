// Package stats provides small statistical helpers used by genomics analyses:
// a 2x2 Fisher exact test and Phred/log2 conversions.
package stats

import "math"

// FisherExact computes p-values for 2x2 contingency tables, caching log-factorial
// values across calls. It ports ngsutilsj support.stats.FisherExact. A
// FisherExact is not safe for concurrent use.
type FisherExact struct {
	logVals []float64
}

// NewFisherExact returns a ready-to-use FisherExact.
func NewFisherExact() *FisherExact { return &FisherExact{} }

func (f *FisherExact) recalc(n int) {
	if f.logVals != nil && len(f.logVals) > n {
		return
	}
	start := 1
	if f.logVals == nil {
		f.logVals = make([]float64, n+1)
		f.logVals[0] = 0.0
	} else {
		start = len(f.logVals)
		tmp := make([]float64, n+1)
		copy(tmp, f.logVals)
		f.logVals = tmp
	}
	for i := start; i <= n; i++ {
		f.logVals[i] = f.logVals[i-1] + math.Log(float64(i))
	}
}

// Pvalue returns the exact hypergeometric p-value for the 2x2 table
// [[a,b],[c,d]].
func (f *FisherExact) Pvalue(a, b, c, d int) float64 {
	n := a + b + c + d
	f.recalc(n)
	return math.Exp((f.logVals[a+b] + f.logVals[c+d] + f.logVals[a+c] + f.logVals[b+d]) -
		(f.logVals[a] + f.logVals[b] + f.logVals[c] + f.logVals[d] + f.logVals[n]))
}

// TwoTailedPvalue returns the two-tailed p-value for the 2x2 table, summing the
// exact table and all less-likely tables on either tail.
func (f *FisherExact) TwoTailedPvalue(a, b, c, d int) float64 {
	exact := f.Pvalue(a, b, c, d)
	p := exact
	p += f.innerLeftTail(a, b, c, d, exact)
	p += f.innerRightTail(a, b, c, d, exact)
	return p
}

func (f *FisherExact) innerLeftTail(a, b, c, d int, thres float64) float64 {
	acc := 0.0
	for a > 0 && d > 0 {
		a, b, c, d = a-1, b+1, c+1, d-1
		p := f.Pvalue(a, b, c, d)
		if p < thres {
			acc += p
		}
	}
	return acc
}

func (f *FisherExact) innerRightTail(a, b, c, d int, thres float64) float64 {
	acc := 0.0
	for b > 0 && c > 0 {
		a, b, c, d = a+1, b-1, c-1, d+1
		p := f.Pvalue(a, b, c, d)
		if p < thres {
			acc += p
		}
	}
	return acc
}

// Phred converts a probability to a Phred-scaled score (-10*log10), clamped to
// [0, 255]. A value <= 0 returns 255; a value >= 1 returns 0.
func Phred(val float64) float64 {
	if val <= 0.0 {
		return 255.0
	}
	if val >= 1.0 {
		return 0.0
	}
	return -10 * math.Log10(val)
}

// Log2 returns the base-2 logarithm of val.
func Log2(val float64) float64 { return math.Log(val) / math.Ln2 }
