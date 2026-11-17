package varstore

import "strings"

// Detaching the strings a buffered row carries from wherever they came from.
//
// A VCF reader hands out CHROM, REF, ALT and FORMAT values as SUBSTRINGS OF THE
// RECORD'S LINE, and a Go substring keeps its whole backing string alive. A line
// carries every sample's genotype fields, so in a wide cohort it is megabytes
// long. The Writer buffers batchSize rows of each table before parquet-go copies
// their bytes, so until each flush every buffered row pinned an entire line:
// 8,192 sites of a 50,000-sample cohort held 8-24 GB, which is how a conversion
// whose real working set is a few hundred megabytes was OOM-killed at 8 GB. Memory
// grew with cohort width times batch size, never with the input's length, which
// is why small fixtures could not show it.
//
// So every string a row keeps is copied as it is buffered. The repeated ones --
// chromosome and sample -- are interned, so a run of rows on one chromosome or
// a roster of samples costs one copy each rather than one per row; REF, ALT and
// GT are short and copied outright.
type owned struct {
	lastChrom string
	samples   map[string]string
}

// chrom returns a detached copy of a chromosome name, reusing the last one when
// it repeats -- which is nearly always, since input is coordinate sorted.
func (o *owned) chrom(s string) string {
	if s != o.lastChrom {
		o.lastChrom = strings.Clone(s)
	}
	return o.lastChrom
}

// sample returns a detached, interned copy of a sample ID.
func (o *owned) sample(s string) string {
	if v, ok := o.samples[s]; ok {
		return v
	}
	if o.samples == nil {
		o.samples = map[string]string{}
	}
	v := strings.Clone(s)
	o.samples[v] = v
	return v
}

func (o *owned) site(s *Site) {
	s.Chrom = o.chrom(s.Chrom)
	s.Ref = strings.Clone(s.Ref)
	s.Alt = strings.Clone(s.Alt)
}

func (o *owned) call(c *Call) {
	c.SampleID = o.sample(c.SampleID)
	c.Chrom = o.chrom(c.Chrom)
	c.Ref = strings.Clone(c.Ref)
	c.Alt = strings.Clone(c.Alt)
	c.GT = strings.Clone(c.GT)
}
