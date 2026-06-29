package vcf

// ContigConverter resolves contig names from any supported naming scheme
// (UCSC "chr1", Ensembl "1", NCBI RefSeq "NC_000001.11") to the spelling used by
// a target namespace, matching by canonical chromosome identity (see
// [CanonicalContig]). It is built from the target's contig names — an annotation
// file's RefNames(), a BAM header's @SQ names, a .fai, or any literal list — so
// resolution is bidirectional and scheme-agnostic: whatever the target calls a
// contig, an incoming name that shares its canonical identity resolves to it.
type ContigConverter struct {
	exact       map[string]bool   // target names, for the exact-match fast path
	byCanonical map[string]string // canonical key -> target name
}

// NewContigConverter indexes the target contig names by canonical identity. When
// two target names share a canonical key (unusual), the first one wins.
func NewContigConverter(targetNames []string) *ContigConverter {
	c := &ContigConverter{
		exact:       make(map[string]bool, len(targetNames)),
		byCanonical: make(map[string]string, len(targetNames)),
	}
	for _, name := range targetNames {
		c.exact[name] = true
		if key, ok := CanonicalContig(name); ok {
			if _, dup := c.byCanonical[key]; !dup {
				c.byCanonical[key] = name
			}
		}
	}
	return c
}

// Resolve returns the target name equivalent to chrom. An exact match wins
// (zero-cost, no behavior change when names already agree); otherwise the name
// is matched by canonical identity. ok is false when no target contig shares
// chrom's canonical identity, in which case the caller should skip the lookup.
func (c *ContigConverter) Resolve(chrom string) (name string, ok bool) {
	if c.exact[chrom] {
		return chrom, true
	}
	if key, kok := CanonicalContig(chrom); kok {
		if target, tok := c.byCanonical[key]; tok {
			return target, true
		}
	}
	return "", false
}
