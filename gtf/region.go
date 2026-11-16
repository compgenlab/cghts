package gtf

import "github.com/compgenlab/cghts/bed"

// GenicRegion classifies a position (or region/read) relative to gene
// annotations: which part of a gene it falls in — coding exon, UTR, intron,
// non-coding exon, splice junction — plus intergenic and mitochondrial, each
// with a sense and an anti-sense (anti_*) variant. It ports ngsutilsj
// GenicRegion. The values and their priority order (lower ord = higher
// priority, used as the tie-breaker in FindGenicRegionForRegion) match the Java
// enum's declaration order exactly.
type GenicRegion struct {
	Name     string
	Code     string
	IsGene   bool
	IsExon   bool
	IsCoding bool
	IsSense  bool
	ord      int
}

// Description returns the human-readable label, appending " (anti-sense)" for
// anti-sense gene regions. Ports GenicRegion.getDescription().
func (r GenicRegion) Description() string {
	if r.IsGene && !r.IsSense {
		return r.Name + " (anti-sense)"
	}
	return r.Name
}

// The genic regions, declared in ngsutilsj priority order. Sense block first
// (Junction is highest priority), then Intergenic/Mitochondrial, then the
// anti-sense block. Codes are byte-for-byte identical to GenicRegion.java.
var (
	Junction      = GenicRegion{Name: "Junction", Code: "junction", IsGene: true, IsExon: true, IsCoding: true, IsSense: true}
	Coding        = GenicRegion{Name: "Coding", Code: "coding_exon", IsGene: true, IsExon: true, IsCoding: true, IsSense: true}
	UTR5          = GenicRegion{Name: "5'UTR", Code: "5_utr", IsGene: true, IsExon: true, IsCoding: false, IsSense: true}
	UTR3          = GenicRegion{Name: "3'UTR", Code: "3_utr", IsGene: true, IsExon: true, IsCoding: false, IsSense: true}
	NCJunction    = GenicRegion{Name: "Non-coding junction", Code: "nc_junction", IsGene: true, IsExon: true, IsCoding: false, IsSense: true}
	NCExon        = GenicRegion{Name: "Non-coding exon", Code: "nc_exon", IsGene: true, IsExon: true, IsCoding: false, IsSense: true}
	CodingIntron  = GenicRegion{Name: "Coding intron", Code: "coding_intron", IsGene: true, IsExon: false, IsCoding: true, IsSense: true}
	UTR5Intron    = GenicRegion{Name: "5'UTR intron", Code: "5_utr_intron", IsGene: true, IsExon: false, IsCoding: false, IsSense: true}
	UTR3Intron    = GenicRegion{Name: "3'UTR intron", Code: "3_utr_intron", IsGene: true, IsExon: false, IsCoding: false, IsSense: true}
	NCIntron      = GenicRegion{Name: "Non-coding intron", Code: "nc_intron", IsGene: true, IsExon: false, IsCoding: false, IsSense: true}
	Intergenic    = GenicRegion{Name: "Intergenic", Code: "intergenic", IsGene: false, IsExon: false, IsCoding: false, IsSense: false}
	Mitochondrial = GenicRegion{Name: "Mitochondrial", Code: "mitochondrial", IsGene: false, IsExon: false, IsCoding: false, IsSense: false}

	JunctionAnti     = GenicRegion{Name: "Junction", Code: "anti_junction", IsGene: true, IsExon: true, IsCoding: true, IsSense: false}
	CodingAnti       = GenicRegion{Name: "Coding", Code: "anti_coding_exon", IsGene: true, IsExon: true, IsCoding: true, IsSense: false}
	UTR5Anti         = GenicRegion{Name: "5'UTR", Code: "anti_5_utr", IsGene: true, IsExon: true, IsCoding: false, IsSense: false}
	UTR3Anti         = GenicRegion{Name: "3'UTR", Code: "anti_3_utr", IsGene: true, IsExon: true, IsCoding: false, IsSense: false}
	NCJunctionAnti   = GenicRegion{Name: "Non-coding junction", Code: "anti_nc_junction", IsGene: true, IsExon: true, IsCoding: false, IsSense: false}
	NCExonAnti       = GenicRegion{Name: "Non-coding exon", Code: "anti_nc_exon", IsGene: true, IsExon: true, IsCoding: false, IsSense: false}
	CodingIntronAnti = GenicRegion{Name: "Coding intron", Code: "anti_coding_intron", IsGene: true, IsExon: false, IsCoding: true, IsSense: false}
	UTR5IntronAnti   = GenicRegion{Name: "5'UTR intron", Code: "anti_5_utr_intron", IsGene: true, IsExon: false, IsCoding: false, IsSense: false}
	UTR3IntronAnti   = GenicRegion{Name: "3'UTR intron", Code: "anti_3_utr_intron", IsGene: true, IsExon: false, IsCoding: false, IsSense: false}
	NCIntronAnti     = GenicRegion{Name: "Non-coding intron", Code: "anti_nc_intron", IsGene: true, IsExon: false, IsCoding: false, IsSense: false}
)

// genicRegionOrder lists every region in priority order. init() stamps each
// region's ord from its index here so equal regions compare equal (ord is part
// of the struct) and FindGenicRegionForRegion can break ties by priority.
var genicRegionOrder = []*GenicRegion{
	&Junction, &Coding, &UTR5, &UTR3, &NCJunction, &NCExon,
	&CodingIntron, &UTR5Intron, &UTR3Intron, &NCIntron,
	&Intergenic, &Mitochondrial,
	&JunctionAnti, &CodingAnti, &UTR5Anti, &UTR3Anti, &NCJunctionAnti, &NCExonAnti,
	&CodingIntronAnti, &UTR5IntronAnti, &UTR3IntronAnti, &NCIntronAnti,
}

func init() {
	for i, r := range genicRegionOrder {
		r.ord = i
	}
}

// GenicRegions returns every region in priority order. Useful for initializing
// a per-region counter (e.g. the future BAM read counter / sam-stats).
func GenicRegions() []GenicRegion {
	out := make([]GenicRegion, len(genicRegionOrder))
	for i, r := range genicRegionOrder {
		out[i] = *r
	}
	return out
}

// FindGenicRegionForPos classifies a single 0-based position on ref relative to
// the gene model. strand is the query's strand: pass bed.StrandNone for an
// unstranded query (e.g. a variant), in which case only sense codes are
// returned; with a real strand, overlaps on the opposite strand return anti_*
// codes. When geneID is non-empty, classification is restricted to that gene
// (used to report a per-gene region for a position overlapping several genes).
// See classifyAt for the rules.
func (s *AnnotationSource) FindGenicRegionForPos(ref string, pos int, strand bed.Strand, geneID string) GenicRegion {
	genes := s.FindGenes(ref, pos, pos+1)
	if geneID != "" {
		var only []*Gene
		for _, g := range genes {
			if g.GeneID == geneID {
				only = append(only, g)
			}
		}
		genes = only
	}
	return classifyAt(genes, ref, pos, strand)
}

// classifyAt classifies a 0-based position on ref given the genes overlapping
// it. It follows ngsutilsj GTFAnnotationSource.findGenicRegionForPos -- the same
// categories, the same precedence (coding, 5'UTR, 3'UTR, non-coding exon, then
// coding/5'UTR/3'UTR/non-coding intron, then intergenic) and the same rule that
// a category is anti-sense when any gene contributing to it is -- except that
// the answer does not depend on the order genes and transcripts are visited.
//
// The rules:
//   - exonic if any transcript of any overlapping gene has an exon or CDS here;
//   - coding if any transcript has CDS here;
//   - otherwise, if exonic, a coding transcript makes it 5'UTR or 3'UTR by
//     where it lies relative to that transcript's CDS, 5'UTR taking precedence;
//     with neither it is a non-coding exon;
//   - if not exonic, a coding transcript makes it a coding, 5'UTR or 3'UTR
//     intron by the same test; with none of those it is a non-coding intron;
//   - intergenic where no gene overlaps, mitochondrial on chrM/M.
//
// Java keeps its isExon/isCoding flags across the whole loop, so a transcript is
// judged partly by the ones visited before it: whether a coding transcript reads
// a position as UTR or intron depends on whether an exon hit was seen earlier.
// It visits transcripts in java.util.HashMap order and overlapping genes in
// HashSet order (JVM identity hash codes), and cghts used to visit them in
// sorted-ID order, so wherever transcripts disagreed the answer was an accident
// of naming -- a renamed transcript could turn a 5'UTR into a non-coding exon.
// Deciding "exonic" over all transcripts first keeps Java's evident intent, an
// exonic position being UTR by the coding transcripts' CDS, without the order.
// See compgenlab/cghts#61.
//
// Kept as in Java, quirks included: the coding-sequence bounds are the CDS rows'
// (the stop codon lies outside them, so its first base is neither CDS nor UTR),
// a coding transcript classifies a position its exons do not cover by its CDS
// bounds even outside its own span, and "exonic" spans genes, so another gene's
// exon makes a coding gene's position a UTR rather than an intron.
func classifyAt(genes []*Gene, ref string, pos int, strand bed.Strand) GenicRegion {
	if ref == "chrM" || ref == "M" {
		return Mitochondrial
	}

	var (
		isGene, isExon, isCoding, isUtr3, isUtr5                bool
		isCodingIntron, isUtr3Intron, isUtr5Intron              bool
		isGeneRev, isExonRev, isCodingRev, isUtr3Rev, isUtr5Rev bool
	)

	// A 1-bp unstranded probe; the antisense decision is made against each
	// gene's strand so that genes on either strand are found.
	probe := span{ref: ref, start: pos, end: pos + 1, strand: bed.StrandNone}
	covers := func(xs []*Exon, ref string) bool {
		for _, x := range xs {
			if x.region(ref).contains(probe) {
				return true
			}
		}
		return false
	}

	// Pass 1: is the position exonic in any transcript of any overlapping gene?
	// Java answers this with a flag that grows as it walks, so a transcript is
	// judged by whichever ones it happened to visit first; asking it of all of
	// them up front gives the same answer whatever the order.
	exonic := false
	for _, gene := range genes {
		for _, txpt := range gene.SortedTranscripts() {
			if covers(txpt.CDS, gene.Ref) || covers(txpt.Exons, gene.Ref) {
				exonic = true
			}
		}
	}

	// Pass 2: each transcript's own findings, with "exonic" deciding -- as in
	// Java -- whether a coding transcript reads the position as UTR or intron.
	for _, gene := range genes {
		isGene = true
		anti := strand != bed.StrandNone && gene.Strand != strand
		if anti {
			isGeneRev = true
		}
		// Sorted only for determinism of iteration; with nothing carried
		// between transcripts the order cannot change the answer.
		for _, txpt := range gene.SortedTranscripts() {
			if !txpt.HasCDS() {
				if covers(txpt.Exons, gene.Ref) {
					isExon = true
					isExonRev = isExonRev || anti
				}
				continue
			}

			if covers(txpt.CDS, gene.Ref) {
				isExon, isCoding = true, true
				isCodingRev = isCodingRev || anti
				continue // inside [CdsStart, CdsEnd): neither UTR nor intron
			}
			if covers(txpt.Exons, gene.Ref) {
				isExon = true
				isExonRev = isExonRev || anti
			}

			// 5' and 3' swap with the gene's strand.
			before, after := &isUtr5, &isUtr3
			beforeRev, afterRev := &isUtr5Rev, &isUtr3Rev
			beforeIntron, afterIntron := &isUtr5Intron, &isUtr3Intron
			if gene.Strand != bed.StrandPlus {
				before, after = after, before
				beforeRev, afterRev = afterRev, beforeRev
				beforeIntron, afterIntron = afterIntron, beforeIntron
			}

			if exonic {
				switch {
				case pos < txpt.CdsStart:
					*before = true
					*beforeRev = *beforeRev || anti
				case pos > txpt.CdsEnd:
					*after = true
					*afterRev = *afterRev || anti
				}
				continue
			}
			switch {
			case pos < txpt.CdsStart:
				*beforeIntron = true
			case pos > txpt.CdsEnd:
				*afterIntron = true
			default:
				isCodingIntron = true
			}
		}
	}

	pick := func(rev bool, anti, sense GenicRegion) GenicRegion {
		if rev {
			return anti
		}
		return sense
	}

	switch {
	case isCoding:
		return pick(isCodingRev, CodingAnti, Coding)
	case isUtr5:
		return pick(isUtr5Rev, UTR5Anti, UTR5)
	case isUtr3:
		return pick(isUtr3Rev, UTR3Anti, UTR3)
	case isExon:
		return pick(isExonRev, NCExonAnti, NCExon)
	case isGene:
		if isGeneRev {
			switch {
			case isCodingIntron:
				return CodingIntronAnti
			case isUtr5Intron:
				return UTR5IntronAnti
			case isUtr3Intron:
				return UTR3IntronAnti
			default:
				return NCIntronAnti
			}
		}
		switch {
		case isCodingIntron:
			return CodingIntron
		case isUtr5Intron:
			return UTR5Intron
		case isUtr3Intron:
			return UTR3Intron
		default:
			return NCIntron
		}
	default:
		return Intergenic
	}
}

// FindGenicRegionForRegion classifies a [start, end) region by classifying both
// ends and reconciling them. A region with one exonic and one non-exonic end is
// crossing a splice junction. Ports findGenicRegionForRegion.
func (s *AnnotationSource) FindGenicRegionForRegion(ref string, start, end int, strand bed.Strand) GenicRegion {
	genStart := s.FindGenicRegionForPos(ref, start, strand, "")
	genEnd := s.FindGenicRegionForPos(ref, end-1, strand, "")

	if genStart == genEnd {
		return genStart
	}
	if genStart.IsGene && !genEnd.IsGene {
		return genStart
	}
	if !genStart.IsGene && genEnd.IsGene {
		return genEnd
	}
	if genStart.IsExon && !genEnd.IsExon {
		return junctionFor(genStart)
	}
	if !genStart.IsExon && genEnd.IsExon {
		return junctionFor(genEnd)
	}
	if genStart.IsCoding && !genEnd.IsCoding {
		return genStart
	}
	if !genStart.IsCoding && genEnd.IsCoding {
		return genEnd
	}
	if genStart.ord < genEnd.ord {
		return genStart
	}
	return genEnd
}

// Junctionize upgrades a gene region to its junction code, for a spliced read
// or region that crosses an intron (the caller detects the splice, e.g. a CIGAR
// N operator). Non-gene regions pass through unchanged. Ports the junction
// promotion in findGenicRegion.
func Junctionize(reg GenicRegion) GenicRegion {
	if !reg.IsGene {
		return reg
	}
	return junctionFor(reg)
}

func junctionFor(reg GenicRegion) GenicRegion {
	switch {
	case reg.IsCoding && reg.IsSense:
		return Junction
	case reg.IsCoding:
		return JunctionAnti
	case reg.IsSense:
		return NCJunction
	default:
		return NCJunctionAnti
	}
}

func pick(rev bool, anti, sense GenicRegion) GenicRegion {
	if rev {
		return anti
	}
	return sense
}
