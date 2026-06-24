package annotate

import (
	"fmt"
	"strconv"
	"strings"

	"github.com/compgenlab/hts/support/stats"
	"github.com/compgenlab/hts/vcf"
)

// requireSAC verifies the header declares an SAC FORMAT field (Number=".",
// Type=Integer), matching the GATK strand-allele-counts field the strand
// annotators consume.
func requireSAC(h *vcf.VcfHeader) error {
	d, ok := h.FormatDef("SAC")
	if !ok || d.Number != "." || d.Type != "Integer" {
		return fmt.Errorf("annotate: \"SAC\" FORMAT annotation missing")
	}
	return nil
}

// parseSAC parses a comma-separated SAC value (ref+,ref-,alt1+,alt1-,...) into
// ints.
func parseSAC(s string) ([]int, error) {
	parts := strings.Split(s, ",")
	out := make([]int, len(parts))
	for i, p := range parts {
		n, err := strconv.Atoi(p)
		if err != nil {
			return nil, fmt.Errorf("annotate: invalid SAC value %q: %w", p, err)
		}
		out[i] = n
	}
	return out, nil
}

// round formats a float with the given number of decimal places, matching
// ngsutilsj's String.format("%.Nf", val).
func round(val float64, places int) string {
	return strconv.FormatFloat(val, 'f', places, 64)
}

// VariantAlleleFrequency adds FORMAT CG_VAF (per-alt allele frequency from SAC).
type VariantAlleleFrequency struct{ closeNoop }

// NewVAF returns a variant-allele-frequency annotator (--vaf).
func NewVAF() *VariantAlleleFrequency { return &VariantAlleleFrequency{} }

// SetupHeader requires SAC and declares CG_VAF.
func (a *VariantAlleleFrequency) SetupHeader(h *vcf.VcfHeader) error {
	if err := requireSAC(h); err != nil {
		return err
	}
	h.AddFormat(formatDef("CG_VAF", "A", "Float", "Allele frequency for all alt-alleles"))
	return nil
}

// Annotate computes per-alt VAF for each sample with SAC.
func (a *VariantAlleleFrequency) Annotate(rec *vcf.VcfRecord) error {
	for i := 0; i < rec.NumSamples(); i++ {
		s, err := rec.Sample(i)
		if err != nil {
			return err
		}
		v, ok := s.Get("SAC")
		if !ok {
			continue
		}
		sac, err := parseSAC(v.String())
		if err != nil {
			return err
		}
		total := 0
		for j := 0; j+1 < len(sac); j += 2 {
			total += sac[j] + sac[j+1]
		}
		var outs []string
		for j := 2; j+1 < len(sac); j += 2 {
			pm := sac[j] + sac[j+1]
			if pm == 0 {
				outs = append(outs, "0.0")
			} else {
				outs = append(outs, round(float64(pm)/float64(total), 3))
			}
		}
		if len(outs) == 0 {
			if err := rec.AddFormat(i, "CG_VAF", "."); err != nil {
				return err
			}
		} else if err := rec.AddFormat(i, "CG_VAF", strings.Join(outs, ",")); err != nil {
			return err
		}
	}
	return nil
}

// MinorStrandPct adds FORMAT CG_SBPCT (percent of alt reads on the minor strand).
type MinorStrandPct struct{ closeNoop }

// NewMinorStrand returns a minor-strand-percentage annotator (--minor-strand).
func NewMinorStrand() *MinorStrandPct { return &MinorStrandPct{} }

// SetupHeader requires SAC and declares CG_SBPCT.
func (a *MinorStrandPct) SetupHeader(h *vcf.VcfHeader) error {
	if err := requireSAC(h); err != nil {
		return err
	}
	h.AddFormat(formatDef("CG_SBPCT", "A", "Float", "Percent of alt-allele reads on the minor strand"))
	return nil
}

// Annotate computes per-alt minor-strand percentage for each sample with SAC.
func (a *MinorStrandPct) Annotate(rec *vcf.VcfRecord) error {
	for i := 0; i < rec.NumSamples(); i++ {
		s, err := rec.Sample(i)
		if err != nil {
			return err
		}
		v, ok := s.Get("SAC")
		if !ok {
			continue
		}
		sac, err := parseSAC(v.String())
		if err != nil {
			return err
		}
		var outs []string
		for j := 2; j+1 < len(sac); j += 2 {
			plus, minus := sac[j], sac[j+1]
			switch {
			case plus+minus == 0:
				outs = append(outs, "0.0")
			case plus > minus:
				outs = append(outs, round(float64(minus)/float64(plus+minus), 3))
			default:
				outs = append(outs, round(float64(plus)/float64(plus+minus), 3))
			}
		}
		if len(outs) == 0 {
			if err := rec.AddFormat(i, "CG_SBPCT", "."); err != nil {
				return err
			}
		} else if err := rec.AddFormat(i, "CG_SBPCT", strings.Join(outs, ",")); err != nil {
			return err
		}
	}
	return nil
}

// FisherStrandBias adds FORMAT CG_FSB (Phred-scaled Fisher strand-bias p-value
// per alt, against a theoretical 50/50 split).
type FisherStrandBias struct {
	closeNoop
	fisher *stats.FisherExact
}

// NewFisherSB returns a Fisher strand-bias annotator (--fisher-sb).
func NewFisherSB() *FisherStrandBias {
	return &FisherStrandBias{fisher: stats.NewFisherExact()}
}

// SetupHeader requires SAC and declares CG_FSB.
func (a *FisherStrandBias) SetupHeader(h *vcf.VcfHeader) error {
	if err := requireSAC(h); err != nil {
		return err
	}
	h.AddFormat(formatDef("CG_FSB", "A", "Float", "Sample-based Fisher Strand Bias for alt alleles (Phred-scale)"))
	return nil
}

// Annotate computes per-alt Fisher strand bias for each sample with SAC.
func (a *FisherStrandBias) Annotate(rec *vcf.VcfRecord) error {
	for i := 0; i < rec.NumSamples(); i++ {
		s, err := rec.Sample(i)
		if err != nil {
			return err
		}
		v, ok := s.Get("SAC")
		if !ok {
			continue
		}
		sac, err := parseSAC(v.String())
		if err != nil {
			return err
		}
		var outs []string
		for j := 2; j+1 < len(sac); j += 2 {
			plus, minus := sac[j], sac[j+1]
			total := plus + minus
			half := total / 2 // rounds down for odd totals
			p := a.fisher.TwoTailedPvalue(half, half, plus, minus)
			outs = append(outs, round(stats.Phred(p), 3))
		}
		if len(outs) == 0 {
			if err := rec.AddFormat(i, "CG_FSB", "."); err != nil {
				return err
			}
		} else if err := rec.AddFormat(i, "CG_FSB", strings.Join(outs, ",")); err != nil {
			return err
		}
	}
	return nil
}

// CopyNumberLogRatio adds INFO CG_CNLR, the log2 ratio of somatic vs germline
// allelic depth (AD) at the variant, optionally normalized by total counts.
type CopyNumberLogRatio struct {
	closeNoop
	germlineSample string
	somaticSample  string
	germlineCount  int64
	somaticCount   int64
	hasTotals      bool

	germTotalLog float64
	somTotalLog  float64
	germlineIdx  int
	somaticIdx   int
}

// NewCopyLogRatio returns a copy-number log-ratio annotator (--copy-logratio).
// When somaticCount/germlineCount are > 0 the ratio is normalized by their log2.
func NewCopyLogRatio(somaticSample, germlineSample string, somaticCount, germlineCount int64) *CopyNumberLogRatio {
	c := &CopyNumberLogRatio{
		germlineSample: germlineSample,
		somaticSample:  somaticSample,
		germlineCount:  germlineCount,
		somaticCount:   somaticCount,
		hasTotals:      somaticCount > 0 && germlineCount > 0,
	}
	if germlineCount >= 0 && somaticCount >= 0 && (germlineCount > 0 || somaticCount > 0) {
		c.germTotalLog = stats.Log2(float64(germlineCount))
		c.somTotalLog = stats.Log2(float64(somaticCount))
	}
	return c
}

// SetupHeader requires AD, resolves the two samples, and declares CG_CNLR.
func (a *CopyNumberLogRatio) SetupHeader(h *vcf.VcfHeader) error {
	d, ok := h.FormatDef("AD")
	if !ok || d.Number != "R" || d.Type != "Integer" {
		return fmt.Errorf("annotate: \"AD\" FORMAT annotation missing")
	}
	if a.hasTotals {
		h.AddInfo(infoDef("CG_CNLR", "1", "Float", fmt.Sprintf("Copy number (log2-ratio); Germline-total:%d, Somatic-total:%d", a.germlineCount, a.somaticCount)))
	} else {
		h.AddInfo(infoDef("CG_CNLR", "1", "Float", "Copy number (log2-ratio)"))
	}
	a.germlineIdx = h.SampleIndex(a.germlineSample)
	a.somaticIdx = h.SampleIndex(a.somaticSample)
	if a.germlineIdx < 0 {
		return fmt.Errorf("annotate: can't find germline sample: %s", a.germlineSample)
	}
	if a.somaticIdx < 0 {
		return fmt.Errorf("annotate: can't find somatic sample: %s", a.somaticSample)
	}
	return nil
}

// Annotate computes CG_CNLR from the two samples' AD sums.
func (a *CopyNumberLogRatio) Annotate(rec *vcf.VcfRecord) error {
	germ, err := rec.Sample(a.germlineIdx)
	if err != nil {
		return err
	}
	som, err := rec.Sample(a.somaticIdx)
	if err != nil {
		return err
	}
	gv, ok := germ.Get("AD")
	if !ok {
		return nil
	}
	sv, ok := som.Get("AD")
	if !ok {
		return nil
	}
	germAcc, err := gv.FloatFor("sum")
	if err != nil {
		return err
	}
	somAcc, err := sv.FloatFor("sum")
	if err != nil {
		return err
	}
	germLog := stats.Log2(germAcc)
	somLog := stats.Log2(somAcc)
	if a.hasTotals {
		rec.AddInfo("CG_CNLR", round((somLog-a.somTotalLog)-(germLog-a.germTotalLog), 6))
	} else {
		rec.AddInfo("CG_CNLR", round(somLog-germLog, 6))
	}
	return nil
}
