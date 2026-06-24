// Package annotate provides a composable framework for annotating VCF variants:
// each [Annotator] reads a [vcf.VcfRecord] (or a bare tuple built with
// [vcf.NewRecord]) and writes INFO/FORMAT/ID fields back onto the same record.
// Annotators are reusable from any Go code, not just the cgio CLI.
//
// A [Pipeline] runs a record source through a set of annotators in order,
// mutating one record pointer in place (no copies). Most annotators implement
// the per-record [Annotator] interface; the few that need neighbor look-ahead
// (e.g. variant distance) implement [StreamWrapper] instead.
//
// Annotators split into two kinds by what they read. Locus annotators use only
// CHROM/POS/REF/ALT/INFO/ID (e.g. Indel, TsTv, AutoID, ConstantTag,
// VariantDistance) and work on any record, including bare tuples built with
// [vcf.NewRecord]. Sample-based annotators read and write per-sample FORMAT
// values (e.g. Dosage, VAF, MinorStrand, FisherStrandBias, CopyNumberLogRatio);
// on a record with no samples (such as a tuple-built record) they iterate zero
// samples and add nothing — a silent no-op.
package annotate

import (
	"github.com/compgenlab/hts/vcf"
)

// Annotator reads a variant record and writes annotations (INFO/FORMAT/ID) back
// onto it. SetupHeader declares the ##INFO/##FORMAT defs the annotator adds.
type Annotator interface {
	SetupHeader(h *vcf.VcfHeader) error
	Annotate(rec *vcf.VcfRecord) error
	Close() error
}

// StreamWrapper is implemented by annotators that need to control the record
// stream (look-ahead/buffering), such as nearest-variant distance. Wrap returns
// a new source that pulls from next and emits annotated records.
type StreamWrapper interface {
	SetupHeader(h *vcf.VcfHeader) error
	Wrap(next Source) Source
	Close() error
}

// Source is a pull-based record stream: it returns io.EOF when exhausted.
type Source func() (*vcf.VcfRecord, error)

// CoordAware is implemented by annotators that resolve their query coordinates
// from the record (optionally overridden by INFO fields). The CLI applies the
// global --alt-chrom/--alt-pos/--end-pos settings to every annotator that
// implements it.
type CoordAware interface {
	SetAltChrom(key string)
	SetAltPos(key string)
	SetEndPos(key string)
}

// closeNoop is embedded by annotators with no resources to release.
type closeNoop struct{}

// Close implements the no-op Close for embedding annotators.
func (closeNoop) Close() error { return nil }

// infoDef builds an ##INFO definition.
func infoDef(id, number, typ, desc string) *vcf.AnnotationDef {
	return &vcf.AnnotationDef{IsInfo: true, ID: id, Number: number, Type: typ, Description: desc}
}

// formatDef builds an ##FORMAT definition.
func formatDef(id, number, typ, desc string) *vcf.AnnotationDef {
	return &vcf.AnnotationDef{IsInfo: false, ID: id, Number: number, Type: typ, Description: desc}
}

// infoDefSrc builds an ##INFO definition with a Source (the originating file).
func infoDefSrc(id, number, typ, desc, source string) *vcf.AnnotationDef {
	d := infoDef(id, number, typ, desc)
	d.Source = source
	return d
}

// formatDefSrc builds an ##FORMAT definition with a Source (the originating file).
func formatDefSrc(id, number, typ, desc, source string) *vcf.AnnotationDef {
	d := formatDef(id, number, typ, desc)
	d.Source = source
	return d
}

// headerSetuper is the common header step of Annotator and StreamWrapper.
type headerSetuper interface {
	SetupHeader(h *vcf.VcfHeader) error
}

// Pipeline composes a set of annotators over a record source. Per-record
// annotators run in the order added; stream wrappers (look-ahead) wrap the
// source. SetupHeaders runs every annotator's header step in add order.
type Pipeline struct {
	perRecord []Annotator
	streams   []StreamWrapper
	order     []headerSetuper
}

// NewPipeline returns an empty pipeline.
func NewPipeline() *Pipeline { return &Pipeline{} }

// Add appends a per-record annotator.
func (p *Pipeline) Add(a Annotator) {
	p.perRecord = append(p.perRecord, a)
	p.order = append(p.order, a)
}

// AddStream appends a look-ahead (stream) annotator.
func (p *Pipeline) AddStream(s StreamWrapper) {
	p.streams = append(p.streams, s)
	p.order = append(p.order, s)
}

// Len reports how many annotators (of both kinds) have been added.
func (p *Pipeline) Len() int { return len(p.order) }

// SetupHeaders runs each annotator's header step, in the order they were added.
func (p *Pipeline) SetupHeaders(h *vcf.VcfHeader) error {
	for _, x := range p.order {
		if err := x.SetupHeader(h); err != nil {
			return err
		}
	}
	return nil
}

// Build wires the annotators over source and returns the composed source. Each
// record flows through every per-record annotator (mutated in place), then
// through the stream wrappers.
func (p *Pipeline) Build(source Source) Source {
	next := source
	if len(p.perRecord) > 0 {
		inner := next
		anns := p.perRecord
		next = func() (*vcf.VcfRecord, error) {
			rec, err := inner()
			if err != nil {
				return nil, err
			}
			for _, a := range anns {
				if err := a.Annotate(rec); err != nil {
					return nil, err
				}
			}
			return rec, nil
		}
	}
	for _, s := range p.streams {
		next = s.Wrap(next)
	}
	return next
}

// Close releases every annotator, returning the first error.
func (p *Pipeline) Close() error {
	var firstErr error
	for _, a := range p.perRecord {
		if err := a.Close(); err != nil && firstErr == nil {
			firstErr = err
		}
	}
	for _, s := range p.streams {
		if err := s.Close(); err != nil && firstErr == nil {
			firstErr = err
		}
	}
	return firstErr
}
