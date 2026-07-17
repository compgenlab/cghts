package gtf

import "github.com/compgenlab/cghts/bed"

// span is a genomic interval, 0-based half-open [start, end), with an optional
// strand. It ports the overlap/contains semantics of ngsutilsj GenomeSpan
// (annotation/GenomeSpan.java) — including its half-open boundary asymmetry
// ("there be off-by-one dragons here"). A negative start (or a query start of
// -1) denotes a whole-chromosome span.
type span struct {
	ref    string
	start  int
	end    int
	strand bed.Strand
}

// contains reports whether s fully contains q (q lies within s). Ports
// GenomeSpan.contains(...) (onlyWithin=true).
func (s span) contains(q span) bool { return s.match(q, true) }

// overlaps reports whether s overlaps q at all. Ports GenomeSpan.overlaps(...)
// (onlyWithin=false).
func (s span) overlaps(q span) bool { return s.match(q, false) }

// match is the shared core of contains/overlaps, a direct port of the protected
// GenomeSpan.contains(qref, qstart, qend, qstrand, onlyWithin). The <=/< split
// on the start/end boundaries is deliberate and must not be "tidied up".
func (s span) match(q span, onlyWithin bool) bool {
	if s.ref != q.ref {
		return false
	}
	// A strand of None matches either strand.
	if q.strand != bed.StrandNone && s.strand != bed.StrandNone && q.strand != s.strand {
		return false
	}
	if s.start < 0 || q.start == -1 {
		return true // whole-chromosome span: ref match is enough
	}

	startWithin := s.start <= q.start && q.start < s.end
	endWithin := s.start < q.end && q.end <= s.end

	if onlyWithin {
		return startWithin && endWithin
	}
	if startWithin || endWithin {
		return true
	}
	// The query region spans the entirety of s.
	return q.start <= s.start && s.start <= q.end && q.start <= s.end && s.end <= q.end
}
