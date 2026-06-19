package bed

import (
	"strconv"
	"strings"
)

// Strand represents the orientation of a BED record.
type Strand string

const (
	// StrandPlus is the forward ("+") strand.
	StrandPlus Strand = "+"
	// StrandMinus is the reverse ("-") strand.
	StrandMinus Strand = "-"
	// StrandNone is an unspecified strand. It is never emitted directly: the
	// writer forces it to "+" when a strand column is required, matching the
	// reference ngsutilsj behavior.
	StrandNone Strand = "."
)

// ParseStrand maps a BED strand field to a Strand. "+" and "-" map to
// [StrandPlus] and [StrandMinus]; anything else (including ".") maps to
// [StrandNone].
func ParseStrand(s string) Strand {
	switch s {
	case "+":
		return StrandPlus
	case "-":
		return StrandMinus
	default:
		return StrandNone
	}
}

// BedRecord is a single BED interval. Coordinates are 0-based half-open
// [Start, End). Name, Score, and Strand are the BED6 fields; Extras holds any
// columns beyond the sixth as opaque strings. When HasName is false the record
// is a bare BED3 interval and the name/score/strand fields are suppressed on
// auto-mode output.
//
// Columns records how many tab-separated columns the source line had (3..N).
// HasName is true even for a BED3 line, so Columns is the way to tell whether a
// record actually carried a strand (Columns >= 6).
type BedRecord struct {
	Ref     string
	Start   int // 0-based
	End     int // 0-based, exclusive
	Name    string
	HasName bool
	Score   float64
	Strand  Strand
	Extras  []string
	Columns int
}

// NewBed3 returns a bare BED3 record (no name/score/strand, no extras).
func NewBed3(ref string, start, end int) *BedRecord {
	return &BedRecord{Ref: ref, Start: start, End: end, HasName: false, Strand: StrandNone, Columns: 3}
}

// NewBed6 returns a full BED6 record with no extra columns.
func NewBed6(ref string, start, end int, name string, score float64, strand Strand) *BedRecord {
	return &BedRecord{
		Ref:     ref,
		Start:   start,
		End:     end,
		Name:    name,
		HasName: true,
		Score:   score,
		Strand:  strand,
		Columns: 6,
	}
}

// Length returns the span length (End - Start).
func (r *BedRecord) Length() int {
	return r.End - r.Start
}

// Extend5 returns a copy of the record extended by n bases on the 5' end
// (strand-aware). On the plus or unspecified strand this moves Start left; on
// the minus strand it moves End right. A negative n shrinks the region. When
// Start moves it is clamped at 0 (matching ngsutilsj GenomeSpan.extend5).
func (r *BedRecord) Extend5(n int) *BedRecord {
	out := r.clone()
	if r.Strand == StrandMinus {
		out.End += n
	} else {
		out.Start -= n
		if out.Start < 0 {
			out.Start = 0
		}
	}
	return out
}

// Extend3 returns a copy of the record extended by n bases on the 3' end
// (strand-aware). On the plus or unspecified strand this moves End right; on
// the minus strand it moves Start left. A negative n shrinks the region. When
// Start moves it is clamped at 0 (matching ngsutilsj GenomeSpan.extend3).
func (r *BedRecord) Extend3(n int) *BedRecord {
	out := r.clone()
	if r.Strand == StrandMinus {
		out.Start -= n
		if out.Start < 0 {
			out.Start = 0
		}
	} else {
		out.End += n
	}
	return out
}

func (r *BedRecord) clone() *BedRecord {
	cp := *r
	if r.Extras != nil {
		cp.Extras = append([]string(nil), r.Extras...)
	}
	return &cp
}

// formatScoreInt renders a score coerced to an integer, truncating toward zero
// (matching a Java (int) cast).
func formatScoreInt(score float64) string {
	return strconv.Itoa(int(score))
}

// formatScoreFloat renders a score as a plain decimal with the shortest
// round-trip representation, dropping a trailing ".0" for integral values
// (e.g. 1.5 -> "1.5", 1.0 -> "1", 0 -> "0"). This matches the reference
// ngsutilsj output for all non-scientific magnitudes; it does not reproduce
// Java Double.toString scientific notation for very large or very small values.
func formatScoreFloat(score float64) string {
	s := strconv.FormatFloat(score, 'f', -1, 64)
	return strings.TrimSuffix(s, ".0")
}

// parseBedLine parses a single BED line shared by the streaming and indexed
// readers. It returns ok=false for lines that should be skipped (blank,
// comment, or fewer than 3 columns). A non-nil error indicates a malformed
// coordinate or score and aborts iteration, matching the reference behavior.
func parseBedLine(line string, ignoreStrand bool) (rec *BedRecord, ok bool, err error) {
	line = strings.TrimRight(line, " \t\r\n")
	if line == "" || line[0] == '#' {
		return nil, false, nil
	}

	cols := strings.Split(line, "\t")
	if len(cols) < 3 {
		return nil, false, nil
	}

	start, err := strconv.Atoi(cols[1])
	if err != nil {
		return nil, false, err
	}
	end, err := strconv.Atoi(cols[2])
	if err != nil {
		return nil, false, err
	}

	out := &BedRecord{
		Ref:     cols[0],
		Start:   start,
		End:     end,
		HasName: true,
		Strand:  StrandNone,
		Columns: len(cols),
	}

	if len(cols) > 3 {
		out.Name = cols[3]
	}
	if len(cols) > 4 {
		score, err := strconv.ParseFloat(cols[4], 64)
		if err != nil {
			return nil, false, err
		}
		out.Score = score
	}
	if !ignoreStrand && len(cols) > 5 {
		out.Strand = ParseStrand(cols[5])
	}
	if len(cols) > 6 {
		out.Extras = append([]string(nil), cols[6:]...)
	}

	return out, true, nil
}
