package varstore

import (
	"fmt"
	"math"
	"regexp"
	"sort"
	"strings"

	"github.com/parquet-go/parquet-go"
)

// Capturing INFO fields from the source VCF into the sites catalog.
//
// WHY THIS IS NOT ANNOTATION, which is the question it always raises. An INFO
// field is a property of this variant *in this file*: imputation R2 depends on
// which reference panel was used, VQSLOD on which model was trained, and a
// different VCF of the same cohort gives different numbers at the same locus.
// Annotation answers "what does this variant do", which is a property of the
// variant itself and belongs to an annotation service. This is a property of the
// OBSERVATION, which is exactly what a store records.
//
// WHY DECLARED COLUMNS RATHER THAN A BLOB. Serializing INFO into one string or
// JSON column is the flexible-looking option and it gives up both halves of why
// this is Parquet at all: nothing can prune a row group on a field inside a
// blob, and a typed float column compresses far better than its decimal text.
// Declaring the fields costs the caller one flag and keeps both.
//
// WHY THE SOURCE KEY IS PRESERVED. Imputation quality is not one field --
// minimac writes R2, IMPUTE2 writes INFO, Beagle writes DR2 and AR2 -- and the
// estimators are not numerically interchangeable. A store that flattened all of
// them into a column called "r2" would let two parts quote different estimators
// under one name, with nothing to catch it. So `info_dr2` and `info_r2` stay
// distinct, and a caller comparing across stores has to see that they differ.

// InfoType is a VCF INFO field's declared Type.
type InfoType string

const (
	InfoInteger InfoType = "Integer"
	InfoFloat   InfoType = "Float"
	InfoFlag    InfoType = "Flag"
	InfoString  InfoType = "String"
)

// InfoField is one captured INFO field, as declared by the source VCF's header.
//
// Number is kept verbatim because it is the reason a field is or is not
// capturable, and a reader asking "why is there no info_ac column" deserves the
// answer rather than silence.
type InfoField struct {
	// Name is the VCF INFO key, e.g. "R2".
	Name string `json:"name"`
	// Column is where it lands in catalog.parquet, e.g. "info_r2".
	Column string   `json:"column"`
	Type   InfoType `json:"type"`
	// Number is the VCF header's Number: "1" (one value per site) or "A" (one
	// per ALT). Flags carry "0".
	Number string `json:"number"`
}

// InfoPrefix namespaces captured columns.
//
// NOT COSMETIC. catalog.parquet already has `ac`, `an` and `n_called`, and those
// are RECOMPUTED over the samples in the store rather than copied from the
// source's INFO -- which is the whole reason they can be trusted after a subset
// or a merge. An unprefixed `--info AC` would overwrite the computed column with
// the source's own claim and silently undo that, with the result looking
// entirely normal. The prefix makes the collision impossible rather than
// forbidden.
const InfoPrefix = "info_"

var infoKeyOK = regexp.MustCompile(`^[A-Za-z_][0-9A-Za-z_.]*$`)

// InfoColumn is the catalog.parquet column an INFO key is captured into.
func InfoColumn(key string) string {
	return InfoPrefix + strings.ToLower(key)
}

// ValidateInfo checks a set of captured fields before a byte is written.
//
// Every refusal here NAMES THE FIELD and says what about it cannot be stored,
// because the caller's remedy is always to drop that one field rather than to
// abandon the conversion.
func ValidateInfo(fields []InfoField) error {
	reserved := map[string]bool{}
	for _, f := range parquet.SchemaOf(Site{}).Fields() {
		reserved[f.Name()] = true
	}

	seenCol := map[string]string{}
	for _, f := range fields {
		if f.Name == "" {
			return fmt.Errorf("info field with no name")
		}
		if !infoKeyOK.MatchString(f.Name) {
			return fmt.Errorf("info %s: not a valid VCF INFO key", f.Name)
		}
		switch f.Type {
		case InfoInteger, InfoFloat, InfoFlag, InfoString:
		default:
			return fmt.Errorf("info %s: unsupported Type=%s", f.Name, f.Type)
		}

		// Number is what decides capturability, and the reason is the shape of
		// catalog.parquet rather than taste. A site row is one (chrom,pos,ref,alt)
		// -- one ALT -- so Number=A maps to it exactly, and Number=1 repeats
		// across a multi-allelic record's rows, which is correct. Number=R
		// carries a reference value with nowhere to go, and Number=G and
		// Number=. are variable-length: keeping "the first value" of either
		// would be a number that means something different per row, which is
		// worse than not having it.
		switch f.Type {
		case InfoFlag:
			if f.Number != "0" && f.Number != "" {
				return fmt.Errorf("info %s: Flag must be Number=0, header says Number=%s", f.Name, f.Number)
			}
		default:
			if f.Number != "1" && f.Number != "A" {
				return fmt.Errorf(
					"info %s: only Number=1 and Number=A can be captured, header says Number=%s "+
						"(a site row is one ALT, so there is nowhere for the other values to go)",
					f.Name, f.Number)
			}
		}

		col := f.Column
		if col == "" {
			col = InfoColumn(f.Name)
		}
		if reserved[col] {
			return fmt.Errorf("info %s: column %s is one of the store's own", f.Name, col)
		}
		if prev, dup := seenCol[col]; dup {
			return fmt.Errorf("info %s and %s both want column %s", prev, f.Name, col)
		}
		seenCol[col] = f.Name
	}
	return nil
}

// normalizeInfo fills in the derived column names and sorts for a stable
// manifest, after validation has established the set is coherent.
func normalizeInfo(fields []InfoField) []InfoField {
	if len(fields) == 0 {
		return nil
	}
	out := make([]InfoField, len(fields))
	copy(out, fields)
	for i := range out {
		if out[i].Column == "" {
			out[i].Column = InfoColumn(out[i].Name)
		}
		if out[i].Type == InfoFlag && out[i].Number == "" {
			out[i].Number = "0"
		}
	}
	sort.Slice(out, func(i, j int) bool { return out[i].Column < out[j].Column })
	return out
}

// infoNode is the parquet leaf an INFO field is stored as.
//
// OPTIONAL FOR EVERYTHING BUT A FLAG, because absent and zero are different
// claims. A site where the caller's imputation program emitted no R2 is not a
// site with R2 = 0, and a required column has no way to say so. A Flag is the
// one case where absence genuinely means false, which is what the VCF spec says
// a Flag is.
func infoNode(f InfoField) parquet.Node {
	var leaf parquet.Node
	switch f.Type {
	case InfoInteger:
		leaf = parquet.Int(32)
	case InfoFloat:
		// Double rather than Float even though VCF declares single precision:
		// parsing "0.87" straight to float64 round-trips the decimal the file
		// actually contained, where going via float32 stores 0.8700000047.
		leaf = parquet.Leaf(parquet.DoubleType)
	case InfoFlag:
		return parquet.Leaf(parquet.BooleanType)
	default:
		// Dictionary encoded: captured strings are overwhelmingly repeated
		// labels rather than free text.
		leaf = parquet.Encoded(parquet.String(), &parquet.RLEDictionary)
	}
	return parquet.Optional(leaf)
}

// siteSchemaWith builds the sites schema carrying these captured columns.
//
// Site's own fields come from its struct so the two can never drift; adding a
// column to Site needs no change here.
func siteSchemaWith(fields []InfoField) *parquet.Schema {
	g := parquet.Group{}
	for _, f := range parquet.SchemaOf(Site{}).Fields() {
		g[f.Name()] = f
	}
	for _, f := range fields {
		g[f.Column] = infoNode(f)
	}
	return parquet.NewSchema("Site", g)
}

// InfoRow reads captured values out of one site row.
//
// Backed by the scan's reusable buffer rather than a map per row: a whole-genome
// catalog is tens of millions of sites, and an allocation each would dominate
// the walk that pre-annotation does.
type InfoRow struct {
	row  parquet.Row
	cols map[string]int
}

// Present reports whether this row carried a value for the field.
//
// The distinction a null column exists to preserve: Value returns zero both for
// "absent" and for a genuine zero, and only this can tell them apart.
func (r InfoRow) Present(column string) bool {
	i, ok := r.cols[column]
	if !ok || i >= len(r.row) {
		return false
	}
	return !r.row[i].IsNull()
}

// Float returns a captured Float, and false when the row had no value.
func (r InfoRow) Float(column string) (float64, bool) {
	v, ok := r.value(column)
	if !ok {
		return 0, false
	}
	return v.Double(), true
}

// Int returns a captured Integer, and false when the row had no value.
func (r InfoRow) Int(column string) (int32, bool) {
	v, ok := r.value(column)
	if !ok {
		return 0, false
	}
	return v.Int32(), true
}

// Flag returns a captured Flag. Absence is false, which is what a Flag means.
func (r InfoRow) Flag(column string) bool {
	v, ok := r.value(column)
	return ok && v.Boolean()
}

// String returns a captured String, and false when the row had no value.
func (r InfoRow) String(column string) (string, bool) {
	v, ok := r.value(column)
	if !ok {
		return "", false
	}
	return v.Clone().String(), true
}

func (r InfoRow) value(column string) (parquet.Value, bool) {
	i, ok := r.cols[column]
	if !ok || i >= len(r.row) {
		return parquet.Value{}, false
	}
	if r.row[i].IsNull() {
		return parquet.Value{}, false
	}
	return r.row[i], true
}

// Depth banding for callable runs.
//
// THE PROBLEM IT SOLVES. A run means "called at or above the gate throughout",
// and recording the lowest depth inside it gives a Ref call a confidence it
// otherwise has none of. But one poorly covered base drags a whole run's bound
// down: a megabase at 60x containing a single site at 11x reports MinDP 11, and
// every reference call across that span inherits the worst moment in it.
//
// Banding is the standard answer, and it is what a gVCF does. Break the run
// when depth crosses a boundary, so each run spans a depth CLASS rather than an
// arbitrary stretch. Bounds stay tight, and run counts grow only where coverage
// is genuinely ragged -- for 30x WGS most of the genome sits inside one band and
// the extra breaks are few.
//
// The boundaries are recorded in the manifest, because two stores banded
// differently do not mean the same thing by a run and a consumer comparing
// their MinDP values across parts needs to know.

// DefaultDepthRatio is the ratio between successive depth bands a conversion
// uses unless given an explicit list.
//
// WHY GEOMETRIC. A band boundary is a query gate the store answers exactly: a
// run never straddles one, so a gate set on a boundary admits or refuses each
// run whole. A gate strictly inside a band is answered conservatively -- a run
// whose floor sits below the gate is refused even where the site the query asks
// about was deeper -- and the loss is bounded by the next boundary up. Spacing
// boundaries by a fixed ratio makes that bound one sentence, at every depth:
// the store never withholds a reference call at a site whose depth is at least
// ratio x the gate.
//
// The cost of a boundary is the runs it breaks, which happens wherever a
// sample's depth crosses it between consecutive sites. Depth noise grows with
// depth, so geometric spacing keeps the crossing rate roughly even across the
// range rather than piling boundaries up where coverage is high.
const DefaultDepthRatio = 1.5

// MaxDepthBand caps a generated band list. Boundaries above any real depth are
// free -- nothing crosses them -- so the cap only bounds the list's length.
const MaxDepthBand int32 = 1000

// GeometricDepthBands returns boundaries floor(start * ratio^k), strictly
// increasing and not above max. start is the conversion's --min-dp: below it
// nothing is callable, so a boundary there would never break a run.
//
// The list is recorded in the manifest as generated, never regenerated on read,
// so a reader does not depend on this function's floating-point behaviour
// staying the same.
func GeometricDepthBands(start int32, ratio float64, max int32) ([]int32, error) {
	if start < 1 {
		return nil, fmt.Errorf("depth bands must start at 1 or more, got %d", start)
	}
	if !(ratio > 1) || math.IsInf(ratio, 0) {
		return nil, fmt.Errorf("depth ratio must be greater than 1, got %v", ratio)
	}
	if max < start {
		return nil, fmt.Errorf("depth band cap %d is below the starting depth %d", max, start)
	}
	bands := []int32{start}
	for k := 1; ; k++ {
		v := math.Floor(float64(start) * math.Pow(ratio, float64(k)))
		if v > float64(max) {
			break
		}
		b := int32(v)
		// A ratio close to 1 at a low start floors to the same integer more
		// than once; a repeated boundary would mean nothing.
		if b > bands[len(bands)-1] {
			bands = append(bands, b)
		}
	}
	return bands, nil
}

// ValidateDepthBands checks an explicit, caller-supplied band list. Boundaries
// must be positive and strictly increasing: a repeat or a reversal would make
// "the band a depth falls in" ambiguous.
func ValidateDepthBands(bands []int32) error {
	for i, b := range bands {
		if b < 1 {
			return fmt.Errorf("depth band %d is not positive", b)
		}
		if i > 0 && b <= bands[i-1] {
			return fmt.Errorf("depth bands must be strictly increasing; %d follows %d", b, bands[i-1])
		}
	}
	return nil
}

// DepthBand returns the index of the band a depth falls in, or -1 when bands
// are not in use. Depths below the first boundary share band 0 with it, since
// they are already excluded by the gate.
func DepthBand(bands []int32, dp int32) int {
	if len(bands) == 0 {
		return -1
	}
	// The first boundary strictly above dp; the band is the one before it.
	i := sort.Search(len(bands), func(i int) bool { return bands[i] > dp })
	if i == 0 {
		return 0
	}
	return i - 1
}

// BandContaining returns the boundaries either side of a gate that falls
// strictly inside a band -- lo < gate < hi -- and ok=false when the gate sits on
// a boundary, below the first, or bands are not in use. hi is 0 when the gate is
// above the last boundary, where the band has no upper edge.
//
// It is what a query reports when its answer is a conservative subset rather
// than exact: reference calls are guaranteed only at sites with depth >= hi.
func BandContaining(bands []int32, gate int32) (lo, hi int32, ok bool) {
	if len(bands) == 0 || gate <= bands[0] {
		return 0, 0, false
	}
	i := DepthBand(bands, gate)
	if bands[i] == gate {
		return 0, 0, false
	}
	lo = bands[i]
	if i+1 < len(bands) {
		hi = bands[i+1]
	}
	return lo, hi, true
}
