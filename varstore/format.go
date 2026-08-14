package varstore

import (
	"fmt"
	"sort"
	"strings"

	"github.com/parquet-go/parquet-go"
)

// Capturing per-sample FORMAT fields onto the ALT calls.
//
// The counterpart to captured INFO, and a better fit for the shape of the data:
// a call row IS one sample at one ALT, so a Number=A field maps to it exactly
// and a Number=1 field simply repeats down the ALTs of a multiallelic record.
//
// WHAT IT COSTS, which is the difference from INFO and the reason this defaults
// to nothing. calls.parquet is the large member -- 4.7 million rows against
// 50,000 sites on a 500-sample fixture, a ratio of 94 -- so a column here costs
// roughly a hundred times what the same column costs on the sites catalog. The
// whole member compresses to about four bytes per call across its twelve
// existing columns, because sample and chromosome are dictionary encoded and
// positions delta encoded; one added integer column is 15-35% on top of that,
// depending on how well the values compress. Worth it for something a query
// will actually use, and not worth it by default.
//
// AND ONLY FOR CARRIERS. A 0/0 is never stored, so a captured FORMAT value
// exists for ALT calls and nowhere else. Depth for reference calls comes from
// the callable runs, which carry MinDP per depth band -- so "DP for a carrier"
// and "DP for a non-carrier" are answered by two different members, by design.

// FormatField is one captured FORMAT field.
type FormatField struct {
	// Name is the VCF FORMAT key, e.g. "PID".
	Name string `json:"name"`
	// Column is where it lands in calls.parquet, e.g. "fmt_pid".
	Column string   `json:"column"`
	Type   InfoType `json:"type"`
	// Number is the VCF header's Number: "1" or "A".
	Number string `json:"number"`
}

// FormatPrefix namespaces captured columns, for the same reason InfoPrefix does:
// a call row already carries dp, gq, ad_ref and ad_alt, and a capture landing in
// one of those would overwrite a value the store computed or normalised.
const FormatPrefix = "fmt_"

// FormatColumn is the calls.parquet column a FORMAT key is captured into.
func FormatColumn(key string) string {
	return FormatPrefix + strings.ToLower(key)
}

// ValidateFormat checks captured FORMAT fields before anything is written.
func ValidateFormat(fields []FormatField) error {
	reserved := map[string]bool{}
	for _, f := range parquet.SchemaOf(Call{}).Fields() {
		reserved[f.Name()] = true
	}
	seen := map[string]string{}
	for _, f := range fields {
		if f.Name == "" {
			return fmt.Errorf("format field with no name")
		}
		if !infoKeyOK.MatchString(f.Name) {
			return fmt.Errorf("format %s: not a valid VCF FORMAT key", f.Name)
		}
		switch f.Type {
		case InfoInteger, InfoFloat, InfoString:
		case InfoFlag:
			// The VCF spec has no FORMAT flags: a per-sample field with no value
			// cannot be told from one that is absent.
			return fmt.Errorf("format %s: FORMAT fields cannot be flags", f.Name)
		default:
			return fmt.Errorf("format %s: unsupported Type=%s", f.Name, f.Type)
		}

		switch f.Number {
		case "1", "A":
		case "R":
			// AD is the Number=R field anybody actually wants, and the store
			// already splits it into ad_ref and ad_alt. Naming the one that
			// exists is more useful than a second, differently-shaped copy.
			return fmt.Errorf(
				"format %s: Number=R needs one value per allele including the reference, "+
					"which a call row has no room for -- AD is already stored as ad_ref and ad_alt",
				f.Name)
		default:
			return fmt.Errorf(
				"format %s: only Number=1 and Number=A can be captured, header says Number=%s "+
					"(a call row is one sample at one ALT, so there is nowhere for the other values to go)",
				f.Name, f.Number)
		}

		col := f.Column
		if col == "" {
			col = FormatColumn(f.Name)
		}
		if reserved[col] {
			return fmt.Errorf("format %s: column %s is one of the call's own", f.Name, col)
		}
		if prev, dup := seen[col]; dup {
			return fmt.Errorf("format %s and %s both want column %s", prev, f.Name, col)
		}
		seen[col] = f.Name
	}
	return nil
}

func normalizeFormat(fields []FormatField) []FormatField {
	if len(fields) == 0 {
		return nil
	}
	out := make([]FormatField, len(fields))
	copy(out, fields)
	for i := range out {
		if out[i].Column == "" {
			out[i].Column = FormatColumn(out[i].Name)
		}
	}
	sort.Slice(out, func(i, j int) bool { return out[i].Column < out[j].Column })
	return out
}

// callSchemaWith builds the calls schema carrying these captured columns.
func callSchemaWith(fields []FormatField) *parquet.Schema {
	g := parquet.Group{}
	for _, f := range parquet.SchemaOf(Call{}).Fields() {
		g[f.Name()] = f
	}
	for _, f := range fields {
		g[f.Column] = infoNode(InfoField{Type: f.Type})
	}
	return parquet.NewSchema("Call", g)
}

// FormatRow reads captured values out of one call row, over the scan's reusable
// buffer rather than a map per row.
type FormatRow struct {
	row  parquet.Row
	cols map[string]int
}

func (r FormatRow) Present(column string) bool {
	i, ok := r.cols[column]
	return ok && i < len(r.row) && !r.row[i].IsNull()
}

func (r FormatRow) Int(column string) (int32, bool) {
	v, ok := r.value(column)
	if !ok {
		return 0, false
	}
	return v.Int32(), true
}

func (r FormatRow) Float(column string) (float64, bool) {
	v, ok := r.value(column)
	if !ok {
		return 0, false
	}
	return v.Double(), true
}

func (r FormatRow) String(column string) (string, bool) {
	v, ok := r.value(column)
	if !ok {
		return "", false
	}
	return v.Clone().String(), true
}

func (r FormatRow) value(column string) (parquet.Value, bool) {
	i, ok := r.cols[column]
	if !ok || i >= len(r.row) || r.row[i].IsNull() {
		return parquet.Value{}, false
	}
	return r.row[i], true
}

// formatAsInfo reuses the INFO scratch machinery, whose only interest in a
// field is its type.
func formatAsInfo(fields []FormatField) []InfoField {
	out := make([]InfoField, len(fields))
	for i, f := range fields {
		out[i] = InfoField{Name: f.Name, Column: f.Column, Type: f.Type, Number: f.Number}
	}
	return out
}
