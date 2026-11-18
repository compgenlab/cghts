package gtstore

import (
	"errors"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/compgenlab/cghts/vcf"
	"github.com/parquet-go/parquet-go"
)

// Captured INFO and FORMAT values in queries: filtering sites on an INFO
// field's presence (Query.ExcludeInfo / RequireInfo) and attaching values to
// the calls a query returns (Query.WithInfo / WithFormat).
//
// Both backends must give the same answer, and a VCF holds text where a store
// holds typed columns, so a value is read the way conversion captured it --
// InfoValue and FormatValue are that reading, shared with the converter -- and
// rendered through one function, ValueText. "0.870" in a VCF and 0.87 in a
// store therefore both come back as "0.87".

// ErrNotCaptured reports a query naming an INFO or FORMAT field the source
// does not have: for a store, one conversion did not capture; for a VCF, one
// its header does not declare. An error rather than "absent everywhere": a
// filter on a field that was never recorded would silently keep or drop every
// site.
var ErrNotCaptured = errors.New("field not captured")

// CallExtra carries the captured values a query asked for, rendered as text.
// A key is absent when the site or sample had no value. Nil on every call when
// the query asked for nothing, which is why it is a pointer: Call stays
// comparable, and a caller that never asks pays nothing.
type CallExtra struct {
	// Info is keyed by VCF INFO key. Site-level: every call at a site, ALT or
	// reference, carries the same values.
	Info map[string]string
	// Format is keyed by VCF FORMAT key. Per call. A reference call
	// reconstructed from a store has none -- a store keeps no FORMAT for a
	// genotype it never wrote down -- where the same call from a VCF has the
	// sample's recorded values, the same asymmetry DP/AD/GQ already have.
	Format map[string]string
}

// InfoValue reads one INFO field of a record for the ALT at altIdx (0-based),
// exactly as conversion captures it: a Flag is present or absent; a Number=A
// field takes that ALT's element; "." and an empty or unparseable value are
// absent. The result is int32, float64, string or bool by f.Type.
func InfoValue(rec *vcf.VcfRecord, f InfoField, altIdx int) (any, bool) {
	val, ok := rec.InfoValue(f.Name)
	if !ok {
		return nil, false
	}
	if f.Type == InfoFlag {
		return true, true
	}
	if val.IsMissing() {
		return nil, false
	}
	return parseCaptured(f.Type, f.Number, val.String(), altIdx)
}

// FormatValue is InfoValue for one sample's FORMAT field: col is the sample's
// column in the record, altIdx the 0-based ALT for a Number=A field.
func FormatValue(rec *vcf.VcfRecord, col int, f FormatField, altIdx int) (any, bool) {
	k := rec.FormatIndex(f.Name)
	if k < 0 {
		return nil, false
	}
	return parseCaptured(f.Type, f.Number, rec.SampleField(col, k), altIdx)
}

func parseCaptured(t InfoType, number, text string, altIdx int) (any, bool) {
	if text == "" || text == "." {
		return nil, false
	}
	if number == "A" {
		parts := strings.Split(text, ",")
		if altIdx < 0 || altIdx >= len(parts) {
			// Fewer values than ALTs: another allele's number would mean
			// something else.
			return nil, false
		}
		text = parts[altIdx]
		if text == "" || text == "." {
			return nil, false
		}
	}
	switch t {
	case InfoInteger:
		n, err := strconv.ParseInt(text, 10, 32)
		if err != nil {
			return nil, false
		}
		return int32(n), true
	case InfoFloat:
		x, err := strconv.ParseFloat(text, 64)
		if err != nil {
			return nil, false
		}
		return x, true
	case InfoFlag:
		return true, true
	default:
		return text, true
	}
}

// ValueText renders a captured value the one way both backends report it:
// integers in decimal, floats in their shortest round-tripping form, strings
// verbatim, and a present Flag as "1".
func ValueText(v any) string {
	switch x := v.(type) {
	case int32:
		return strconv.FormatInt(int64(x), 10)
	case float64:
		return strconv.FormatFloat(x, 'g', -1, 64)
	case bool:
		if x {
			return "1"
		}
		return "0"
	case string:
		return x
	}
	return fmt.Sprint(v)
}

// infoRowValue reads a captured INFO value from a store row. A Flag is
// written false when absent, so only true is presence.
func infoRowValue(r InfoRow, f InfoField) (any, bool) {
	switch f.Type {
	case InfoFlag:
		if r.Flag(f.Column) {
			return true, true
		}
		return nil, false
	case InfoInteger:
		return toAny(r.Int(f.Column))
	case InfoFloat:
		return toAny(r.Float(f.Column))
	default:
		return toAny(r.String(f.Column))
	}
}

func formatRowValue(r FormatRow, f FormatField) (any, bool) {
	switch f.Type {
	case InfoInteger:
		return toAny(r.Int(f.Column))
	case InfoFloat:
		return toAny(r.Float(f.Column))
	default:
		return toAny(r.String(f.Column))
	}
}

func toAny[T any](v T, ok bool) (any, bool) {
	if !ok {
		return nil, false
	}
	return v, true
}

// siteSelect is a query's INFO filter and value requests, resolved against the
// fields a source actually has.
type siteSelect struct {
	exclude, require []InfoField
	info             []InfoField
	format           []FormatField
}

func (s *siteSelect) filters() bool { return len(s.exclude)+len(s.require) > 0 }
func (s *siteSelect) wantsInfo() bool {
	return s.filters() || len(s.info) > 0
}

// excluded applies the presence filter to one site's values.
func (s *siteSelect) excluded(get func(InfoField) (any, bool)) bool {
	for _, f := range s.exclude {
		if _, ok := get(f); ok {
			return true
		}
	}
	for _, f := range s.require {
		if _, ok := get(f); !ok {
			return true
		}
	}
	return false
}

// infoText renders the requested INFO values for one site; nil when none was
// requested, an empty map when requested but absent.
func (s *siteSelect) infoText(get func(InfoField) (any, bool)) map[string]string {
	if len(s.info) == 0 {
		return nil
	}
	out := make(map[string]string, len(s.info))
	for _, f := range s.info {
		if v, ok := get(f); ok {
			out[f.Name] = ValueText(v)
		}
	}
	return out
}

func (s *siteSelect) formatText(get func(FormatField) (any, bool)) map[string]string {
	if len(s.format) == 0 {
		return nil
	}
	out := make(map[string]string, len(s.format))
	for _, f := range s.format {
		if v, ok := get(f); ok {
			out[f.Name] = ValueText(v)
		}
	}
	return out
}

// extra combines site INFO text and call FORMAT text into a CallExtra, or nil
// when the query asked for neither.
func (s *siteSelect) extra(info, format map[string]string) *CallExtra {
	if len(s.info) == 0 && len(s.format) == 0 {
		return nil
	}
	if format == nil && len(s.format) > 0 {
		format = map[string]string{}
	}
	if info == nil && len(s.info) > 0 {
		info = map[string]string{}
	}
	return &CallExtra{Info: info, Format: format}
}

// resolveSelect matches a query's keys against the INFO/FORMAT fields a source
// has. Keys are VCF keys, case-sensitive as VCF is. "*" in WithInfo/WithFormat
// means every captured field.
func resolveSelect(q Query, info []InfoField, format []FormatField) (*siteSelect, error) {
	byInfo := map[string]InfoField{}
	for _, f := range info {
		byInfo[f.Name] = f
	}
	byFmt := map[string]FormatField{}
	for _, f := range format {
		byFmt[f.Name] = f
	}
	pickInfo := func(what string, keys []string, star bool) ([]InfoField, error) {
		var out []InfoField
		seen := map[string]bool{}
		for _, k := range keys {
			if star && k == "*" {
				for _, f := range info {
					if !seen[f.Name] {
						seen[f.Name] = true
						out = append(out, f)
					}
				}
				continue
			}
			f, ok := byInfo[k]
			if !ok {
				return nil, fmt.Errorf("%s INFO %s: %w (captured: %s)", what, k, ErrNotCaptured, capturedInfoNames(info))
			}
			if !seen[k] {
				seen[k] = true
				out = append(out, f)
			}
		}
		return out, nil
	}
	var s siteSelect
	var err error
	if s.exclude, err = pickInfo("exclude", q.ExcludeInfo, false); err != nil {
		return nil, err
	}
	if s.require, err = pickInfo("require", q.RequireInfo, false); err != nil {
		return nil, err
	}
	if s.info, err = pickInfo("show", q.WithInfo, true); err != nil {
		return nil, err
	}
	seen := map[string]bool{}
	for _, k := range q.WithFormat {
		if k == "*" {
			for _, f := range format {
				if !seen[f.Name] {
					seen[f.Name] = true
					s.format = append(s.format, f)
				}
			}
			continue
		}
		f, ok := byFmt[k]
		if !ok {
			names := make([]string, 0, len(format))
			for _, f := range format {
				names = append(names, f.Name)
			}
			return nil, fmt.Errorf("show FORMAT %s: %w (captured: %s)", k, ErrNotCaptured, listOrNone(names))
		}
		if !seen[k] {
			seen[k] = true
			s.format = append(s.format, f)
		}
	}
	return &s, nil
}

func capturedInfoNames(fields []InfoField) string {
	names := make([]string, 0, len(fields))
	for _, f := range fields {
		names = append(names, f.Name)
	}
	return listOrNone(names)
}

func listOrNone(names []string) string {
	if len(names) == 0 {
		return "none"
	}
	sort.Strings(names)
	return strings.Join(names, ", ")
}

// HeaderInfoFields describes every INFO field a VCF header declares in the
// terms a store uses, so the VCF backend can resolve a query the way a store
// does. Number=R/G and fixed counts are kept: they cannot be captured, but a
// presence filter or a verbatim value still makes sense for them, so they are
// read as verbatim text.
func HeaderInfoFields(h *vcf.VcfHeader) []InfoField {
	var out []InfoField
	for _, id := range h.InfoIDs() {
		d, ok := h.InfoDef(id)
		if !ok {
			continue
		}
		out = append(out, headerField(d))
	}
	return out
}

// HeaderFormatFields is HeaderInfoFields for FORMAT, without the fields every
// call already carries in its own columns -- GT, DP, GQ and AD -- so "*" means
// the same thing against a VCF as against a store: everything beyond the
// standard columns. (A store's captured fields never include them unless a
// conversion asked for them by name.)
func HeaderFormatFields(h *vcf.VcfHeader) []FormatField {
	var out []FormatField
	for _, id := range h.FormatIDs() {
		switch id {
		case "GT", "DP", "GQ", "AD":
			continue
		}
		d, ok := h.FormatDef(id)
		if !ok {
			continue
		}
		f := headerField(d)
		out = append(out, FormatField{Name: f.Name, Column: FormatColumn(f.Name), Type: f.Type, Number: f.Number})
	}
	return out
}

func headerField(d *vcf.AnnotationDef) InfoField {
	typ, number := d.Type, d.Number
	switch {
	case typ == string(InfoFlag):
		number = "0"
	case number == "1" || number == "A" || number == NumberAny:
	default:
		// R, G or a fixed count: no per-ALT split this package trusts, so the
		// value is carried verbatim, as Number=. is.
		number = NumberAny
	}
	if typ != string(InfoInteger) && typ != string(InfoFloat) && typ != string(InfoFlag) {
		typ = string(InfoString)
	}
	return InfoFieldFor(d.ID, typ, number)
}

// scanRowsPruned is scanParquetPruned for dynamically read rows, so captured
// INFO/FORMAT columns -- which a typed Site or Call cannot hold -- come back
// with the rest. Pruning is the same: shards from the manifest, row groups from
// their statistics. The row views the scan's buffer: copy what you keep.
func scanRowsPruned(set *shardSet, keep scanFilter, fn func(row parquet.Row, cols map[string]int) bool) error {
	for _, sh := range set.shards {
		if sh.bounds && keep.shard != nil && !keep.shard(sh.info) {
			continue
		}
		pf, err := sh.m.parsed()
		if err != nil {
			return err
		}
		schema := pf.Schema()
		cols := map[string]int{}
		for _, path := range schema.Columns() {
			if lc, ok := schema.Lookup(path...); ok {
				cols[strings.Join(path, ".")] = lc.ColumnIndex
			}
		}
		stop, err := func() (bool, error) {
			r := parquet.NewGenericReader[any](pf, schema)
			defer r.Close()
			buf := make([]parquet.Row, 512)
			var offset int64
			for _, rg := range pf.RowGroups() {
				n := rg.NumRows()
				if keep.rowGroup != nil && !keep.rowGroup(rg) {
					offset += n
					continue
				}
				if err := r.SeekToRow(offset); err != nil {
					return false, err
				}
				for remaining := n; remaining > 0; {
					want := int64(len(buf))
					if remaining < want {
						want = remaining
					}
					got, err := r.ReadRows(buf[:want])
					for i := 0; i < got; i++ {
						if !fn(buf[i], cols) {
							return true, nil
						}
					}
					remaining -= int64(got)
					if err == io.EOF || got == 0 {
						break
					}
					if err != nil {
						return false, err
					}
				}
				offset += n
			}
			return false, nil
		}()
		if err != nil {
			return err
		}
		if stop {
			return nil
		}
	}
	return nil
}

// siteMeta is what a store says about one catalog site beyond its fixed
// columns: whether the query's INFO filter excludes it, and the INFO values
// the query asked to see.
type siteMeta struct {
	excluded bool
	info     map[string]string
}

// selectFor resolves a query's INFO/FORMAT requests against this volume.
//
// A volume with an empty catalog -- the blocks-only volume attach creates for a
// chromosome with no variant site -- has no site for a filter to act on and no
// value to show, so a field it never captured is not an error there; anywhere
// else it is, since filtering on a field that was never recorded would silently
// keep or drop everything.
func (s *ParquetVolume) selectFor(q Query) (*siteSelect, error) {
	if !q.wantsExtras() {
		return &siteSelect{}, nil
	}
	sel, err := resolveSelect(q, s.InfoFields(), s.FormatFields())
	if err != nil && errors.Is(err, ErrNotCaptured) && s.catalogEmpty() {
		return &siteSelect{}, nil
	}
	return sel, err
}

func (s *ParquetVolume) catalogEmpty() bool {
	if s.manifest == nil {
		return false
	}
	ti, ok := s.manifest.Tables[CatalogTable]
	return ok && ti.Rows == 0
}

// siteMetas reads the catalog sites a query touches with the INFO columns its
// filter and value requests need, keyed by canonical locus. Nil when the query
// asked for neither. Pruned like any catalog scan, so a locus query reads one
// row group, and bounded by the query, not the store.
func (s *ParquetVolume) siteMetas(p *plan, sel *siteSelect) (map[Locus]siteMeta, error) {
	if !sel.wantsInfo() {
		return nil, nil
	}
	if !s.hasSites {
		return nil, fmt.Errorf("%s is missing, so captured INFO cannot be read", CatalogPath(s.base))
	}
	out := map[Locus]siteMeta{}
	err := scanRowsPruned(s.sites, callsFilter(p.q), func(row parquet.Row, cols map[string]int) bool {
		site := siteFromRow(row, cols)
		loc := site.Locus()
		if !p.wantsRecord(loc, site.RefEnd) {
			return true
		}
		ir := InfoRow{row: row, cols: cols}
		get := func(f InfoField) (any, bool) { return infoRowValue(ir, f) }
		out[canonLocus(loc)] = siteMeta{excluded: sel.excluded(get), info: sel.infoText(get)}
		return true
	})
	return out, err
}

// scanCalls walks the calls a query touches, reading the captured FORMAT
// columns only when the query asked for them: the typed scan is the fast path,
// and calls is the table worth keeping fast.
func (s *ParquetVolume) scanCalls(keep scanFilter, sel *siteSelect, fn func(Call, map[string]string) bool) error {
	if len(sel.format) == 0 {
		return scanParquetPruned(s.calls, keep, func(c Call) bool { return fn(c, nil) })
	}
	return scanRowsPruned(s.calls, keep, func(row parquet.Row, cols map[string]int) bool {
		fr := FormatRow{row: row, cols: cols}
		return fn(callFromRow(row, cols), sel.formatText(func(f FormatField) (any, bool) { return formatRowValue(fr, f) }))
	})
}
