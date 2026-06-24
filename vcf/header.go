package vcf

import (
	"fmt"
	"io"
	"strconv"
	"strings"
)

// AnnotationDef is a parsed ##INFO or ##FORMAT header definition.
type AnnotationDef struct {
	IsInfo      bool
	ID          string
	Number      string // "A", "R", "G", ".", or an integer
	Type        string // Integer, Float, Flag (INFO only), Character, String
	Description string
	Source      string
	Version     string
	Extras      map[string]string
	OrigLine    string // verbatim source line, used for round-trip output
}

// String renders the definition as a ##INFO/##FORMAT line. The original line is
// returned verbatim when available.
func (d *AnnotationDef) String() string {
	if d.OrigLine != "" {
		return d.OrigLine
	}
	var b strings.Builder
	if d.IsInfo {
		b.WriteString("##INFO=<")
	} else {
		b.WriteString("##FORMAT=<")
	}
	b.WriteString("ID=" + d.ID)
	b.WriteString(",Number=" + d.Number)
	b.WriteString(",Type=" + d.Type)
	b.WriteString(",Description=\"" + quoteString(d.Description) + "\"")
	if d.Source != "" {
		b.WriteString(",Source=\"" + quoteString(d.Source) + "\"")
	}
	if d.Version != "" {
		b.WriteString(",Version=\"" + quoteString(d.Version) + "\"")
	}
	b.WriteString(">")
	return b.String()
}

// FilterDef is a parsed ##FILTER header definition.
type FilterDef struct {
	ID          string
	Description string
	OrigLine    string
}

// String renders the definition as a ##FILTER line.
func (d *FilterDef) String() string {
	if d.OrigLine != "" {
		return d.OrigLine
	}
	return "##FILTER=<ID=" + d.ID + ",Description=\"" + quoteString(d.Description) + "\">"
}

// ContigDef is a parsed ##contig header definition. Length is -1 when absent.
type ContigDef struct {
	ID       string
	Length   int64
	OrigLine string
}

// String renders the definition as a ##contig line.
func (d *ContigDef) String() string {
	if d.OrigLine != "" {
		return d.OrigLine
	}
	if d.Length > 0 {
		return "##contig=<ID=" + d.ID + ",length=" + strconv.FormatInt(d.Length, 10) + ">"
	}
	return "##contig=<ID=" + d.ID + ">"
}

// AltDef is a parsed ##ALT header definition.
type AltDef struct {
	ID          string
	Description string
	OrigLine    string
}

// String renders the definition as a ##ALT line.
func (d *AltDef) String() string {
	if d.OrigLine != "" {
		return d.OrigLine
	}
	if d.Description != "" {
		return "##ALT=<ID=" + d.ID + ",Description=\"" + d.Description + "\">"
	}
	return "##ALT=<ID=" + d.ID + ">"
}

// VcfHeader holds the parsed metadata and sample names from a VCF header. It
// preserves the relative order of definitions for faithful round-trip output.
type VcfHeader struct {
	FileFormat string // the full "##fileformat=..." line

	infoOrder  []string
	infoDefs   map[string]*AnnotationDef
	formatOrd  []string
	formatDefs map[string]*AnnotationDef
	filterOrd  []string
	filterDefs map[string]*FilterDef
	contigOrd  []string
	contigDefs map[string]*ContigDef
	altOrd     []string
	altDefs    map[string]*AltDef

	otherLines []string // any other ## line, in original order

	samples   []string
	sampleIdx map[string]int
}

func newHeader() *VcfHeader {
	return &VcfHeader{
		FileFormat: "##fileformat=VCFv4.2",
		infoDefs:   map[string]*AnnotationDef{},
		formatDefs: map[string]*AnnotationDef{},
		filterDefs: map[string]*FilterDef{},
		contigDefs: map[string]*ContigDef{},
		altDefs:    map[string]*AltDef{},
		sampleIdx:  map[string]int{},
	}
}

// Samples returns the sample names in column order.
func (h *VcfHeader) Samples() []string { return h.samples }

// SampleIndex returns the 0-based index of a sample by name. If name is not a
// sample name but parses as a 1-based number ("1", "2", ...), that is used
// instead. It returns -1 when no sample matches.
func (h *VcfHeader) SampleIndex(name string) int {
	if i, ok := h.sampleIdx[name]; ok {
		return i
	}
	if n, err := strconv.Atoi(name); err == nil {
		return n - 1
	}
	return -1
}

// InfoDef returns the ##INFO definition for id, if present.
func (h *VcfHeader) InfoDef(id string) (*AnnotationDef, bool) {
	d, ok := h.infoDefs[id]
	return d, ok
}

// FormatDef returns the ##FORMAT definition for id, if present.
func (h *VcfHeader) FormatDef(id string) (*AnnotationDef, bool) {
	d, ok := h.formatDefs[id]
	return d, ok
}

// FilterDef returns the ##FILTER definition for id, if present.
func (h *VcfHeader) FilterDef(id string) (*FilterDef, bool) {
	d, ok := h.filterDefs[id]
	return d, ok
}

// ContigDef returns the ##contig definition for id, if present.
func (h *VcfHeader) ContigDef(id string) (*ContigDef, bool) {
	d, ok := h.contigDefs[id]
	return d, ok
}

// ContigNames returns the contig IDs in header order.
func (h *VcfHeader) ContigNames() []string { return h.contigOrd }

// InfoIDs returns the ##INFO field IDs in header order.
func (h *VcfHeader) InfoIDs() []string { return h.infoOrder }

// FormatIDs returns the ##FORMAT field IDs in header order.
func (h *VcfHeader) FormatIDs() []string { return h.formatOrd }

// MatchInfoIDs returns the ##INFO IDs matching the given glob (* and ?), in
// header order.
func (h *VcfHeader) MatchInfoIDs(glob string) []string {
	var out []string
	for _, id := range h.infoOrder {
		if globMatch(id, glob) {
			out = append(out, id)
		}
	}
	return out
}

// MatchFormatIDs returns the ##FORMAT IDs matching the given glob (* and ?), in
// header order.
func (h *VcfHeader) MatchFormatIDs(glob string) []string {
	var out []string
	for _, id := range h.formatOrd {
		if globMatch(id, glob) {
			out = append(out, id)
		}
	}
	return out
}

// MetaLines returns the metadata definition lines (INFO, FILTER, FORMAT, contig,
// ALT, then other preserved lines), excluding the fileformat line and the
// #CHROM column line. It matches ngsutilsj VCFHeader.write(out, false).
func (h *VcfHeader) MetaLines() []string {
	lines := h.Lines()
	// Lines() places the fileformat line first; drop it.
	if len(lines) > 0 && strings.HasPrefix(lines[0], "##fileformat=") {
		return lines[1:]
	}
	return lines
}

// AddInfo registers (or replaces) an ##INFO definition.
func (h *VcfHeader) AddInfo(d *AnnotationDef) {
	if _, ok := h.infoDefs[d.ID]; !ok {
		h.infoOrder = append(h.infoOrder, d.ID)
	}
	h.infoDefs[d.ID] = d
}

// AddFormat registers (or replaces) an ##FORMAT definition.
func (h *VcfHeader) AddFormat(d *AnnotationDef) {
	if _, ok := h.formatDefs[d.ID]; !ok {
		h.formatOrd = append(h.formatOrd, d.ID)
	}
	h.formatDefs[d.ID] = d
}

// AddFilter registers (or replaces) an ##FILTER definition.
func (h *VcfHeader) AddFilter(d *FilterDef) {
	if _, ok := h.filterDefs[d.ID]; !ok {
		h.filterOrd = append(h.filterOrd, d.ID)
	}
	h.filterDefs[d.ID] = d
}

// AddContig registers (or replaces) an ##contig definition.
func (h *VcfHeader) AddContig(d *ContigDef) {
	if _, ok := h.contigDefs[d.ID]; !ok {
		h.contigOrd = append(h.contigOrd, d.ID)
	}
	h.contigDefs[d.ID] = d
}

// RemoveContig drops an ##contig definition.
func (h *VcfHeader) RemoveContig(id string) {
	if _, ok := h.contigDefs[id]; !ok {
		return
	}
	delete(h.contigDefs, id)
	for i, c := range h.contigOrd {
		if c == id {
			h.contigOrd = append(h.contigOrd[:i], h.contigOrd[i+1:]...)
			break
		}
	}
}

// AddAlt registers (or replaces) an ##ALT definition.
func (h *VcfHeader) AddAlt(d *AltDef) {
	if _, ok := h.altDefs[d.ID]; !ok {
		h.altOrd = append(h.altOrd, d.ID)
	}
	h.altDefs[d.ID] = d
}

// AddLine appends a raw "##" metadata line, preserved verbatim on output.
func (h *VcfHeader) AddLine(line string) {
	h.otherLines = append(h.otherLines, line)
}

// SetFileDate replaces any existing ##fileDate line with one for the given date
// (in YYYYMMDD form).
func (h *VcfHeader) SetFileDate(date string) {
	kept := h.otherLines[:0]
	for _, l := range h.otherLines {
		if !strings.HasPrefix(l, "##fileDate=") {
			kept = append(kept, l)
		}
	}
	h.otherLines = append(kept, "##fileDate="+date)
}

// SetSamples replaces the sample list (and rebuilds the name index).
func (h *VcfHeader) SetSamples(samples []string) {
	h.samples = append([]string(nil), samples...)
	h.sampleIdx = make(map[string]int, len(samples))
	for i, s := range h.samples {
		h.sampleIdx[s] = i
	}
}

func (h *VcfHeader) addSample(name string) {
	h.sampleIdx[name] = len(h.samples)
	h.samples = append(h.samples, name)
}

// filterAllowed reports whether a FILTER code should be retained when parsing a
// record. It is a placeholder for future keep/remove-filter support and
// currently always returns true.
func (h *VcfHeader) filterAllowed(string) bool { return true }

// Lines returns the metadata lines (fileformat first, then INFO, FILTER,
// FORMAT, contig, ALT, then any other preserved lines), without the #CHROM
// column line. The ordering matches the reference ngsutilsj output.
func (h *VcfHeader) Lines() []string {
	out := []string{h.FileFormat}
	for _, id := range h.infoOrder {
		out = append(out, h.infoDefs[id].String())
	}
	for _, id := range h.filterOrd {
		out = append(out, h.filterDefs[id].String())
	}
	for _, id := range h.formatOrd {
		out = append(out, h.formatDefs[id].String())
	}
	for _, id := range h.contigOrd {
		out = append(out, h.contigDefs[id].String())
	}
	for _, id := range h.altOrd {
		out = append(out, h.altDefs[id].String())
	}
	out = append(out, h.otherLines...)
	return out
}

// ChromLine returns the #CHROM column line, including FORMAT and sample columns
// when samples are present.
func (h *VcfHeader) ChromLine() string {
	line := "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
	if len(h.samples) > 0 {
		line += "\tFORMAT\t" + strings.Join(h.samples, "\t")
	}
	return line
}

// WriteTo writes the full header (metadata lines plus the #CHROM line) to w.
func (h *VcfHeader) WriteTo(w io.Writer) (int64, error) {
	var total int64
	for _, line := range h.Lines() {
		n, err := io.WriteString(w, line+"\n")
		total += int64(n)
		if err != nil {
			return total, err
		}
	}
	n, err := io.WriteString(w, h.ChromLine()+"\n")
	total += int64(n)
	return total, err
}

// parseHeaderLines builds a VcfHeader from the collected "##" metadata lines and
// the "#CHROM" column line.
func parseHeaderLines(metaLines []string, chromLine string) (*VcfHeader, error) {
	h := newHeader()
	h.FileFormat = ""

	for _, line := range metaLines {
		switch {
		case strings.HasPrefix(line, "##fileformat="):
			h.FileFormat = line
		case strings.HasPrefix(line, "##INFO="):
			d, err := parseAnnotationDef(line)
			if err != nil {
				return nil, err
			}
			h.AddInfo(d)
		case strings.HasPrefix(line, "##FORMAT="):
			d, err := parseAnnotationDef(line)
			if err != nil {
				return nil, err
			}
			h.AddFormat(d)
		case strings.HasPrefix(line, "##FILTER="):
			d, err := parseFilterDef(line)
			if err != nil {
				return nil, err
			}
			h.AddFilter(d)
		case strings.HasPrefix(line, "##contig="):
			d, err := parseContigDef(line)
			if err != nil {
				return nil, err
			}
			h.AddContig(d)
		case strings.HasPrefix(line, "##ALT="):
			d, err := parseAltDef(line)
			if err != nil {
				return nil, err
			}
			h.AddAlt(d)
		default:
			h.otherLines = append(h.otherLines, line)
		}
	}

	if h.FileFormat == "" {
		h.FileFormat = "##fileformat=VCFv4.2"
	}

	cols := strings.Split(chromLine, "\t")
	if len(cols) > 9 {
		for _, s := range cols[9:] {
			h.addSample(s)
		}
	}
	return h, nil
}

func parseAnnotationDef(line string) (*AnnotationDef, error) {
	var isInfo bool
	var body string
	switch {
	case strings.HasPrefix(line, "##INFO=<") && strings.HasSuffix(line, ">"):
		isInfo = true
		body = line[len("##INFO=<") : len(line)-1]
	case strings.HasPrefix(line, "##FORMAT=<") && strings.HasSuffix(line, ">"):
		body = line[len("##FORMAT=<") : len(line)-1]
	default:
		return nil, fmt.Errorf("vcf: cannot parse header line: %s", line)
	}

	vals, err := parseQuotedLine(body)
	if err != nil {
		return nil, err
	}
	d := &AnnotationDef{
		IsInfo:      isInfo,
		ID:          vals["ID"],
		Number:      vals["Number"],
		Type:        vals["Type"],
		Description: vals["Description"],
		Source:      vals["Source"],
		Version:     vals["Version"],
		OrigLine:    line,
	}
	delete(vals, "ID")
	delete(vals, "Number")
	delete(vals, "Type")
	delete(vals, "Description")
	delete(vals, "Source")
	delete(vals, "Version")
	if len(vals) > 0 {
		d.Extras = vals
	}
	if d.ID == "" || d.Type == "" || d.Number == "" {
		return nil, fmt.Errorf("vcf: invalid INFO/FORMAT line: %s", line)
	}
	return d, nil
}

func parseFilterDef(line string) (*FilterDef, error) {
	if !strings.HasPrefix(line, "##FILTER=<") || !strings.HasSuffix(line, ">") {
		return nil, fmt.Errorf("vcf: cannot parse header line: %s", line)
	}
	vals, err := parseQuotedLine(line[len("##FILTER=<") : len(line)-1])
	if err != nil {
		return nil, err
	}
	return &FilterDef{ID: vals["ID"], Description: vals["Description"], OrigLine: line}, nil
}

func parseContigDef(line string) (*ContigDef, error) {
	if !strings.HasPrefix(line, "##contig=<") || !strings.HasSuffix(line, ">") {
		return nil, fmt.Errorf("vcf: cannot parse header line: %s", line)
	}
	vals, err := parseQuotedLine(line[len("##contig=<") : len(line)-1])
	if err != nil {
		return nil, err
	}
	d := &ContigDef{ID: vals["ID"], Length: -1, OrigLine: line}
	if l, ok := vals["length"]; ok {
		n, err := strconv.ParseInt(l, 10, 64)
		if err != nil {
			return nil, fmt.Errorf("vcf: invalid contig length in %q: %w", line, err)
		}
		d.Length = n
	}
	return d, nil
}

func parseAltDef(line string) (*AltDef, error) {
	if !strings.HasPrefix(line, "##ALT=<") || !strings.HasSuffix(line, ">") {
		return nil, fmt.Errorf("vcf: cannot parse header line: %s", line)
	}
	vals, err := parseQuotedLine(line[len("##ALT=<") : len(line)-1])
	if err != nil {
		return nil, err
	}
	return &AltDef{ID: vals["ID"], Description: vals["Description"], OrigLine: line}, nil
}

// parseQuotedLine parses a comma-separated key=value list where values may be
// double-quoted (and may themselves contain commas or escaped quotes). It ports
// VCFHeader.parseQuotedLine from ngsutilsj.
func parseQuotedLine(s string) (map[string]string, error) {
	values := map[string]string{}
	var key string
	var acc strings.Builder
	haveKey := false
	inquote := false

	for i := 0; i < len(s); i++ {
		c := s[i]
		if !haveKey {
			if c == '=' {
				key = acc.String()
				haveKey = true
				inquote = false
				acc.Reset()
			} else {
				acc.WriteByte(c)
			}
		} else {
			switch {
			case !inquote && c == '"':
				inquote = true
			case inquote && c == '"':
				inquote = false
			case !inquote && c == ',':
				values[key] = acc.String()
				haveKey = false
				acc.Reset()
			case inquote && c == '\\' && i < len(s)-1:
				acc.WriteByte(s[i+1])
				i++
			default:
				acc.WriteByte(c)
			}
		}
	}
	if haveKey {
		values[key] = acc.String()
	}
	return values, nil
}

// quoteString escapes backslashes and double quotes for a ##header value.
func quoteString(s string) string {
	s = strings.ReplaceAll(s, "\\", "\\\\")
	s = strings.ReplaceAll(s, "\"", "\\\"")
	return s
}
