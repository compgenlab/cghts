package annotate

import (
	"bufio"
	"fmt"
	"os"
	"strings"

	"github.com/compgenlab/cghts/vcf"
)

// InfoFileOptions configures an [InfoInFile]. The file is a plain list of keys
// (one per line, '#' comments skipped); with Col set it is tab-delimited and the
// first column is the key while column Col supplies a value.
type InfoFileOptions struct {
	Filename  string // lookup file
	Tag       string // record INFO field whose value is looked up
	FlagName  string // INFO key to add
	Delimiter string // if the record value may hold multiple entries (e.g. ","), the separator
	Col       int    // 1-based column supplying the added value; 0 = presence flag
}

// InfoInFile flags (or annotates) a record when the value of one of its INFO
// fields is present in a plain text lookup file loaded into memory. It ports
// ngsutilsj InfoInFile (--in-file).
type InfoInFile struct {
	closeNoop
	opts     InfoFileOptions
	values   map[string]string // key -> value (value used only when Col > 0)
	hasValue bool
}

// NewInfoInFile loads the lookup file and returns the annotator.
func NewInfoInFile(opts InfoFileOptions) (*InfoInFile, error) {
	f, err := os.Open(opts.Filename)
	if err != nil {
		return nil, fmt.Errorf("annotate: open %s: %w", opts.Filename, err)
	}
	defer f.Close()

	values := map[string]string{}
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 1024*1024), 10*1024*1024)
	for sc.Scan() {
		line := sc.Text()
		if len(line) > 0 && line[0] == '#' {
			continue
		}
		if opts.Col > 0 {
			cols := strings.Split(strings.TrimSpace(line), "\t")
			key := cols[0]
			val := ""
			if opts.Col-1 < len(cols) {
				val = cols[opts.Col-1]
			}
			values[key] = val
		} else {
			key := strings.TrimSpace(line)
			if key == "" {
				continue
			}
			values[key] = ""
		}
	}
	if err := sc.Err(); err != nil {
		return nil, fmt.Errorf("annotate: reading %s: %w", opts.Filename, err)
	}

	return &InfoInFile{opts: opts, values: values, hasValue: opts.Col > 0}, nil
}

// SetupHeader declares the added INFO field.
func (a *InfoInFile) SetupHeader(h *vcf.VcfHeader) error {
	h.AddInfo(infoDefSrc(a.opts.FlagName, "0", "Flag", "Is value "+a.opts.Tag+" present in file", a.opts.Filename))
	return nil
}

// Annotate looks up the record's INFO[Tag] value(s) in the file.
func (a *InfoInFile) Annotate(rec *vcf.VcfRecord) error {
	v, ok := rec.InfoValue(a.opts.Tag)
	if !ok {
		return nil
	}
	val := v.String()
	if val == "." {
		return nil
	}
	if a.opts.Delimiter == "" {
		a.apply(rec, val)
	} else {
		for _, piece := range strings.Split(val, a.opts.Delimiter) {
			a.apply(rec, piece)
		}
	}
	return nil
}

func (a *InfoInFile) apply(rec *vcf.VcfRecord, key string) {
	mapped, found := a.values[strings.TrimSpace(key)]
	if !found {
		return
	}
	if !a.hasValue {
		rec.AddInfoFlag(a.opts.FlagName)
	} else if mapped != "" {
		rec.AddInfo(a.opts.FlagName, mapped)
	}
}
