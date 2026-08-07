package seqio

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/compgenlab/cghts/htsio/bgzf"
)

// BuildFastaIndex writes the indexes a random-access reader needs: <path>.fai,
// and <path>.gzi when path is BGZF. It returns the entries written.
//
// The counterpart to NewIndexedFastaReader and bgzf.LoadGZIndex, which could
// read both formats but not produce them — so anything provisioning a reference
// had to shell out to samtools. Output is byte-identical to htslib's.
//
// Plain gzip is refused rather than indexed. A .fai addresses uncompressed
// offsets and only a container that can seek to one can use it: BGZF can, plain
// gzip cannot, so the index would be unusable rather than merely imperfect.
func BuildFastaIndex(path string) ([]FaiEntry, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	magic := make([]byte, 2)
	if _, err := io.ReadFull(f, magic); err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	if _, err := f.Seek(0, io.SeekStart); err != nil {
		return nil, err
	}

	var entries []FaiEntry
	if magic[0] == 0x1f && magic[1] == 0x8b {
		// Streamed, not buffered: a human genome is ~3 GB uncompressed, and
		// holding that to find block boundaries turns indexing into an OOM kill.
		ir := bgzf.NewIndexingReader(f)
		if entries, err = indexFasta(ir); err != nil {
			return nil, fmt.Errorf("%s: %w", path, err)
		}
		if err := bgzf.WriteGZIndex(path+".gzi", ir.Blocks()); err != nil {
			return nil, err
		}
	} else if entries, err = indexFasta(f); err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}

	if err := writeFai(path+".fai", entries); err != nil {
		return nil, err
	}
	return entries, nil
}

// indexFasta walks a FASTA and records where each sequence starts.
//
// A record's lines must be uniform apart from its last. The format stores one
// line length per sequence and computes offsets arithmetically from it, so a
// ragged record cannot be described — and would yield wrong coordinates rather
// than an error, which is the failure worth refusing.
func indexFasta(r io.Reader) ([]FaiEntry, error) {
	var (
		out    []FaiEntry
		cur    *FaiEntry
		pos    int64
		lineNo int64
		short  bool // a short line was seen; only legal as the record's last
	)
	br := bufio.NewReaderSize(r, 1<<20)
	for {
		line, err := br.ReadBytes('\n')
		// A read error is not end-of-input. Treating the two alike would turn a
		// truncated or corrupt stream into a short index that looks complete —
		// and, when the reader is decompressing, would swallow "this is not
		// BGZF" entirely.
		if err != nil && err != io.EOF {
			return nil, err
		}
		if len(line) == 0 && err != nil {
			break
		}
		lineNo++
		width := int64(len(line))
		trimmed := bytes.TrimRight(line, "\r\n")

		switch {
		case len(trimmed) > 0 && trimmed[0] == '>':
			if cur != nil {
				out = append(out, *cur)
			}
			name := ""
			if fields := strings.Fields(string(trimmed[1:])); len(fields) > 0 {
				name = fields[0]
			}
			cur = &FaiEntry{Name: name, Offset: pos + width}
			short = false
		case cur == nil:
			// Anything before the first header belongs to no sequence.
		default:
			bases := len(trimmed)
			// A short line ends the record. Anything after one is ragged
			// whatever its length, so the check is "did a short line already
			// happen", not "is this line a different length" — 8/3/8 passes the
			// second and fails the first.
			if short {
				return nil, fmt.Errorf("sequence %q has uneven line lengths at line %d; "+
					"a .fai cannot describe that", cur.Name, lineNo)
			}
			if cur.LineBases == 0 {
				cur.LineBases, cur.LineWidth = bases, int(width)
			} else if bases != cur.LineBases {
				if bases > cur.LineBases {
					return nil, fmt.Errorf("sequence %q has a long line at line %d; "+
						"a .fai cannot describe that", cur.Name, lineNo)
				}
				short = true
			}
			cur.Length += bases
		}
		pos += width
		if err != nil {
			break
		}
	}
	if cur != nil {
		out = append(out, *cur)
	}
	return out, nil
}

func writeFai(path string, entries []FaiEntry) error {
	var b strings.Builder
	for _, e := range entries {
		fmt.Fprintf(&b, "%s\t%d\t%d\t%d\t%d\n",
			e.Name, e.Length, e.Offset, e.LineBases, e.LineWidth)
	}
	return os.WriteFile(path, []byte(b.String()), 0o644)
}
