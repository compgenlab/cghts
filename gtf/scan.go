package gtf

import (
	"bufio"
	"io"
	"strings"
)

// GeneRef is a gene's identity without its structure: the pair a membership
// question is asked in terms of.
type GeneRef struct {
	GeneID   string
	GeneName string
}

// ScanGenes streams the distinct genes out of a GTF, calling visit once per
// gene_id in the order each is first seen. Stops and returns visit's error.
//
// Separate from NewAnnotationSource because the two want different things from
// the same file. Building the model means holding every transcript, exon and CDS
// of every gene in memory to answer positional queries; a caller that only wants
// to know which genes exist pays gigabytes for an answer that is a few megabytes
// of strings. GENCODE v48 is ~1.7M rows and ~78k genes.
//
// Only the identity columns are read, so this does not apply the requiredTags
// filter or any of the biotype derivation — a gene present in the file is
// reported, whatever a downstream annotator would do with its features.
//
// Deduplication is by gene_id, not by (ref, gene_id) as the model uses: a gene
// listed on two contigs is one gene to anybody asking whether a symbol exists,
// and reporting it twice would only make callers dedupe again.
func ScanGenes(r io.Reader, visit func(GeneRef) error) error {
	seen := make(map[string]bool)

	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
	for sc.Scan() {
		line := sc.Text()
		if len(line) == 0 || line[0] == '#' {
			continue
		}
		cols := strings.Split(line, "\t")
		if len(cols) < 9 {
			continue
		}
		a := parseAttributes(cols[8])
		if a.geneID == "" {
			continue
		}
		if seen[a.geneID] {
			continue
		}
		seen[a.geneID] = true

		name := a.geneName
		if name == "" {
			name = a.gene // RefSeq GTFs use "gene" rather than "gene_name"
		}
		if err := visit(GeneRef{GeneID: a.geneID, GeneName: name}); err != nil {
			return err
		}
	}
	return sc.Err()
}

// TrimGeneIDVersion removes the version suffix an Ensembl or GENCODE gene_id
// carries: ENSG00000141510.17 → ENSG00000141510.
//
// The suffix counts the revisions of a gene's model, not the gene, so two
// releases name the same gene differently and a list written against one stops
// matching the other. Nothing that asks "is this gene in my set" wants that.
//
// The rule is a literal ".", one or more digits, and then either end-of-string
// or "_" — the second case being GENCODE's pseudoautosomal ids, where the suffix
// is in the middle: ENSG00000182378.14_PAR_Y → ENSG00000182378_PAR_Y.
//
// Deliberately not restricted to ids that look Ensembl. An id this over-trims is
// over-trimmed identically wherever the rule is applied, so both sides of a
// comparison still agree — which is the property that matters, and one a
// prefix test would give up in exchange for nothing.
func TrimGeneIDVersion(id string) string {
	for i := 0; i < len(id); i++ {
		if id[i] != '.' {
			continue
		}
		j := i + 1
		for j < len(id) && id[j] >= '0' && id[j] <= '9' {
			j++
		}
		if j == i+1 {
			continue // "." with no digits after it
		}
		if j == len(id) || id[j] == '_' {
			return id[:i] + id[j:]
		}
	}
	return id
}
