package vcf

import "testing"

func TestCanonicalContig(t *testing.T) {
	cases := []struct {
		in   string
		want string
		ok   bool
	}{
		// Ensembl
		{"1", "1", true},
		{"22", "22", true},
		{"X", "X", true},
		{"Y", "Y", true},
		{"MT", "MT", true},
		{"M", "MT", true},
		// UCSC
		{"chr1", "1", true},
		{"chrX", "X", true},
		{"chrM", "MT", true},
		{"chrMT", "MT", true},
		// NCBI RefSeq (version ignored; assembly-independent)
		{"NC_000001.11", "1", true},
		{"NC_000001.10", "1", true},
		{"NC_000022.11", "22", true},
		{"NC_000023.11", "X", true},
		{"NC_000024.10", "Y", true},
		{"NC_012920.1", "MT", true},
		// Non-resolving: scaffolds, alt loci, non-human/other accessions
		{"chr1_KI270706v1_random", "", false},
		{"GL000009.2", "", false},
		{"NT_167214.1", "", false},
		{"NC_045512.2", "", false}, // SARS-CoV-2: NC_ but outside the human primary set
		{"NC_000025.1", "", false}, // no chromosome 25
		{"", "", false},
	}
	for _, c := range cases {
		got, ok := CanonicalContig(c.in)
		if got != c.want || ok != c.ok {
			t.Errorf("CanonicalContig(%q) = (%q, %v), want (%q, %v)", c.in, got, ok, c.want, c.ok)
		}
	}
}

func TestContigConverter(t *testing.T) {
	// Target uses UCSC naming.
	conv := NewContigConverter([]string{"chr1", "chr2", "chrX", "chrM", "weird_contig"})

	cases := []struct {
		in   string
		want string
		ok   bool
	}{
		{"1", "chr1", true},                    // Ensembl -> UCSC
		{"chr1", "chr1", true},                 // exact-match fast path
		{"NC_000001.11", "chr1", true},         // NCBI -> UCSC
		{"MT", "chrM", true},                   // Ensembl mito -> UCSC
		{"chrM", "chrM", true},                 // exact
		{"weird_contig", "weird_contig", true}, // exact match for an unresolvable name
		{"X", "chrX", true},
		{"3", "", false},    // no chr3 in target
		{"chrY", "", false}, // no Y in target
	}
	for _, c := range cases {
		got, ok := conv.Resolve(c.in)
		if got != c.want || ok != c.ok {
			t.Errorf("Resolve(%q) = (%q, %v), want (%q, %v)", c.in, got, ok, c.want, c.ok)
		}
	}
}

func TestContigConverterReverse(t *testing.T) {
	// Target uses Ensembl naming; a UCSC/NCBI query must still resolve.
	conv := NewContigConverter([]string{"1", "2", "MT"})
	for in, want := range map[string]string{
		"chr1":         "1",
		"NC_000001.11": "1",
		"chrM":         "MT",
	} {
		got, ok := conv.Resolve(in)
		if !ok || got != want {
			t.Errorf("Resolve(%q) = (%q, %v), want (%q, true)", in, got, ok, want)
		}
	}
}
