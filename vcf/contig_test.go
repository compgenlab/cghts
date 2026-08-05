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
		// Case. The "chr" prefix used to be trimmed before the upper-casing, and
		// TrimPrefix is case-sensitive, so every spelling but the lowercase one
		// fell through to the accession parser and resolved to nothing.
		{"Chr1", "1", true},
		{"CHR1", "1", true},
		{"CHRX", "X", true},
		{"chrx", "X", true},
		{"x", "X", true},
		{"mt", "MT", true},
		{"ChrM", "MT", true},
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

// IsPrimaryHuman and CanonicalContig share one normalization, so they cannot
// disagree about what "primary" means for a given spelling. They used to: the
// former did not normalize the mitochondrion at all, and neither handled case.
func TestIsPrimaryHumanAgreesWithCanonicalContig(t *testing.T) {
	for _, name := range []string{
		"1", "22", "X", "Y", "M", "MT",
		"chr1", "Chr1", "CHR1", "chrX", "chrx", "chrM", "chrMT", "ChrM",
		"chr1_KI270706v1_random", "GL000009.2", "",
	} {
		_, canon := CanonicalContig(name)
		primary := IsPrimaryHuman(name)
		// CanonicalContig additionally resolves NCBI accessions, which are not
		// "chr"-named, so agreement is asserted only over the naming forms
		// IsPrimaryHuman claims to handle.
		if primary != canon {
			t.Errorf("%q: IsPrimaryHuman = %v but CanonicalContig ok = %v", name, primary, canon)
		}
	}
}
