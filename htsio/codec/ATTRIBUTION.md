# Provenance and attribution

The entropy codecs in this package are clean-room Go implementations of the
CRAM-related codecs from the HTS ecosystem. They were written from the published
formats and reference sources below; no upstream source code was copied.

## Sources

- **rANS (4x8 and Nx16), arithmetic/range coding, fqzcomp quality model,
  name tokenizer, simple adaptive model** — based on the
  [htscodecs](https://github.com/samtools/htscodecs) reference implementations
  (`rANS_static.c`, `rANS_static32x16pr.c`, `arith_dynamic.c`,
  `fqzcomp_qual.c`, `tokenise_name3.c`, `c_simple_model.h`).
  htscodecs is distributed under the **MIT/BSD license**.

- **Range coder** — based on Eugene Shelwien's carryless range coder, which is
  released into the **public domain** and is the same coder used by htscodecs.

The relevant upstream license terms (MIT/BSD for htscodecs; public domain for
the Shelwien range coder) permit reimplementation and redistribution. Individual
files carry `// Based on ...` comments pointing at the specific reference they
follow.

## Compatibility intent

These implementations aim for bit-exact interoperability with `samtools`/`htslib`
CRAM output. samtools is treated as the reference: any input it accepts must
decode here, and our encoded output is validated against it. The fuzz targets in
`fuzz_test.go` exercise both the decode-no-panic and encode/decode round-trip
properties.
