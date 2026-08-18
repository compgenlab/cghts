#!/bin/sh
# Summarise the reference-block structure of a gVCF.
#   usage: gvcf-blocks.sh sample.g.vcf.gz [chrom]
# Reads the whole file (or one chromosome); prints a few lines. No dependencies
# beyond zcat and awk.
f="$1"; chrom="$2"
echo "== banding declared in the header =="
zcat "$f" | sed -n '/^#CHROM/q;/^##GVCFBlock/p' | head -20
echo
zcat "$f" | awk -F'\t' -v want="$chrom" '
  /^#/ { next }
  want != "" && $1 != want { next }
  {
    end = ""
    if (match($8, /(^|;)END=[0-9]+/)) { end = substr($8, RSTART, RLENGTH); sub(/.*END=/, "", end) }
    if (end == "") { vars++; next }
    blocks++
    len = end - $2 + 1
    total += len
    if (len > max) max = len
    if      (len <    10) h1++
    else if (len <   100) h2++
    else if (len <  1000) h3++
    else if (len < 10000) h4++
    else                  h5++
  }
  END {
    printf "ref_blocks        %d\n", blocks
    printf "variant_records   %d\n", vars
    printf "bases_in_blocks   %d\n", total
    printf "mean_block_len    %.1f\n", (blocks ? total/blocks : 0)
    printf "max_block_len     %d\n", max
    printf "block_len  <10    %d\n", h1
    printf "block_len  <100   %d\n", h2
    printf "block_len  <1k    %d\n", h3
    printf "block_len  <10k   %d\n", h4
    printf "block_len  >=10k  %d\n", h5
  }'
