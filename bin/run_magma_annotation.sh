#!/bin/bash

# Annotate SNPs to genes with MAGMA (step 1 of the gene-level analysis).
# Maps each SNP to nearby genes using an up/downstream window, producing the
# "<OUT_PREFIX>.genes.annot" file that the gene test later consumes.
#
# Positional args:
#   $1 SNP_LOC     - SNP location file (chr/pos), e.g. the reference .bim
#   $2 GENE_LOC    - gene location file (gene id, chr, start, end)
#   $3 WINDOW_UP   - upstream window size (kb)
#   $4 WINDOW_DOWN - downstream window size (kb)
#   $5 OUT_PREFIX  - output prefix for the .genes.annot file
set -euo pipefail

SNP_LOC=$1

GENE_LOC=$2

WINDOW_UP=$3

WINDOW_DOWN=$4

OUT_PREFIX=$5

magma \
  --annotate window=$WINDOW_UP,$WINDOW_DOWN \
  --snp-loc $SNP_LOC \
  --gene-loc $GENE_LOC \
  --out $OUT_PREFIX