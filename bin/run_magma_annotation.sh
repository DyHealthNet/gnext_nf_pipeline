#!/bin/bash
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