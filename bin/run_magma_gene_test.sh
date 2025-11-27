#!/bin/bash
set -euo pipefail

mkdir -p temp_annot

REF_FILE=$1

ANNO_LOC=$2

SNP_PVAL_FILE=$3

MAX_WORKERS=$4

N_PARAM=$5

OUT_PREFIX=$6

echo "MAGMA parameters:"
echo "Reference file: $REF_FILE"
echo "Annotation file: $ANNO_LOC"
echo "SNP p-value file: $SNP_PVAL_FILE"
echo "Max workers: $MAX_WORKERS"
echo "Sample size parameter: $N_PARAM"
echo "Output prefix: $OUT_PREFIX"



# run magma
# N_PARAM can be either "N=12345" or "ncol=12"
magma \
    --bfile "$REF_FILE" synonyms=0 \
    --gene-annot "$ANNO_LOC" \
    --gene-model snp-wise=mean \
    --pval "$SNP_PVAL_FILE" $N_PARAM \
    --out "${OUT_PREFIX}" \
    --genes-only

echo "MAGMA gene analysis completed."


# run magma in parallel

#seq 1 "$MAX_WORKERS" | parallel -j "$MAX_WORKERS" \
#  magma \
#    --batch {} "$MAX_WORKERS" \
#    --bfile "$REF_FILE" synonyms=0 \
#    --gene-annot "$ANNO_LOC" \
#    --gene-model snp-wise=mean \
#    --pval "$SNP_PVAL_FILE" N="$NR_SAMPLES" \
#    --out "temp_annot/${OUT_PREFIX}" \
#    --genes-only

# merge intermediate results generated under the temp_annot files and merge to one single file set

#magma \
#  --merge temp_annot/$OUT_PREFIX \
#  --out temp_annot/$OUT_PREFIX

# extract merged files for subsequent analysis
#cp temp_annot/$OUT_PREFIX.genes.* .

# remove the temporary folder
#rm -r temp_annot