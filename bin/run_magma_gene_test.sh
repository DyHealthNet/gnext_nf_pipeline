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
echo "Reference file/prefix: $REF_FILE"
echo "Annotation file: $ANNO_LOC"
echo "SNP p-value file: $SNP_PVAL_FILE"
echo "Max workers: $MAX_WORKERS"
echo "Sample size parameter: $N_PARAM"
echo "Output prefix: $OUT_PREFIX"

MODE=""
REF_FILE_CHR=""
declare -a CHR_LABELS=()

# ---------------------------------------
# Detect single vs chromosome-split PLINK input
# ---------------------------------------

if [[ -f "${REF_FILE}.bed" && -f "${REF_FILE}.bim" && -f "${REF_FILE}.fam" ]]; then
  MODE="single"

else
  echo "[INFO] Searching for chromosome-split PLINK files"

  shopt -s nullglob
  for bim in "${REF_FILE}"*.bim; do
    prefix="${bim%.bim}"
    suffix="${prefix#$REF_FILE}"

    if [[ "$suffix" =~ ^([1-9]|1[0-9]|2[0-5]|X|Y|M|MT)$ ]]; then
      if [[ -f "${prefix}.bed" && -f "${prefix}.fam" ]]; then
        CHR_LABELS+=( "$suffix" )
      fi
    fi
  done
  shopt -u nullglob

  if [[ ${#CHR_LABELS[@]} -eq 0 ]]; then
    echo "ERROR: could not detect valid PLINK input from prefix: $REF_FILE" >&2
    echo "Expected either:" >&2
    echo "  single dataset: ${REF_FILE}.bed/.bim/.fam" >&2
    echo "or one or more complete triplets matching:" >&2
    echo "  ${REF_FILE}{1..25}.bed/.bim/.fam" >&2
    echo "  ${REF_FILE}X.bed/.bim/.fam" >&2
    echo "  ${REF_FILE}Y.bed/.bim/.fam" >&2
    echo "  ${REF_FILE}M.bed/.bim/.fam" >&2
    echo "  ${REF_FILE}MT.bed/.bim/.fam" >&2
    exit 1
  fi

  MODE="chr"
  REF_FILE_CHR="${REF_FILE}#CHR#"

  mapfile -t CHR_LABELS < <(printf '%s\n' "${CHR_LABELS[@]}" | sort -V)
fi

echo "[INFO] Detected mode: $MODE"

if [[ "$MODE" == "chr" ]]; then
  echo "[INFO] Chromosome template: $REF_FILE_CHR"
  echo "[INFO] Detected chromosome labels: ${CHR_LABELS[*]}"

  for chr in "${CHR_LABELS[@]}"; do
    chr_prefix="${REF_FILE_CHR//#CHR#/$chr}"
    if [[ ! -f "${chr_prefix}.bed" || ! -f "${chr_prefix}.bim" || ! -f "${chr_prefix}.fam" ]]; then
      echo "ERROR: missing PLINK files for '${chr}': ${chr_prefix}.bed/.bim/.fam" >&2
      exit 1
    fi
  done
fi

# ---------------------------------------
# Run MAGMA
# ---------------------------------------

if [[ "$MODE" == "single" ]]; then
  echo "[INFO] Running single-dataset batch mode..."

  seq 1 "$MAX_WORKERS" | parallel -j "$MAX_WORKERS" \
    magma \
      --batch {} "$MAX_WORKERS" \
      --bfile "$REF_FILE" synonyms=0 \
      --gene-annot "$ANNO_LOC" \
      --gene-model snp-wise=mean \
      --pval "$SNP_PVAL_FILE" "$N_PARAM" \
      --out "temp_annot/${OUT_PREFIX}" \
      --genes-only

else
  echo "[INFO] Running chromosome batch mode..."

  printf '%s\n' "${CHR_LABELS[@]}" | parallel -j "$MAX_WORKERS" \
    magma \
      --batch {} chr \
      --bfile "$REF_FILE_CHR" synonyms=0 \
      --gene-annot "$ANNO_LOC" \
      --gene-model snp-wise=mean \
      --pval "$SNP_PVAL_FILE" "$N_PARAM" \
      --out "temp_annot/${OUT_PREFIX}" \
      --genes-only
fi

echo "[INFO] Merging MAGMA batches..."

magma \
  --merge "temp_annot/${OUT_PREFIX}" \
  --out "temp_annot/${OUT_PREFIX}"

echo "[INFO] Copying final .genes.* files..."
cp "temp_annot/${OUT_PREFIX}".genes.* .

rm -r temp_annot

echo "[INFO] Done."
echo "[INFO] Final results copied to current directory:"
echo "       ${OUT_PREFIX}.genes.*"