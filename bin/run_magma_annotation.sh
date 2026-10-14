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

SNP_PREFIX=$1
GENE_LOC=$2
WINDOW_UP=$3
WINDOW_DOWN=$4
OUT_PREFIX=$5

TMP_BIM="${OUT_PREFIX}.merged.bim"
FINAL_BIM=""

echo "[INFO] SNP prefix: $SNP_PREFIX"

# ---------------------------------------
# Detect single vs chromosome-split BIM
# ---------------------------------------

if [[ -f "${SNP_PREFIX}.bim" ]]; then
  echo "[INFO] Detected single BIM file"
  FINAL_BIM="${SNP_PREFIX}.bim"

else
  echo "[INFO] Searching for chromosome-split BIM files"

  rm -f "$TMP_BIM"
  bim_files=()

  shopt -s nullglob
  for f in "${SNP_PREFIX}"*.bim; do
    stem="${f%.bim}"
    suffix="${stem#${SNP_PREFIX}}"

    if [[ "$suffix" =~ ^([1-9]|1[0-9]|2[0-5]|X|Y|M|MT)$ ]]; then
      bim_files+=( "$f" )
    fi
  done
  shopt -u nullglob

  if [[ ${#bim_files[@]} -eq 0 ]]; then
    echo "ERROR: could not detect BIM input from prefix: $SNP_PREFIX" >&2
    echo "Expected either:" >&2
    echo "  ${SNP_PREFIX}.bim" >&2
    echo "or one or more files matching:" >&2
    echo "  ${SNP_PREFIX}{1..25}.bim" >&2
    echo "  ${SNP_PREFIX}X.bim" >&2
    echo "  ${SNP_PREFIX}Y.bim" >&2
    echo "  ${SNP_PREFIX}M.bim" >&2
    echo "  ${SNP_PREFIX}MT.bim" >&2
    exit 1
  fi

  mapfile -t bim_files < <(printf '%s\n' "${bim_files[@]}" | sort -V)

  for chr_file in "${bim_files[@]}"; do
    echo "[INFO] Adding $chr_file"
    cat "$chr_file" >> "$TMP_BIM"
  done

  FINAL_BIM="$TMP_BIM"
fi

# ---------------------------------------
# Run MAGMA annotation
# ---------------------------------------

magma \
  --annotate window=${WINDOW_UP},${WINDOW_DOWN} \
  --snp-loc "$FINAL_BIM" \
  --gene-loc "$GENE_LOC" \
  --out "$OUT_PREFIX"

# ---------------------------------------
# Cleanup
# ---------------------------------------

if [[ "$FINAL_BIM" == "$TMP_BIM" ]]; then
  rm -f "$TMP_BIM"
fi