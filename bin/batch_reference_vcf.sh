#!/usr/bin/env bash
set -euo pipefail

BATCH_LIST="$1"
OUTPUT_VCF="$2"
NUM_JOBS="${3:-16}"

{
  printf '##fileformat=VCFv4.2\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
} > "$OUTPUT_VCF"

TMP_DIR=$(dirname "$OUTPUT_VCF")/tmp_extract_$$
mkdir -p "$TMP_DIR"
trap "rm -rf $TMP_DIR" EXIT

echo "Starting extraction of $(wc -l < "$BATCH_LIST") files with $NUM_JOBS jobs..." >&2

# Extract in parallel with optimized processing
cat "$BATCH_LIST" \
  | parallel -j"$NUM_JOBS" --no-notice --line-buffer \
      "gzip -cd {} \
      | tail -n +2 \
      | awk -F'\t' 'BEGIN{OFS=\"\t\"} {print \$1, \$2, \".\", \$3, \$4, \".\", \".\", \".\"}' \
      && echo 'Done: {/}' >&2
  " \
  | sort -k1,1 -k2,2n -S10G --parallel="$NUM_JOBS" --compress-program=gzip -T "$TMP_DIR" \
  | uniq \
  >> "$OUTPUT_VCF"

echo "Batch extraction complete: $OUTPUT_VCF" >&2