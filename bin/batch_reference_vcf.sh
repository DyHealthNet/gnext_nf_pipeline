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

# Define a function to process each file
process_file() {
  local file="$1"
  gzip -cd "$file" \
    | tail -n +2 \
    | awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, ".", $4, $5, ".", ".", "."}' \
    && echo "Done: $(basename "$file")" >&2
}
export -f process_file

# Extract in parallel with optimized processing
cat "$BATCH_LIST" \
  | parallel -j"$NUM_JOBS" --no-notice --line-buffer process_file {} \
  | sort -k1,1 -k2,2n -S10G --parallel="$NUM_JOBS" --compress-program=gzip -T "$TMP_DIR" \
  | uniq \
  >> "$OUTPUT_VCF"

echo "Batch extraction complete: $OUTPUT_VCF" >&2