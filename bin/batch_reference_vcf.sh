#!/usr/bin/env bash

# Build a minimal, sorted, de-duplicated VCF of variant sites from a batch of
# normalized GWAS files. Each input contributes its CHROM/POS/REF/ALT columns;
# the union of sites becomes one VCF that is later merged into the reference VCF
# used for VEP annotation. Extraction is parallelized across files.
#
# Positional args:
#   $1 BATCH_LIST - file listing the normalized .gz inputs (one path per line)
#   $2 OUTPUT_VCF - destination VCF path
#   $3 NUM_JOBS   - parallel jobs / sort threads (default 16)
set -euo pipefail

BATCH_LIST="$1"
OUTPUT_VCF="$2"
NUM_JOBS="${3:-16}"

# Write a minimal VCF header (sites-only, no sample/genotype columns).
{
  printf '##fileformat=VCFv4.2\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
} > "$OUTPUT_VCF"

# Per-process temp dir for sort spill files; cleaned up on exit.
TMP_DIR=$(dirname "$OUTPUT_VCF")/tmp_extract_$$
mkdir -p "$TMP_DIR"
trap "rm -rf $TMP_DIR" EXIT

echo "Starting extraction of $(wc -l < "$BATCH_LIST") files with $NUM_JOBS jobs..." >&2

# Convert one normalized file into bare VCF site rows: drop the header line and
# remap columns to CHROM POS . REF ALT . . . (ID/QUAL/FILTER/INFO left empty).
process_file() {
  local file="$1"
  gzip -cd "$file" \
    | tail -n +2 \
    | awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, ".", $4, $5, ".", ".", "."}' \
    && echo "Done: $(basename "$file")" >&2
}
export -f process_file

# Extract all files in parallel, then sort by chrom/pos and collapse duplicate
# sites before appending to the VCF body.
cat "$BATCH_LIST" \
  | parallel -j"$NUM_JOBS" --no-notice --line-buffer process_file {} \
  | sort -k1,1 -k2,2n -S10G --parallel="$NUM_JOBS" --compress-program=gzip -T "$TMP_DIR" \
  | uniq \
  >> "$OUTPUT_VCF"

echo "Batch extraction complete: $OUTPUT_VCF" >&2