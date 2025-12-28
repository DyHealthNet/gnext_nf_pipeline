#!/usr/bin/env bash
set -euo pipefail

# merge_vcfs.sh: Merge, sort, and deduplicate multiple VCF files
# Usage: merge_vcfs.sh OUTPUT_VCF BATCH_VCF1 [BATCH_VCF2 ...]
#   OUTPUT_VCF: final merged VCF file
#   BATCH_VCF*: input VCF files to merge

OUTPUT_VCF="$1"
shift
BATCH_VCFS=("$@")

# header (from first file)
head -n 2 "${BATCH_VCFS[0]}" > "$OUTPUT_VCF"

TMP_DIR=$(dirname "$OUTPUT_VCF")/tmp_merge_$$
mkdir -p "$TMP_DIR"
trap "rm -rf $TMP_DIR" EXIT

cat "${BATCH_VCFS[@]}" | grep -v '^#' \
  | sort -k1,1 -k2,2n -S10G --temporary-directory="$TMP_DIR" \
  | uniq \
  >> "$OUTPUT_VCF"

echo "VCF merge complete: $OUTPUT_VCF"

# Compress and index the merged VCF
echo "Compressing and indexing final VCF with bgzip/tabix..."
bgzip -f "$OUTPUT_VCF"
tabix -p vcf "$OUTPUT_VCF.gz"
echo "VCF compressed and indexed: $OUTPUT_VCF.gz (.tbi)"