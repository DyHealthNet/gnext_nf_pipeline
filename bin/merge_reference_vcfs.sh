#!/usr/bin/env bash
set -euo pipefail

# merge_vcfs.sh: Merge, sort, and deduplicate multiple VCF files
# Usage: merge_vcfs.sh MANIFEST_FILE OUTPUT_VCF [CPUS]
#   MANIFEST_FILE: text file with one VCF path per line
#   OUTPUT_VCF: final merged VCF file
#   CPUS: number of parallel threads (optional, default 8)

MANIFEST="$1"
OUTPUT_VCF="$2"
CPUS="${3:-8}"

if [ ! -f "$MANIFEST" ]; then
    echo "Error: Manifest file not found: $MANIFEST" >&2
    exit 1
fi

# Read first VCF from manifest for header
FIRST_VCF=$(head -n 1 "$MANIFEST")
if [ ! -f "$FIRST_VCF" ]; then
    echo "Error: First VCF file not found: $FIRST_VCF" >&2
    exit 1
fi

# Write header from first file
head -n 2 "$FIRST_VCF" > "$OUTPUT_VCF"

TMP_DIR=$(dirname "$OUTPUT_VCF")/tmp_merge_$$
mkdir -p "$TMP_DIR"
trap "rm -rf $TMP_DIR" EXIT

NUM_FILES=$(wc -l < "$MANIFEST")
echo "Merging $NUM_FILES VCF files with $CPUS threads..." >&2

# Extract variants from all files (skip headers), sort and deduplicate
while IFS= read -r vcf; do
    tail -n +3 "$vcf"  # Skip first 2 header lines
done < "$MANIFEST" \
  | sort -k1,1 -k2,2n -S10G --parallel="$CPUS" --compress-program=gzip --temporary-directory="$TMP_DIR" \
  | uniq \
  >> "$OUTPUT_VCF"

echo "VCF merge complete: $OUTPUT_VCF" >&2

# Compress and index the merged VCF
echo "Compressing and indexing final VCF with bgzip/tabix..." >&2
bgzip -f "$OUTPUT_VCF"
tabix -p vcf "${OUTPUT_VCF}.gz"
echo "VCF compressed and indexed: ${OUTPUT_VCF}.gz" >&2