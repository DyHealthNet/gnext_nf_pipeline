#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

OUTPUT_FILE="$1"
NUM_JOBS="$2"
shift 2
FILES=("$@")

echo "Starting VCF generation at $(date)"
echo "Output file: $OUTPUT_FILE"
echo "Number of parallel jobs: $NUM_JOBS"
echo "Input files: ${#FILES[@]} files"

# header
{
  printf '##fileformat=VCFv4.2\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
} > "$OUTPUT_FILE"

# tmp directory - use work directory not system temp
WORK_DIR=$(dirname "$OUTPUT_FILE")
TMP_DIR="$WORK_DIR/tmp_vcf_$$"
mkdir -p "$TMP_DIR"
trap "rm -rf $TMP_DIR" EXIT

# Set TMPDIR environment variable for all child processes
export TMPDIR="$TMP_DIR"


echo "Using GNU parallel for efficient processing..."
echo "Processing ${#FILES[@]} files with $NUM_JOBS parallel jobs..."
echo "Streaming pre-sorted files through merge..."

# Stream extraction through sort -m directly (no intermediate files needed)
printf '%s\n' "${FILES[@]}" \
  | parallel -j"$NUM_JOBS" --tmpdir "$TMP_DIR" --compress --no-notice --line-buffer \
      'gzip -cd {} \
      | awk -F"\t" -v OFS="\t" "
          NR==1 {
            for (i=1; i<=NF; i++) {
              h=tolower(\$i)
              if (h==\"#chrom\"||h==\"chrom\") c=i
              if (h==\"pos\")                  p=i
              if (h==\"ref\")                  r=i
              if (h==\"alt\")                  a=i
            }
            if (!c||!p||!r||!a) {
              print \"ERROR: header missing columns\" > \"/dev/stderr\"
              exit 1
            }
            next
          }
          c && p && r && a {
            print \$c, \$p, \".\", \$r, \$a, \".\", \".\", \".\"
          }
      " && echo "Completed: {/}" >&2
  ' \
  | sort -m -k1,1 -k2,2n -S10G --temporary-directory="$TMP_DIR" \
  | uniq \
  >> "$OUTPUT_FILE"


chmod 777 "$OUTPUT_FILE"
echo "Done. VCF written to: $OUTPUT_FILE at $(date)"