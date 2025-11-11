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
echo "Input files: ${FILES[*]}"

# header
{
  printf '##fileformat=VCFv4.2\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
} > "$OUTPUT_FILE"

# tmp directory
TMP_DIR=$(mktemp -d)
trap "rm -rf $TMP_DIR" EXIT

# Iterate over all files one by one
for FILE in "${FILES[@]}"; do
  PHENO=$(basename "$FILE" .gz)
  gzip -cd "$FILE" \
  | awk -v pheno="$PHENO" -F"\t" -v OFS="\t" '
      NR==1 {
        for (i=1; i<=NF; i++) {
          h=tolower($i)
          if (h=="#chrom" || h=="chrom") c=i
          if (h=="pos")                  p=i
          if (h=="ref")                  r=i
          if (h=="alt")                  a=i
        }
        if (!c || !p || !r || !a) {
          print "ERROR: header missing columns in " pheno > "/dev/stderr"
          exit 1
        }
        next
      }
      c && p && r && a {
        print $c, $p, ".", $r, $a, ".", ".", "."
      }
    '
done \
| sort -k1,1 -k2,2n -S10G --temporary-directory="$TMP_DIR" \
| uniq \
>> "$OUTPUT_FILE"

chmod 777 "$OUTPUT_FILE"
echo "Done. VCF written to: $OUTPUT_FILE"