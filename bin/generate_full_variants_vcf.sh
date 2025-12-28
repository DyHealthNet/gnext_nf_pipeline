#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

# Usage: generate_full_variants_vcf.sh OUTPUT_FILE NUM_JOBS BATCH_SIZE FILE1.gz FILE2.gz ...
OUTPUT_FILE="$1"
NUM_JOBS="$2"
shift 3
FILES=("$@")

echo "Starting VCF generation at $(date)"
echo "Output file: $OUTPUT_FILE"
echo "Number of parallel jobs: $NUM_JOBS"
echo "Input files: ${#FILES[@]} files"

echo "Batch size: $BATCH_SIZE"

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
export TMPDIR="$TMP_DIR"


# Step 1: Split FILES into batches of size BATCH_SIZE
echo "Splitting input files into batches of $BATCH_SIZE..."
TOTAL_FILES=${#FILES[@]}
NUM_BATCHES=$(( (TOTAL_FILES + BATCH_SIZE - 1) / BATCH_SIZE ))
BATCH_LIST=()
for ((i=0; i<NUM_BATCHES; i++)); do
  BATCH_FILE="$TMP_DIR/batch_$i.list"
  BATCH_LIST+=("$BATCH_FILE")
  start=$((i * BATCH_SIZE))
  end=$((start + BATCH_SIZE))
  if (( end > TOTAL_FILES )); then end=$TOTAL_FILES; fi
  printf '%s\n' "${FILES[@]:start:end-start}" > "$BATCH_FILE"
done

# Step 2: For each batch, process all files, sort, and deduplicate, writing to a batch VCF
echo "Processing batches in parallel..."
parallel -j"$NUM_JOBS" --tmpdir "$TMP_DIR" --compress --no-notice --line-buffer \
  'BATCH_IDX={#}; BATCH_FILE={}; \
   BATCH_OUT="$TMPDIR/batch_variants_$BATCH_IDX.vcf"; \
   xargs -a "$BATCH_FILE" -I{} bash -c "gzip -cd {} | awk -F\"\t\" -v OFS=\"\t\" '
       NR==1 {
         for (i=1; i<=NF; i++) {
           h=tolower($i)
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
         print $c, $p, \".\", $r, $a, \".\", \".\", \".\"
       }
     '"" \
   | sort -k1,1 -k2,2n \
   | uniq \
   > "$BATCH_OUT" && echo "Completed: $BATCH_FILE -> $BATCH_OUT" >&2' ::: "${BATCH_LIST[@]}"

# Step 3: Merge all batch outputs, sort and deduplicate again
echo "Merging all batch VCFs..."
find "$TMP_DIR" -name 'batch_variants_*.vcf' \
  | sort \
  | xargs cat \
  | sort -k1,1 -k2,2n -S10G --temporary-directory="$TMP_DIR" \
  | uniq \
  >> "$OUTPUT_FILE"

# Step 2: Merge all per-file outputs, sort and deduplicate again
echo "Merging all per-file VCFs..."
find "$TMP_DIR" -name '*.vcf' \
  | sort \
  | xargs cat \
  | sort -k1,1 -k2,2n -S10G --temporary-directory="$TMP_DIR" \
  | uniq \
  >> "$OUTPUT_FILE"

chmod 777 "$OUTPUT_FILE"
echo "Done. VCF written to: $OUTPUT_FILE at $(date)"