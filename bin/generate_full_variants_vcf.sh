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

# tmp directory
TMP_DIR=$(mktemp -d)
trap "rm -rf $TMP_DIR" EXIT

# Check if GNU parallel is available
if command -v parallel &>/dev/null; then
  echo "Using GNU parallel for efficient processing..."
  
  # Use GNU parallel with streaming to single sort (like the working version)
  printf '%s\n' "${FILES[@]}" \
    | parallel -j"$NUM_JOBS" --tmpdir "$TMP_DIR" --no-notice '
        gzip -cd {} \
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
        "
    ' \
    | sort -k1,1 -k2,2n -u -S2G --compress-program=gzip -T "$TMP_DIR" \
    >> "$OUTPUT_FILE"
  
else
  echo "GNU parallel not available, using fallback method with reduced parallelism..."
  
  # Fallback: reduce NUM_JOBS to avoid OOM
  MAX_JOBS=8
  if [ "$NUM_JOBS" -gt "$MAX_JOBS" ]; then
    echo "WARNING: Reducing NUM_JOBS from $NUM_JOBS to $MAX_JOBS to prevent OOM" >&2
    NUM_JOBS=$MAX_JOBS
  fi
  
  # Extract and stream to single sort (no per-file sorting)
  IDX=0
  PIDS=()
  for FILE in "${FILES[@]}"; do
    (
      gzip -cd "$FILE" \
      | awk -F"\t" -v OFS="\t" '
          NR==1 {
            for (i=1; i<=NF; i++) {
              h=tolower($i)
              if (h=="#chrom" || h=="chrom") c=i
              if (h=="pos")                  p=i
              if (h=="ref")                  r=i
              if (h=="alt")                  a=i
            }
            next
          }
          c && p && r && a {
            print $c, $p, ".", $r, $a, ".", ".", "."
          }
        ' > "$TMP_DIR/extracted_$IDX.txt"
    ) &
    
    PIDS+=($!)
    IDX=$((IDX + 1))
    
    # Limit concurrent extractions
    if (( IDX % NUM_JOBS == 0 )); then
      for PID in "${PIDS[@]}"; do wait "$PID"; done
      PIDS=()
    fi
  done
  
  # Wait for remaining
  for PID in "${PIDS[@]}"; do wait "$PID"; done
  
  # Single sort of all extracted data
  echo "Sorting and deduplicating all variants..."
  cat "$TMP_DIR"/extracted_*.txt \
    | sort -k1,1 -k2,2n -u -S2G --compress-program=gzip -T "$TMP_DIR" \
    >> "$OUTPUT_FILE"
fi

chmod 777 "$OUTPUT_FILE"
echo "Done. VCF written to: $OUTPUT_FILE at $(date)"