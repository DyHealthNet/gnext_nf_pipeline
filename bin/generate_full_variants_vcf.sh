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
# Don't trap cleanup yet - wait until all subshells complete

# Process files in parallel using a simpler approach
echo "Processing files in parallel (up to $NUM_JOBS at a time)..."
echo "Starting loop over ${#FILES[@]} files..."
IDX=0
PIDS=()
for FILE in "${FILES[@]}"; do
  echo "Launching job $IDX for: $(basename "$FILE")"
  (
    LOCAL_TMP_DIR="$TMP_DIR"
    LOCAL_IDX="$IDX"
    PHENO=$(basename "$FILE" .gz)
    OUTPUT_SORTED="${LOCAL_TMP_DIR}/sorted_${LOCAL_IDX}.txt"
    
    # Test if we can write to temp dir
    if ! touch "$OUTPUT_SORTED" 2>/dev/null; then
      echo "ERROR: Cannot write to $OUTPUT_SORTED" >&2
      exit 1
    fi
    rm -f "$OUTPUT_SORTED"
    
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
          if (!c || !p || !r || !a) {
            print "ERROR: header missing columns" > "/dev/stderr"
            exit 1
          }
          next
        }
        c && p && r && a {
          print $c, $p, ".", $r, $a, ".", ".", "."
        }
      ' \
    | sort -k1,1 -k2,2n -u -S1G -T "$LOCAL_TMP_DIR" > "$OUTPUT_SORTED"
    
    SORT_EXIT=$?
    if [ $SORT_EXIT -ne 0 ]; then
      echo "ERROR: Sort failed for $PHENO (exit: $SORT_EXIT)" >&2
      exit 1
    fi
    
    if [ ! -s "$OUTPUT_SORTED" ]; then
      echo "ERROR: Empty output for $PHENO" >&2
      exit 1
    fi
    
    echo "Completed: $PHENO ($(wc -l < "$OUTPUT_SORTED") variants)" >&2
  ) &
  
  PIDS+=($!)
  IDX=$((IDX + 1))
  
  # Limit concurrent jobs
  if (( IDX % NUM_JOBS == 0 )); then
    echo "Waiting for batch to complete..." >&2
    FAILED=0
    for PID in "${PIDS[@]}"; do
      if ! wait "$PID"; then
        WAIT_EXIT=$?
        echo "ERROR: Job PID $PID failed with exit code $WAIT_EXIT" >&2
        FAILED=1
      fi
    done
    if [ $FAILED -eq 1 ]; then
      echo "ERROR: One or more jobs in batch failed" >&2
      exit 1
    fi
    PIDS=()
  fi
done

# Wait for remaining jobs
echo "Waiting for final batch..." >&2
FAILED=0
for PID in "${PIDS[@]}"; do
  if ! wait "$PID"; then
    WAIT_EXIT=$?
    echo "ERROR: Job PID $PID failed with exit code $WAIT_EXIT" >&2
    FAILED=1
  fi
done
if [ $FAILED -eq 1 ]; then
  echo "ERROR: One or more final jobs failed" >&2
  exit 1
fi

# Now set trap to clean up temp dir after all jobs complete
trap "rm -rf $TMP_DIR" EXIT

# Check if any sorted files were created
SORTED_COUNT=$(ls -1 "$TMP_DIR"/sorted_*.txt 2>/dev/null | wc -l)
echo "Found $SORTED_COUNT pre-sorted files"

if [ "$SORTED_COUNT" -eq 0 ]; then
  echo "ERROR: No sorted files were created" >&2
  exit 1
fi

if [ "$SORTED_COUNT" -ne "${#FILES[@]}" ]; then
  echo "WARNING: Expected ${#FILES[@]} files but only got $SORTED_COUNT" >&2
fi

# Merge all sorted files
echo "Merging $SORTED_COUNT files..."
sort -m -k1,1 -k2,2n -u "$TMP_DIR"/sorted_*.txt >> "$OUTPUT_FILE"

chmod 777 "$OUTPUT_FILE"
echo "Done. VCF written to: $OUTPUT_FILE at $(date)"