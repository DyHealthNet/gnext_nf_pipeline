#!/usr/bin/env python3
import argparse
from asyncio import as_completed
import os
import re
import sys
import traceback
import json
import math
import gzip
import sys, os  
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed 
import time

# Import zorp with retry logic to handle race conditions
max_retries = 5
for attempt in range(max_retries):
    try:
        from locuszoom import manhattan, qq
        from zorp import sniffers, lookups
        break
    except json.JSONDecodeError as e:
        if attempt < max_retries - 1:
            wait_time = 2 ** attempt  # exponential backoff: 1, 2, 4, 8, 16 seconds
            print(f"Zorp import failed (attempt {attempt + 1}/{max_retries}), retrying in {wait_time}s...", file=sys.stderr)
            time.sleep(wait_time)
        else:
            print(f"Failed to import zorp after {max_retries} attempts", file=sys.stderr)
            raise

def generate_manhattan(reader, out_filename: str, manhattan_num_unbinned=500, manhattan_peak_max_count=500, manhattan_peak_pval_threshold=1e-6, manhattan_peak_sprawl_dist=200_000) -> bool:
    binner = manhattan.Binner(
        num_unbinned=manhattan_num_unbinned,
        peak_max_count=manhattan_peak_max_count,
        peak_neg_log_pval_threshold= -math.log10(manhattan_peak_pval_threshold),
        peak_sprawl_dist=manhattan_peak_sprawl_dist
    )
    print("Processing variants for Manhattan plot...", file=sys.stderr)
    for variant in reader:
        binner.process_variant(variant)
    print("Getting results of the binning...", file=sys.stderr)
    manhattan_data = binner.get_result()

    with open(out_filename, 'w') as f:
        json.dump(manhattan_data, f)
    return True

def generate_qq(reader, out_filename: str) -> bool:
    variants = list(qq.augment_variants(reader))

    rv = {}
    if variants:
        if variants[0].maf is not None:
            rv['overall'] = qq.make_qq_unstratified(variants, include_qq=False)
            rv['by_maf'] = qq.make_qq_stratified(variants)
            rv['ci'] = list(qq.get_confidence_intervals(len(variants) / len(rv['by_maf'])))
        else:
            rv['overall'] = qq.make_qq_unstratified(variants, include_qq=True)
            rv['ci'] = list(qq.get_confidence_intervals(len(variants)))

    with open(out_filename, 'w') as f:
        json.dump(rv, f)
    return True

def process_file(gwas_file, phenocode, manhattan_num_unbinned=500, within_pheno_mask_around_peak=500_000, between_pheno_mask_around_peak=1_000_000, manhattan_peak_max_count=500, manhattan_peak_pval_threshold=1e-6, manhattan_peak_sprawl_dist=200_000, manhattan_peak_variant_counting_pval_threshold=5e-8):
    manhattan_file = f"{phenocode}_manhattan.json"
    qq_file = f"{phenocode}_qq.json"
    
    # Read data for Manhattan plot
    reader_for_manhattan = sniffers.guess_gwas_standard(gwas_file).add_filter('neg_log_pvalue')
    start = time.time()
    generate_manhattan(reader_for_manhattan, manhattan_file, manhattan_num_unbinned, manhattan_peak_max_count, manhattan_peak_pval_threshold, manhattan_peak_sprawl_dist)
    print(f"Time for generating Manhattan plot: {time.time() - start} seconds", file=sys.stderr)
    
    # Read data for QQ plot
    reader_for_qq = sniffers.guess_gwas_standard(gwas_file).add_filter('neg_log_pvalue')
    start = time.time()
    generate_qq(reader_for_qq, qq_file)
    print(f"Time for generating QQ plot: {time.time() - start} seconds", file=sys.stderr)
    return manhattan_file, qq_file

def main():
    parser = argparse.ArgumentParser(description="Generate Manhattan and QQ JSON files from a GWAS .gz file")
    parser.add_argument("--input-files", required=True, help="TSV with two columns: <gwas_file.gz> <phenocode>")
    parser.add_argument("--max-workers", type=int, default=4, help="Maximum number of parallel workers")
    parser.add_argument("--manhattan-num-unbinned", type=int, default=500, help="Number of unbinned variants for Manhattan plot")
    parser.add_argument("--manhattan-peak-max-count", type=int, default=500, help="Max count of peaks for Manhattan plot")
    parser.add_argument("--manhattan-peak-pval-threshold", type=float, default= 1e-6, help="P-value threshold for peaks in Manhattan plot")
    parser.add_argument("--manhattan-peak-sprawl-dist", type=int, default=200_000, help="Sprawl distance for peaks in Manhattan plot")

    args = parser.parse_args()

    # Pre-initialize zorp lookups in main process to avoid race condition in workers
    # This forces the manifest download/update to happen once before forking
    try:
        _ = lookups.get_dbsnp_by_rsid()  # trigger asset initialization
        print("Zorp assets pre-initialized successfully", file=sys.stderr)
    except Exception as e:
        print(f"Warning: Could not pre-initialize zorp assets: {e}", file=sys.stderr)

    input_files_dt = pd.read_csv(args.input_files, sep="\t", header=None, names=["file_path", "phenocode"])
    with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        futures = [
            executor.submit(
                process_file,
                row["file_path"],
                row["phenocode"],
                args.manhattan_num_unbinned,
                args.manhattan_peak_max_count,
                args.manhattan_peak_pval_threshold,
                args.manhattan_peak_sprawl_dist
            )
            for index, row in input_files_dt.iterrows()
        ]
        for future in as_completed(futures):
            try:
                manhattan_file, qq_file = future.result()
                print(f"[OK] Generated {manhattan_file}, {qq_file}", file=sys.stderr)
            except Exception as e:
                print(f"[ERROR] {e}", file=sys.stderr)


if __name__ == "__main__":
    main()