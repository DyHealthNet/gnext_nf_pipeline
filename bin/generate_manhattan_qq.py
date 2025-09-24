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

from locuszoom import manhattan, qq
from zorp import sniffers, lookups
import time

def add_rsID_with_lmdb(reader, lmdb_path):
    # Open LMDB once (shared by all threads)
    rsid_finder = lookups.SnpToRsid(lmdb_path, test=False)
    reader.add_lookup('rsid', lambda variant: rsid_finder(variant.chrom, variant.pos, variant.ref, variant.alt))

def generate_manhattan(reader, out_filename: str, manhattan_num_unbinned=500, within_pheno_mask_around_peak=500_000, between_pheno_mask_around_peak=1_000_000, manhattan_peak_max_count=500, manhattan_peak_pval_threshold=1e-6, manhattan_peak_sprawl_dist=200_000, manhattan_peak_variant_counting_pval_threshold=5e-8) -> bool:
    binner = manhattan.Binner(
        manhattan_num_unbinned=manhattan_num_unbinned,
        within_pheno_mask_around_peak=within_pheno_mask_around_peak,
        between_pheno_mask_around_peak=between_pheno_mask_around_peak,
        manhattan_peak_max_count=manhattan_peak_max_count,
        manhattan_peak_pval_threshold=manhattan_peak_pval_threshold,
        manhattan_peak_sprawl_dist=manhattan_peak_sprawl_dist,
        manhattan_peak_variant_counting_pval_threshold=manhattan_peak_variant_counting_pval_threshold
    )
    print("Processing variants for Manhattan plot...", file=sys.stderr)
    for variant in reader:
        binner.process_variant(variant)
        
    print("Getting results of the binning...", file=sys.stderr)
    manhattan_data = binner.get_result()
    #for v_dict in manhattan_data['unbinned_variants']:
    #    if math.isinf(v_dict['neg_log_pvalue']):
    #        v_dict['neg_log_pvalue'] = 'Infinity'

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

def process_file(gwas_file, phenocode, lmdb_path, manhattan_num_unbinned=500, within_pheno_mask_around_peak=500_000, between_pheno_mask_around_peak=1_000_000, manhattan_peak_max_count=500, manhattan_peak_pval_threshold=1e-6, manhattan_peak_sprawl_dist=200_000, manhattan_peak_variant_counting_pval_threshold=5e-8):
    manhattan_file = f"{phenocode}_manhattan.json"
    qq_file = f"{phenocode}_qq.json"
    
    # Read data for Manhattan plot & add rsIDs to variants using the LMDB path
    start = time.time()
    reader_for_manhattan = sniffers.guess_gwas_standard(gwas_file).add_filter('neg_log_pvalue')
    print(f"Time for creating the reader: {time.time() - start} seconds", file=sys.stderr)
    start = time.time()
    if lmdb_path is not None:
        add_rsID_with_lmdb(reader_for_manhattan, lmdb_path)
    print(f"Time for adding rsIDs from LMDB: {time.time() - start} seconds", file=sys.stderr)
    start = time.time()
    generate_manhattan(reader_for_manhattan, manhattan_file, manhattan_num_unbinned, within_pheno_mask_around_peak, between_pheno_mask_around_peak, manhattan_peak_max_count, manhattan_peak_pval_threshold, manhattan_peak_sprawl_dist, manhattan_peak_variant_counting_pval_threshold)
    print(f"Time for generating Manhattan plot: {time.time() - start} seconds", file=sys.stderr)
    
    # Read data for QQ plot
    reader_for_qq = sniffers.guess_gwas_standard(gwas_file).add_filter('neg_log_pvalue')
    generate_qq(reader_for_qq, qq_file)
    return manhattan_file, qq_file

def main():
    parser = argparse.ArgumentParser(description="Generate Manhattan and QQ JSON files from a GWAS .gz file")
    parser.add_argument("--input-files", required=True, help="TSV with two columns: <gwas_file.gz> <phenocode>")
    parser.add_argument("--lmdb", help="Path to the LMDB database for variant ID mapping")
    parser.add_argument("--max-workers", type=int, default=4, help="Maximum number of parallel workers")
    
    parser.add_argument("--manhattan-num-unbinned", type=int, default=500, help="Number of unbinned variants for Manhattan plot")
    parser.add_argument("--within-pheno-mask-around-peak", type=int, default=500_000, help="Within-phenotype mask around peak for Manhattan plot")
    parser.add_argument("--between-pheno-mask-around-peak", type=int, default=1_000_000, help="Between-phenotype mask around peak for Manhattan plot")
    parser.add_argument("--manhattan-peak-max-count", type=int, default=500, help="Max count of peaks for Manhattan plot")
    parser.add_argument("--manhattan-peak-pval-threshold", type=float, default= 1e-6, help="P-value threshold for peaks in Manhattan plot")
    parser.add_argument("--manhattan-peak-sprawl-dist", type=int, default=200_000, help="Sprawl distance for peaks in Manhattan plot")
    parser.add_argument("--manhattan-peak-variant-counting-pval-threshold", type=float, default=5e-8, help="P-value threshold for counting variants in peaks for Manhattan plot")

    args = parser.parse_args()

    input_files_dt = pd.read_csv(args.input_files, sep="\t", header=None, names=["file_path", "phenocode"])
    with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        futures = [
            executor.submit(
                process_file,
                row["file_path"],
                row["phenocode"],
                args.lmdb,
                args.manhattan_num_unbinned,
                args.within_pheno_mask_around_peak,
                args.between_pheno_mask_around_peak,
                args.manhattan_peak_max_count,
                args.manhattan_peak_pval_threshold,
                args.manhattan_peak_sprawl_dist,
                args.manhattan_peak_variant_counting_pval_threshold,
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