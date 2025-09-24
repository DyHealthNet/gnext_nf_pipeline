#!/usr/bin/env python3
import argparse
import os
import re
import sys
import traceback
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

from zorp import sniffers, parsers

def process_and_normalize(file_path,
                          phenocode,
                          chr_col, pos_col, ref_col, alt_col,
                          pval_col, pval_neglog10,
                          beta_col, se_col, af_col):
    """
    Normalize a single GWAS summary statistics file.
    Always writes into the current working directory (CWD).
    """
        
    # Build normalized filename without extensions
    norm_filepath = os.path.join(os.getcwd(), str(phenocode))
    
    parser_options = {
        "chrom_col": chr_col,
        "pos_col": pos_col,
        "ref_col": ref_col,
        "alt_col": alt_col,
        "pval_col": pval_col,
        "is_neg_log_pvalue": bool(pval_neglog10),
        "beta_col": beta_col,
        "stderr_beta_col": se_col,
        "allele_freq_col": af_col,
        "rsid": None,
    }
    
    parser = parsers.GenericGwasLineParser(**parser_options)
    
    reader = sniffers.guess_gwas_generic(file_path, parser=parser, skip_errors=True)

    columns = [
        "chrom", "pos", "rsid", "ref", "alt",
        "neg_log_pvalue", "pvalue", "beta", "stderr_beta",
        "alt_allele_freq"
    ]

    reader.write(norm_filepath, make_tabix=True, columns=columns)

    return norm_filepath + ".gz"


def main():
    parser = argparse.ArgumentParser(description="Normalize a single GWAS summary statistics file")
    parser.add_argument("--input-files", required=True, help="TSV with two columns: <gwas_file> <phenocode>")
    parser.add_argument("--max-workers", type=int, default=4, help="Maximum number of parallel workers")
    parser.add_argument("--chr-col", type=int, required=True)
    parser.add_argument("--pos-col", type=int, required=True)
    parser.add_argument("--ref-col", type=int, required=True)
    parser.add_argument("--alt-col", type=int, required=True)
    parser.add_argument("--pval-col", type=int, required=True)
    parser.add_argument("--pval-neglog10", action="store_true")
    parser.add_argument("--beta-col", type=int, required=True)
    parser.add_argument("--se-col", type=int, required=True)
    parser.add_argument("--af-col", type=int, required=True)

    args = parser.parse_args()
    
    input_files_dt = pd.read_csv(args.input_files, sep="\t", header=None, names=["file_path", "phenocode"])
    with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        futures = [
            executor.submit(process_and_normalize, row['file_path'], row['phenocode'], args.chr_col, args.pos_col, args.ref_col, args.alt_col, args.pval_col, args.pval_neglog10, args.beta_col, args.se_col, args.af_col)
            for idx, row in input_files_dt.iterrows()
        ]
        for future in as_completed(futures):
            try:
                normalized_file = future.result()
                print(f"[OK] Generated {normalized_file}", file=sys.stderr)
            except Exception as e:
                print(f"[ERROR] {e}", file=sys.stderr)


if __name__ == "__main__":
    main()