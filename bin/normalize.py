#!/usr/bin/env python3
import argparse
import os
import re
import sys
import traceback
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

from zorp import sniffers, parsers

# Define missing values inline to avoid import issues
MISSING_VALUES = {'', '.', 'NA', 'N/A', 'nan', 'NaN', 'NULL', 'null', 'None'}


class BasicVariantWithSampleSize(parsers.BasicVariant):
    """Extended BasicVariant that includes sample size"""
    __slots__ = ('n_samples',)
    _fields = ('chrom', 'pos', 'rsid', 'ref', 'alt', 'neg_log_pvalue', 'beta', 'stderr_beta', 'alt_allele_freq', 'n_samples')
    
    def __init__(self, chrom, pos, rsid, ref, alt, neg_log_pvalue, beta, stderr_beta, alt_allele_freq, n_samples=None):
        super().__init__(chrom, pos, rsid, ref, alt, neg_log_pvalue, beta, stderr_beta, alt_allele_freq)
        self.n_samples = n_samples


def create_parser_with_sample_size(sample_size_col, **parser_options):
    """
    Create a custom parser that extracts sample size from GWAS files.
    This wraps GenericGwasLineParser and adds sample size extraction.
    """
    # Create the base parser without sample_size_col (which it doesn't support for storage)
    base_parser = parsers.GenericGwasLineParser(**parser_options)
    
    # Store the column index (convert from 1-based to 0-based)
    n_samples_col_idx = sample_size_col - 1
    
    def wrapper(line):
        """Parse line and add sample size"""
        # Parse with base parser
        variant = base_parser(line)
        
        # Extract sample size from the line
        fields = line.strip().split('\t')
        if n_samples_col_idx < len(fields):
            ns_value = fields[n_samples_col_idx]
            if ns_value not in MISSING_VALUES:
                try:
                    n_samples = int(float(ns_value))
                except (ValueError, TypeError):
                    n_samples = None
            else:
                n_samples = None
        else:
            n_samples = None
        
        # Create new variant with sample size
        variant_with_ns = BasicVariantWithSampleSize(
            variant.chrom, variant.pos, variant.rsid,
            variant.ref, variant.alt, variant.neg_log_pvalue,
            variant.beta, variant.stderr_beta, variant.alt_allele_freq,
            n_samples
        )
        return variant_with_ns
    
    return wrapper

def process_and_normalize(file_path,
                          phenocode,
                          chr_col, pos_col, ref_col, alt_col,
                          pval_col, pval_neglog10,
                          beta_col, se_col, af_col, sample_size_col):
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
        "rsid_col": None,
    }
    
    columns = [
        "chrom", "pos", "rsid", "ref", "alt",
        "neg_log_pvalue", "pvalue", "beta", "stderr_beta",
        "alt_allele_freq"
    ]
    
    if sample_size_col is not None:
        # Use custom parser that extracts and stores sample size
        parser = create_parser_with_sample_size(sample_size_col, **parser_options)
        columns.append("n_samples")
    else:
        # Use standard parser
        parser = parsers.GenericGwasLineParser(**parser_options)
    
    reader = sniffers.guess_gwas_generic(file_path, parser=parser, skip_errors=True)

    reader.write(norm_filepath, make_tabix=False, columns=columns)

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
    parser.add_argument("--sample-size-col", type=int, required=False, help="Column index for sample size (optional)")

    args = parser.parse_args()
    
    input_files_dt = pd.read_csv(args.input_files, sep="\t", header=None, names=["file_path", "phenocode"])
    print(f"[INFO] Loaded {len(input_files_dt)} files to normalize", file=sys.stderr)
    with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        futures = [
            executor.submit(process_and_normalize, row['file_path'], row['phenocode'], args.chr_col, args.pos_col, args.ref_col, args.alt_col, args.pval_col, args.pval_neglog10, args.beta_col, args.se_col, args.af_col, args.sample_size_col)
            for idx, row in input_files_dt.iterrows()
        ]
        for future in as_completed(futures):
            try:
                normalized_file = future.result()
                print(f"[OK] Generated {normalized_file}", file=sys.stderr)
            except Exception as e:
                print(f"[ERROR] {e}", file=sys.stderr)
                exit(1)


if __name__ == "__main__":
    main()