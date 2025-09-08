#!/usr/bin/env python3
import argparse
import os
import re
import sys
import traceback

from zorp import sniffers, parsers


def get_bool(value):
    val = str(value).strip().lower()
    if val in {"true", "1", "yes", "on"}:
        return True
    if val in {"false", "0", "no", "off"}:
        return False
    return None


def process_and_normalize(filename,
                          gwas_dir,
                          phenocode,
                          chr_col, pos_col, ref_col, alt_col,
                          pval_col, pval_neglog10,
                          beta_col, se_col, af_col):
    """
    Normalize a single GWAS summary statistics file.
    Always writes into the current working directory (CWD).
    """

    in_filepath = os.path.join(gwas_dir, filename)

    # Build normalized filename without extensions
    norm_filepath = os.path.join(os.getcwd(), phenocode)

    if os.path.exists(norm_filepath + ".gz"):
        print(f"[INFO] Skipping {filename}, already normalized")
        return norm_filepath + ".gz"

    parser_options = {
        "chrom_col": chr_col,
        "pos_col": pos_col,
        "ref_col": ref_col,
        "alt_col": alt_col,
        "pval_col": pval_col,
        "is_neg_log_pvalue": get_bool(pval_neglog10),
        "beta_col": beta_col,
        "stderr_beta_col": se_col,
        "allele_freq_col": af_col,
        "rsid": None,
    }

    parser = parsers.GenericGwasLineParser(**parser_options)
    reader = sniffers.guess_gwas_generic(in_filepath, parser=parser, skip_errors=True)

    columns = [
        "chrom", "pos", "rsid", "ref", "alt",
        "neg_log_pvalue", "pvalue", "beta", "stderr_beta",
        "alt_allele_freq"
    ]

    reader.write(norm_filepath, make_tabix=True, columns=columns)
    print(f"[INFO] Completed normalization: {filename}")

    return norm_filepath + ".gz"


def main():
    parser = argparse.ArgumentParser(description="Normalize a single GWAS summary statistics file")
    parser.add_argument("--input-dir", required=True, help="Directory containing GWAS files")
    parser.add_argument("--filename", required=True, help="GWAS file (relative to input_dir)")
    parser.add_argument("--phenocode", required=True, help="Phenocode for the GWAS file")
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

    try:
        result = process_and_normalize(
            args.filename,
            args.input_dir,
            args.phenocode,
            args.chr_col, args.pos_col, args.ref_col, args.alt_col,
            args.pval_col, args.pval_neglog10,
            args.beta_col, args.se_col, args.af_col,
        )
        print(f"[INFO] Output written: {result}")
    except Exception:
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()