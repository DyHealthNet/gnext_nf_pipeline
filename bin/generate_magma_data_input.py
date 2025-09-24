#!/usr/bin/env python3
import argparse
import os
import time

from zorp import sniffers, lookups


def main():
    parser = argparse.ArgumentParser(
        description="Generate a MAGMA input file (rsid + pval) from a normalized GWAS file"
    )
    parser.add_argument("--input", required=True, help="Normalized GWAS file (.gz)")
    parser.add_argument("--phenocode", required=True, help="Phenocode for naming output files")
    parser.add_argument("--lmdb",required=True,help="Path to the LMDB directory containing data.mdb")

    args = parser.parse_args()

    # Initialize reader
    reader = sniffers.guess_gwas_standard(args.input).add_filter('neg_log_pvalue')
    rsid_finder = lookups.SnpToRsid(args.lmdb, test=False)
    reader.add_lookup('rsid', lambda variant: rsid_finder(variant.chrom, variant.pos, variant.ref, variant.alt))
    reader.write(args.phenocode + ".txt", columns=["rsid", "pval"], make_tabix=False)

if __name__ == "__main__":
    main()