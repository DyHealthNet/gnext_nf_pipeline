#!/usr/bin/env python3
import argparse
import pandas as pd
import sys
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(
        description="Parse MAGMA gene output and optionally merge with gene location file."
    )
    parser.add_argument(
        "--magma-output", required=True,
        help="Path to MAGMA .genes.out file (e.g. phenotype_magma.genes.out)"
    )
    parser.add_argument(
        "--gene-location", required=True,
        help="Optional TSV file containing gene location information (must have 'GENE' column)"
    )
    args = parser.parse_args()

    magma_path = Path(args.magma_output)

    # Read MAGMA gene results (whitespace-separated)
    try:
        df = pd.read_csv(magma_path, sep=r"\s+", comment="#")
    except Exception as e:
        sys.exit(f"Error reading MAGMA output: {e}")

    print(f"Loaded MAGMA results: {df.shape[0]} genes, {df.shape[1]} columns")

    # Merge with gene-location file
    try:
        loc = pd.read_csv(args.gene_location, sep="\t", header = None, names=["GENE", "CHR", "START", "END", "STRAND", "SYMBOL"])
        loc = loc.loc[:, ["GENE","SYMBOL"]]
    except Exception as e:
        sys.exit(f"Error reading gene location file: {e}")

    if "GENE" not in loc.columns:
        sys.exit("Error: gene location file must contain a 'GENE' column.")

    df = pd.merge(df, loc, on="GENE", how="left")
    print(f"Merged with gene locations: {df.shape[1]} columns total")

    # Save parsed file
    df.to_csv(magma_path, sep="\t", index=False)
    print(f"Saved parsed file to: {magma_path}")

if __name__ == "__main__":
    main()