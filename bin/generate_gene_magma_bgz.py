#!/usr/bin/env python3
"""
Generate a single BGZ file with MAGMA gene p-values across all traits.
Each row represents a gene, with columns for each trait's p-value.
Also creates an LMDB index mapping ENSG_ID to byte offset for O(1) lookups.

Input manifest file format (tab-separated):
    phenocode1<tab>/path/to/trait1_magma.genes.out
    phenocode2<tab>/path/to/trait2_magma.genes.out
    ...

Example:
    pineapple_all	/data/magma/pineapple_all_magma.genes.out
    cinnamon_all	/data/magma/cinnamon_all_magma.genes.out

Usage:
    ./generate_gene_magma_bgz.py \
        --manifest manifest.tsv \
        --gene-file mapped_genes.tsv \
        --output gene_pvalues.tsv.bgz \
        --index gene_magma_index.lmdb
"""
import os
import argparse
import sys
import numpy as np
import pandas as pd
import pysam
import lmdb
from collections import OrderedDict

def load_gene_info(gene_file):
    """
    Load gene information from a unique gene file.
    Expected format: ENSG_ID<tab>Symbol (header or no header)
    Returns DataFrame with GENE and SYMBOL columns.
    """
    # Try reading with header first
    try:
        df = pd.read_csv(gene_file, sep="\t")
    except Exception as e:
        print(f"Error reading gene file: {e}", file=sys.stderr)
        raise
    
    # Keep only GENE and SYMBOL
    gene_info = df[['ensg_id', 'symbol']].copy()
    gene_info.columns = ['GENE', 'SYMBOL']
    print(f"Loaded {len(gene_info)} genes from {gene_file}", file=sys.stderr)
    return gene_info

def load_manifest(manifest_file):
    """
    Load manifest file with format: phenocode\tmagma_file_path
    Returns: (trait_names, magma_files)
    """
    traits = []
    files = []
    
    with open(manifest_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) != 2:
                print(f"Warning: Invalid line format (expected 2 columns): {line}", file=sys.stderr)
                continue
            
            phenocode, magma_file = parts
            traits.append(phenocode)
            files.append(magma_file)
    
    return traits, files

def load_pvalues_matrix(gene_info, trait_names, magma_files):
    """
    Load p-values from all MAGMA files into a matrix.
    Rows = genes (from gene_info), Columns = traits
    """
    ngenes = len(gene_info)
    ntraits = len(trait_names)
    
    # Create empty matrix filled with NaN
    pval_matrix = np.full((ngenes, ntraits), np.nan, dtype=np.float64)
    
    # Create gene index mapping
    gene_to_idx = {gene: i for i, gene in enumerate(gene_info['GENE'])}
    
    # Load p-values for each trait
    for j, (trait, magma_file) in enumerate(zip(trait_names, magma_files)):
        if not os.path.exists(magma_file):
            print(f"Warning: Skipping missing file for trait {trait}", file=sys.stderr)
            continue
        
        try:
            df = pd.read_csv(magma_file, sep="\t")
            
            for _, row in df.iterrows():
                gene = row['GENE']
                pval = row['P']
                
                if gene in gene_to_idx:
                    gene_idx = gene_to_idx[gene]
                    pval_matrix[gene_idx, j] = pval
            
            print(f"Loaded p-values for trait {trait} ({j+1}/{ntraits})", file=sys.stderr)
        except Exception as e:
            print(f"Warning: Error loading {magma_file}: {e}", file=sys.stderr)
            continue
    
    return pval_matrix

def write_bgz_file_with_index(gene_info, trait_names, pval_matrix, output_file, index_file):
    """
    Write the gene p-value matrix to a BGZ file and create LMDB index mapping ENSG_ID to virtual offset.
    Uses pysam's tell() to get BGZF virtual offsets for random access.
    """    
    print(f"Writing BGZ file: {output_file}", file=sys.stderr)
    print(f"Creating LMDB index: {index_file}", file=sys.stderr)
    
    # Create LMDB environment
    # map_size: 1GB should be enough for index (much smaller than data)
    # subdir=False creates lmdb-data.mdb and lmdb-data.mdb-lock files directly
    env = lmdb.open(index_file, map_size=1024*1024*1024, max_dbs=1, subdir=False)
    
    with env.begin(write=True) as txn:
        with pysam.BGZFile(output_file, "w") as bgz:
            # Write header
            header = "#" + "\t".join(["GENE", "SYMBOL"] + trait_names) + "\n"
            header_bytes = header.encode()
            bgz.write(header_bytes)
            
            # Write data rows and build index
            for i in range(len(gene_info)):
                gene = gene_info.iloc[i]['GENE']
                symbol = gene_info.iloc[i]['SYMBOL']
                
                # Get current virtual offset BEFORE writing the line
                # This is the BGZF virtual offset that can be used with seek()
                virtual_offset = bgz.tell()
                
                # Store virtual offset in LMDB (8 bytes, little-endian)
                # Key: ENSG_ID (bytes), Value: virtual offset (8-byte integer)
                txn.put(gene.encode('utf-8'), virtual_offset.to_bytes(8, byteorder='little'))
                
                # Format p-values (use '.' for NaN)
                pvals_str = []
                for pval in pval_matrix[i, :]:
                    if np.isnan(pval):
                        pvals_str.append('.')
                    else:
                        pvals_str.append(f"{pval:.6e}")
                
                line = "\t".join([gene, symbol] + pvals_str) + "\n"
                line_bytes = line.encode()
                bgz.write(line_bytes)
    
    env.sync()
    env.close()
    print(f"Successfully created {output_file} with index {index_file}", file=sys.stderr)
    print(f"Indexed {len(gene_info)} genes", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description="Generate BGZ file with MAGMA gene p-values across all traits"
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="Tab-separated manifest file with format: phenocode<tab>magma_file_path"
    )
    parser.add_argument(
        "--gene-file",
        required=True,
        help="Tab-separated file with gene information: ENSG_ID<tab>Symbol"
    )
    parser.add_argument(
        "--out-bgz-file",
        default="gene_magma_pvalues.tsv.bgz",
        help="Output BGZ file name (default: gene_magma_pvalues.tsv.bgz)"
    )
    parser.add_argument(
        "--out-lmdb-file",
        default="lmdb-data.mdb",
        help="Output LMDB index file name (default: lmdb-data.mdb)"
    )
    
    args = parser.parse_args()
    
    # Load gene information
    print(f"Reading gene file: {args.gene_file}", file=sys.stderr)
    gene_info = load_gene_info(args.gene_file)
    
    # Load manifest file
    print(f"Reading manifest file: {args.manifest}", file=sys.stderr)
    trait_names, magma_paths = load_manifest(args.manifest)
    
    if not magma_paths:
        print("Error: No MAGMA files found in manifest", file=sys.stderr)
        return 1
    
    print(f"Processing {len(magma_paths)} MAGMA files", file=sys.stderr)
    
    # Load p-values into matrix
    pval_matrix = load_pvalues_matrix(gene_info, trait_names, magma_paths)
    
    # Write BGZ file with LMDB index
    write_bgz_file_with_index(gene_info, trait_names, pval_matrix, args.out_bgz_file, args.out_lmdb_file)
    
    print("Done!", file=sys.stderr)
    return 0

if __name__ == "__main__":
    sys.exit(main())
