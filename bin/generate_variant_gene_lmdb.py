#!/usr/bin/env python3
import argparse
import os
import struct
import lmdb
import msgpack
import pysam
import sys
from collections import defaultdict
import time


def load_gene_annotations(gene_file):
    """
    Load gene annotations from TSV file.
    Expected columns: ensg_id, chr, start, end, strand, symbol
    Returns dict: {chr: [gene_objects]}
    """
    genes_by_chr = defaultdict(list)
    
    print(f"[INFO] Loading gene annotations from: {gene_file}", file=sys.stderr)
    
    with open(gene_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 6:
                print(f"[WARNING] Skipping line {line_num}: insufficient columns", file=sys.stderr)
                continue
            
            ensg_id, chrom, start, end, strand, symbol = parts[:6]
            
            # Normalize chromosome
            chrom = chrom.replace('chr', '')
            
            try:
                gene = {
                    'ensg_id': ensg_id,
                    'chr': chrom,
                    'start': int(start),
                    'end': int(end),
                    'strand': strand,
                    'symbol': symbol
                }
                genes_by_chr[chrom].append(gene)
            except ValueError as e:
                print(f"[WARNING] Skipping line {line_num}: {e}", file=sys.stderr)
                continue
    
    # Sort genes by start position for efficient searching
    for chrom in genes_by_chr:
        genes_by_chr[chrom].sort(key=lambda g: g['start'])
    
    total_genes = sum(len(genes) for genes in genes_by_chr.values())
    print(f"[INFO] Loaded {total_genes:,} genes across {len(genes_by_chr)} chromosomes", file=sys.stderr)
    
    return genes_by_chr


def get_gene_window(gene, window_up, window_down):
    """
    Calculate the search window for a gene based on strand.
    For + strand: upstream = before start, downstream = after end
    For - strand: upstream = after end, downstream = before start
    Returns (window_start, window_end)
    """
    if gene['strand'] == '+':
        # Positive strand: upstream is before gene, downstream is after
        window_start = gene['start'] - window_up
        window_end = gene['end'] + window_down
    elif gene['strand'] == '-':
        # Negative strand: upstream is after gene, downstream is before
        window_start = gene['start'] - window_down
        window_end = gene['end'] + window_up
    else:
        # Unknown strand: use symmetric window
        window_start = gene['start'] - max(window_up, window_down)
        window_end = gene['end'] + max(window_up, window_down)
    
    return max(0, window_start), window_end


def find_genes_for_variant(chrom, pos, genes_by_chr, window_up, window_down):
    """
    Find all genes within the window of a variant using binary search.
    Returns list of (ensg_id, symbol, distance) tuples.
    """
    if chrom not in genes_by_chr:
        return []
    
    genes = genes_by_chr[chrom]
    nearby_genes = []
    
    # Binary search to find starting point
    # Search for genes whose end >= pos - max_window
    max_window = max(window_up, window_down)
    search_start = pos - max_window
    
    # Find first gene that could overlap with variant window
    left, right = 0, len(genes)
    while left < right:
        mid = (left + right) // 2
        if genes[mid]['end'] < search_start:
            left = mid + 1
        else:
            right = mid
    
    # Check genes starting from the found position
    for i in range(left, len(genes)):
        gene = genes[i]
        
        # If gene starts too far after variant, we can stop
        if gene['start'] > pos + max_window:
            break
        
        window_start, window_end = get_gene_window(gene, window_up, window_down)
        
        # Check if variant falls within the gene's window
        if window_start <= pos <= window_end:
            # Calculate distance from variant to gene body
            if pos < gene['start']:
                distance = gene['start'] - pos
            elif pos > gene['end']:
                distance = pos - gene['end']
            else:
                distance = 0  # Variant is within gene body
            
            nearby_genes.append((gene['ensg_id'], gene['symbol'], distance))
    
    # Sort by distance (closest first)
    nearby_genes.sort(key=lambda x: x[2])
    
    return nearby_genes


def build_snp_gene_map_lmdb_from_vcf(vcf_path, gene_file, window_up, window_down, out_file, num_chroms=28):
    """
    Build an LMDB mapping of (chrom,pos,ref/alt) -> nearest genes from an annotated VCF file.
    """
    # Load gene annotations
    genes_by_chr = load_gene_annotations(gene_file)
    
    print(f"[INFO] Creating LMDB at: {out_file}", file=sys.stderr)
    print(f"[INFO] Window sizes: upstream={window_up}, downstream={window_down}", file=sys.stderr)
    
    env = lmdb.open(out_file, map_size=1024 ** 4, max_dbs=num_chroms, subdir=False)
    db_handles = {}

    vcf = pysam.VariantFile(vcf_path)
    
    total_variants = 0
    mapped_variants = 0
    start_time = time.time()
    batch_size = 10000
    batch_data = []
    
    # Track unique genes that have been mapped (store full tuples)
    mapped_genes_set = set()
    
    # Process variants and create databases on-demand
    for rec in vcf.fetch():
        chrom = rec.chrom.replace('chr', '')
        
        # Skip chromosomes not in gene annotation
        if chrom not in genes_by_chr:
            continue
        
        # Create database handle on-demand (lazily)
        if chrom not in db_handles:
            with env.begin(write=True) as txn:
                db_handles[chrom] = env.open_db(chrom.encode(), txn=txn)
        
        pos = rec.pos
        
        # Find nearby genes for this position (independent of ref/alt)
        nearby_genes = find_genes_for_variant(chrom, pos, genes_by_chr, window_up, window_down)
        
        total_variants += 1
        
        if nearby_genes:
            # Create key: position as 4-byte big-endian int
            key_bytes = struct.pack(">I", pos)
            # Pack gene data as list of tuples: [(ensg_id, symbol, distance), ...]
            value_bytes = msgpack.packb(nearby_genes, use_bin_type=True)
            batch_data.append((chrom, key_bytes, value_bytes))
            mapped_variants += 1
            
            # Track genes that were mapped - add the tuples
            for ensg_id, symbol, distance in nearby_genes:
                mapped_genes_set.add((ensg_id, symbol, chrom))
        
        # Write batch when it reaches batch_size
        if len(batch_data) >= batch_size:
            with env.begin(write=True) as txn:
                for chrom, key, value in batch_data:
                    txn.put(key, value, db=db_handles[chrom])
            batch_data.clear()
        
        if total_variants % 100000 == 0:
            elapsed = time.time() - start_time
            rate = total_variants / elapsed if elapsed > 0 else 0
            print(f"[PROGRESS] Processed {total_variants:,} variants ({mapped_variants:,} mapped) - {rate:.0f} var/sec", file=sys.stderr)
    
    # Write remaining batch
    if batch_data:
        with env.begin(write=True) as txn:
            for chrom, key, value in batch_data:
                txn.put(key, value, db=db_handles[chrom])

    env.sync()
    env.close()
    
    # Write gene metadata to TSV file
    print(f"[INFO] Writing metadata for {len(mapped_genes_set):,} unique mapped genes", file=sys.stderr)
    gene_metadata_file = "mapped_genes.tsv"
    
    # Create a lookup dict for quick gene info retrieval
    gene_lookup = {}
    for chrom, gene_list in genes_by_chr.items():
        for gene in gene_list:
            gene_lookup[gene['ensg_id']] = gene
    
    with open(gene_metadata_file, 'w') as f:
        # Write header
        f.write("ensg_id\tsymbol\tchr\tstart\tend\tstrand\n")
        
        # Write gene info for all mapped genes
        genes_written = 0
        for ensg_id, symbol, chrom in sorted(mapped_genes_set):
            gene = gene_lookup.get(ensg_id)
            if gene:
                f.write(f"{ensg_id}\t{symbol}\t{chrom}\t{gene['start']}\t{gene['end']}\t{gene['strand']}\n")
                genes_written += 1
    
    elapsed = time.time() - start_time
    pct_mapped = (mapped_variants / total_variants * 100) if total_variants else 0
    
    print(f"[DONE] Processed {total_variants:,} variants in {elapsed:.2f}s", file=sys.stderr)
    print(f"[DONE] Mapped {mapped_variants:,} variants to genes ({pct_mapped:.2f}%)", file=sys.stderr)
    print(f"[DONE] Wrote metadata for {genes_written:,} unique genes to: {gene_metadata_file}", file=sys.stderr)
    print(f"[DONE] LMDB created at: {out_file}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Build variant-gene LMDB mapping from annotated VCF"
    )
    parser.add_argument("--vcf", required=True, help="Path to annotated VCF (bgzipped + indexed)")
    parser.add_argument("--gene-file", required=True, help="Path to gene annotation file (TSV)")
    parser.add_argument("--window-up", type=int, default=10, help="Upstream window size in kb")
    parser.add_argument("--window-down", type=int, default=10, help="Downstream window size in kb")
    parser.add_argument("--out-file", required=True, help="Output directory for LMDB")
    parser.add_argument("--num-chroms", type=int, default=28, help="Max number of chromosome DBs")
    args = parser.parse_args()

    # Convert kb to bp
    window_up_bp = args.window_up * 1000
    window_down_bp = args.window_down * 1000

    build_snp_gene_map_lmdb_from_vcf(args.vcf, args.gene_file, window_up_bp, window_down_bp, args.out_file, num_chroms=args.num_chroms)

if __name__ == "__main__":
    main()