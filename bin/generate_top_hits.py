#!/usr/bin/env python3
"""
Generate a single global "top hits" table across all phenotypes.

For every per-phenotype Manhattan JSON produced upstream, this script keeps the
strongest variant per peak (below a p-value cutoff), annotates it with the trait
metadata and its nearest genes (looked up from the variant->gene LMDB), then merges
everything into one JSON sorted by signal strength. Consumed by the front-end to
show the most significant associations across the whole study.
"""
import os
import re
import argparse
import json
import pandas as pd
import lmdb, struct, msgpack

class LMDBGeneQuery:
    """
    Context manager wrapping a read-only LMDB environment for variant->gene lookups.

    The environment is opened once on __enter__ and reused across many
    get_genes_for_variant() calls, avoiding the cost of reopening the database for
    every variant.
    """
    def __init__(self, lmdb_path):
        self.lmdb_path = lmdb_path
        self.env = None
    
    def __enter__(self):
        self.env = lmdb.open(self.lmdb_path, readonly=True, subdir=False, lock=False, max_dbs=128)
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.env:
            self.env.close()
    
    def get_genes_for_variant(self, chrom, pos):
        """
        Look up the genes mapped to a single variant position.

        Returns an {ensg_id: symbol} dict ordered nearest-first, or {} when the
        position has no mapped genes / the chromosome DB is missing. The environment
        stays open, but per-chromosome DB handles are not cached (opened per call).
        """
        try:
            # Ensure chrom is a string and normalize
            chrom = str(chrom).replace('chr', '')
            
            # Ensure pos is an integer
            pos = int(pos)
            
            # Create key from position (4-byte unsigned big-endian integer)
            key_bytes = struct.pack(">I", pos)
            
            # Open database and query in same transaction (don't cache handles)
            with self.env.begin(write=False) as txn:
                try:
                    db = self.env.open_db(chrom.encode(), txn=txn)
                    value_bytes = txn.get(key_bytes, db=db)
                except lmdb.NotFoundError:
                    logger.warning(f"Chromosome database '{chrom}' not found in LMDB")
                    return []
            
            # Unpack result
            if value_bytes:
                genes = msgpack.unpackb(value_bytes, raw=False)
                # Sort by distance (already sorted, but explicit for clarity)
                genes_sorted = sorted(genes, key=lambda x: x[2])
                # Return dictionary mapping ensg_id to symbol
                return {ensg_id:symbol for ensg_id, symbol, distance in genes_sorted}
            return {}
        except lmdb.Error as e:
            logger.error(f"LMDB error for {chrom}:{pos} - {e}")
            return {}
        except Exception as e:
            logger.error(f"Error getting genes for variant {chrom}:{pos} - {e}")
            return {}


def get_hits(manhattan_file, pval_cutoff, pheno_dt, lmdb_gene_file):
    """
    Collect the best (lowest p-value) variant per peak for one Manhattan JSON file.

    Args:
        manhattan_file: path to a "<phenocode>_manhattan.json" file.
        pval_cutoff: only variants at or below this p-value are considered.
        pheno_dt: DataFrame mapping phenocode -> category/description metadata.
        lmdb_gene_file: variant->gene LMDB used to annotate nearest genes.

    Returns a list of variant dicts, each enriched with trait metadata, a display
    label ("top_variant") and a "nearest_genes" mapping.
    """
    # Recover the phenocode from the filename ("<phenocode>_manhattan.json").
    phenocode = re.sub(r'_manhattan\.json$', '', os.path.basename(manhattan_file))
    match = pheno_dt.loc[pheno_dt['phenocode'] == phenocode]
    if not match.empty:
        trait_group = match['category'].values[0]
        trait_label = match['description'].values[0]
    else:
        trait_group, trait_label = None, None

    with open(manhattan_file) as f:
        variants = json.load(f)["unbinned_variants"]

    # Keep only the single strongest variant per (chrom, pos) peak.
    peak_to_best = {}
    with LMDBGeneQuery(lmdb_gene_file) as gene_query:
        for v in variants:
            # Only consider peak-flagged variants that clear the p-value cutoff.
            if v.get("pvalue", 1.0) <= pval_cutoff and v.get("peak", False):
                key = (v["chrom"], v["pos"])
                best = peak_to_best.get(key)
                # Replace the stored hit only when this variant is more significant.
                if best is None or v["pvalue"] < best["pvalue"]:
                    v["trait_id"] = phenocode
                    v["trait_group"] = trait_group
                    v["trait_label"] = trait_label
                    chrom = v["chrom"]
                    pos = v["pos"]
                    ref = v.get("ref", "")
                    alt = v.get("alt", "")
                    # Build a human-readable variant label, appending the rsID if known.
                    rsid = v.get("rsid")
                    if rsid and rsid != ".":
                        v["top_variant"] = f"{chrom}_{pos}_{ref}/{alt} ({rsid})"
                    else:
                        v["top_variant"] = f"{chrom}_{pos}_{ref}/{alt}"
                    v["nearest_genes"] = gene_query.get_genes_for_variant(chrom, pos)
                    peak_to_best[key] = v

    return list(peak_to_best.values())

def generate_top_hits(manhattan_files, phenocode_file,outpath, pval_cutoff=1e-6, max_limit=10000, lmdb_gene_file=None):
    """
    Aggregate top hits across every Manhattan file and write the merged JSON.

    Collects per-peak best variants from all phenotypes, ranks them by signal
    strength, truncates to max_limit, and dumps the result to outpath.
    """
    hits = []
    # phenocode -> trait metadata table shared across all files.
    pheno_dt = pd.read_csv(phenocode_file, sep=",")

    for mf in manhattan_files:
        hits.extend(get_hits(mf, pval_cutoff, pheno_dt, lmdb_gene_file))

    # Sort by neg_log_pvalue descending (i.e. strongest signals first) and cap the count.
    hits.sort(key=lambda h: h.get("neg_log_pvalue", 0), reverse=True)
    hits = hits[:max_limit]

    with open(outpath, "w") as f:
        json.dump(hits, f, indent=2)

    print(f"Wrote {len(hits)} top hits to {outpath}")

def main():
    """Parse CLI arguments and drive top-hit aggregation."""
    ap = argparse.ArgumentParser(description="Generate global top hits across manhattan JSON files")
    ap.add_argument("--manhattan-files-file", required=True, help="List of manhattan JSON files (phenocode_manhattan.json)")
    ap.add_argument("--phenocode-file", required=True, help="TSV file mapping phenocodes to descriptions")
    ap.add_argument("--out", default="top_hits.json", help="Output JSON file")
    ap.add_argument("--pval-cutoff", type=float, default=1e-6, help="p-value threshold for top hits")
    ap.add_argument("--max-limit", type=int, default=10000, help="Maximum number of hits to keep")
    ap.add_argument("--lmdb-gene-file", required=True, help="Path to variant-gene LMDB file")
    args = ap.parse_args()

    # The list of Manhattan JSON paths is itself passed as a one-column file.
    manhattan_df = pd.read_csv(args.manhattan_files_file, header=None)
    manhattan_files = manhattan_df[0].tolist()

    generate_top_hits(manhattan_files, args.phenocode_file, args.out, args.pval_cutoff, args.max_limit, args.lmdb_gene_file)

if __name__ == "__main__":
    main()