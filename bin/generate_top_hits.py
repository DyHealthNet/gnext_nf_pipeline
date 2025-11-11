#!/usr/bin/env python3
import os
import re
import argparse
import json
import pandas as pd

def get_hits(manhattan_file, pval_cutoff, pheno_dt):
    """
    Collect best variant per peak for a given manhattan JSON file.
    """
    phenocode = re.sub(r'_manhattan\.json$', '', os.path.basename(manhattan_file))
    match = pheno_dt.loc[pheno_dt['phenocode'] == phenocode]
    if not match.empty:
        trait_group = match['category'].values[0]
        trait_label = match['description'].values[0]
    else:
        trait_group, trait_label = None, None

    with open(manhattan_file) as f:
        variants = json.load(f)["unbinned_variants"]

    peak_to_best = {}
    for v in variants:
        if v.get("pvalue", 1.0) <= pval_cutoff and v.get("peak", False):
            key = (v["chrom"], v["pos"])
            best = peak_to_best.get(key)
            if best is None or v["pvalue"] < best["pvalue"]:
                v["trait_id"] = phenocode
                v["trait_group"] = trait_group
                v["trait_label"] = trait_label
                chrom = v["chrom"]
                pos = v["pos"]
                ref = v.get("ref", "")
                alt = v.get("alt", "")
                rsid = v.get("rsid")
                if rsid and rsid != ".":
                    v["top_variant"] = f"{chrom}_{pos}_{ref}/{alt} ({rsid})"
                else:
                    v["top_variant"] = f"{chrom}_{pos}_{ref}/{alt}"
                alt_allele_freq = v.get("alt_allele_freq")
                v["MAF"] = (
                    min(float(alt_allele_freq), 1 - float(alt_allele_freq))
                    if alt_allele_freq not in (None, ".")
                    else None
                )
                peak_to_best[key] = v

    return list(peak_to_best.values())

def generate_top_hits(manhattan_files, phenocode_file,outpath, pval_cutoff=1e-6, max_limit=10000):
    hits = []
    # Read phenocode file
    pheno_dt = pd.read_csv(phenocode_file, sep=",") 
    
    for mf in manhattan_files:
        hits.extend(get_hits(mf, pval_cutoff, pheno_dt))

    # Sort by neg_log_pvalue descending (i.e. strongest signals first)
    hits.sort(key=lambda h: h.get("neg_log_pvalue", 0), reverse=True)
    hits = hits[:max_limit]

    with open(outpath, "w") as f:
        json.dump(hits, f, indent=2)

    print(f"Wrote {len(hits)} top hits to {outpath}")

def main():
    ap = argparse.ArgumentParser(description="Generate global top hits across manhattan JSON files")
    ap.add_argument("--manhattan-files-file", required=True, help="List of manhattan JSON files (phenocode_manhattan.json)")
    ap.add_argument("--phenocode-file", required=True, help="TSV file mapping phenocodes to descriptions")
    ap.add_argument("--out", default="top_hits.json", help="Output JSON file")
    ap.add_argument("--pval-cutoff", type=float, default=1e-6, help="p-value threshold for top hits")
    ap.add_argument("--max-limit", type=int, default=10000, help="Maximum number of hits to keep")
    args = ap.parse_args()
    
    # Read Manhattan files from a file and put as list
    manhattan_df = pd.read_csv(args.manhattan_files_file, header=None)
    manhattan_files = manhattan_df[0].tolist()

    generate_top_hits(manhattan_files, args.phenocode_file, args.out, args.pval_cutoff, args.max_limit)

if __name__ == "__main__":
    main()