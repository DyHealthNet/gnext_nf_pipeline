#!/usr/bin/env python3
import os
import re
import argparse
import json

def get_hits(manhattan_file, pval_cutoff):
    """
    Collect best variant per peak for a given manhattan JSON file.
    """
    phenocode = re.sub(r'_manhattan\.json$', '', os.path.basename(manhattan_file))

    with open(manhattan_file) as f:
        variants = json.load(f)["unbinned_variants"]

    peak_to_best = {}
    for v in variants:
        if v.get("pvalue", 1.0) <= pval_cutoff and v.get("peak", False):
            key = (v["chrom"], v["pos"])
            best = peak_to_best.get(key)
            if best is None or v["pvalue"] < best["pvalue"]:
                v["phenocode"] = phenocode
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

def generate_top_hits(manhattan_files, outpath, pval_cutoff=1e-6, max_limit=10000):
    hits = []
    for mf in manhattan_files:
        hits.extend(get_hits(mf, pval_cutoff))

    # Sort by neg_log_pvalue descending (i.e. strongest signals first)
    hits.sort(key=lambda h: h.get("neg_log_pvalue", 0), reverse=True)
    hits = hits[:max_limit]

    with open(outpath, "w") as f:
        json.dump(hits, f, indent=2)

    print(f"Wrote {len(hits)} top hits to {outpath}")

def main():
    ap = argparse.ArgumentParser(description="Generate global top hits across manhattan JSON files")
    ap.add_argument("--manhattan-files", nargs="+", required=True, help="List of manhattan JSON files (phenocode_manhattan.json)")
    ap.add_argument("--out", default="top_hits.json", help="Output JSON file")
    ap.add_argument("--pval-cutoff", type=float, default=1e-6, help="p-value threshold for top hits")
    ap.add_argument("--max-limit", type=int, default=10000, help="Maximum number of hits to keep")
    args = ap.parse_args()

    generate_top_hits(args.manhattan_files, args.out, args.pval_cutoff, args.max_limit)

if __name__ == "__main__":
    main()