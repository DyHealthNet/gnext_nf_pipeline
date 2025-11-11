#!/usr/bin/env python3
import argparse
import sys
import time
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from zorp import sniffers


# -------------------- Helpers --------------------

def normalize_chr(chr_str: str) -> str:
    """Normalize chromosome names to match .bim convention (1–25)."""
    chr_str = str(chr_str).replace("chr", "")
    mapping = {"X": "23", "Y": "24", "M": "25", "MT": "25"}
    return mapping.get(chr_str, chr_str)

def load_bim_ids(bim_path: str, genome_build: str) -> set:
    """Read .bim file and return unique IDs (chr:pos:allele1:allele2)."""
    df = pd.read_csv(bim_path, sep=r"\s+", header=None, names=["Chromosome", "ID", "CM", "Position", "Allele1", "Allele2"])   
    print(f"[INFO] Loaded {len(df):,} variants from BIM file.")
    # Apply normalization to chromosome column
    df["Chromosome"] = df["Chromosome"].apply(normalize_chr)
    return set(df["ID"]), set(df["Chromosome"])

def meta_id_generation(variant):
    """Construct ID like <chr>:<pos>:<allele1>:<allele2> with alphabetically sorted alleles."""
    alleles = ":".join(sorted([variant.ref, variant.alt]))
    return f"{variant.chrom}:{variant.pos}:{alleles}"

def in_mhc_region(chrom: str, pos: int, genome_build: str) -> bool:
    """Check if a variant lies in the MHC region."""
    if genome_build == "GRCh37":
        start, end = 28_477_797, 33_448_354
    else:
        start, end = 28_510_120, 33_480_577
    return chrom == "6" and start <= pos <= end

# -------------------- Main --------------------

def generate_MAGMA_gwas_input_file(gwas_path: str, phenocode: str, bim_ids: set, bim_chr: set, genome_build: str) -> dict:
    """Filter GWAS variants by presence in reference BIM and outside the MHC region."""
    reader = sniffers.guess_gwas_standard(gwas_path).add_filter("neg_log_pvalue").add_filter("pvalue")

    rows, total, retained, chrom_mismatch = [], 0, 0, 0
    start = time.time()
    print(f"[INFO] Reading GWAS: {gwas_path}")

    for variant in reader:

        # Skip MHC region (based on GWAS)
        chrom = normalize_chr(variant.chrom)
        if in_mhc_region(chrom, variant.pos, genome_build):
            continue
        
        # Check if chromosome matching
        if chrom not in bim_chr:
            chrom_mismatch+=1

        gwas_id = f"{chrom}:{variant.pos}:{':'.join(sorted([variant.ref, variant.alt]))}"
        total += 1
        if gwas_id in bim_ids:
            rows.append((gwas_id, variant.pvalue))
            retained += 1

    pct = (retained / total * 100) if total else 0
    print(f"[DONE] Processed {total:,} variants; retained {retained:,} ({pct:.2f}%)")    
    print(f"[TIME] Elapsed: {time.time() - start:.2f}s")
    df = pd.DataFrame(rows, columns=["SNP", "P"])
    df.to_csv(f"{phenocode}_magma.tsv", sep="\t", index=False, header=False)
    # return stats
    return {
        "phenocode": phenocode,
        "input_file": gwas_path,
        "magma_file": f"{phenocode}_magma.tsv",
        "total_variants": total,
        "retained_variants": retained,
        "retained_pct": round(pct,2),
        "chrom_mismatch": chrom_mismatch,
    }



def main():
    parser = argparse.ArgumentParser(
        description="Annotate GWAS variants against a BIM reference using PyRanges."
    )
    parser.add_argument("--input-files", required=True, help="TSV with two columns: <gwas_file.gz> <phenocode>")
    parser.add_argument("--ref-bim", required=True, help="Path to reference .bim file")
    parser.add_argument("--genome-build", choices=["GRCh37", "GRCh38"], required=True)
    parser.add_argument("--max-workers", type=int, default=4, help="Maximum number of parallel workers")

    args = parser.parse_args()

    # Load data
    bim_ids, bim_chr = load_bim_ids(args.ref_bim, args.genome_build)
    all_stats = []
    input_files_dt = pd.read_csv(args.input_files, sep="\t", header=None, names=["file_path", "phenocode"])
    with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        futures = [
            executor.submit(generate_MAGMA_gwas_input_file, row['file_path'], row['phenocode'], bim_ids, bim_chr, args.genome_build)
            for idx, row in input_files_dt.iterrows()
        ]
        for future in as_completed(futures):
            try:
                stats = future.result()
                all_stats.append(stats)
                print(f"[OK] Generated {stats['magma_file']}", file=sys.stderr)
            except Exception as e:
                print(f"[ERROR] {e}", file=sys.stderr)
    
    # Save stats
    stats_df = pd.DataFrame(all_stats)
    phenocodes = input_files_dt["phenocode"]
    stats_df.to_csv("mapping_summary.tsv", sep="\t", index=False)
    print(f"[SUMMARY] Mapping statistics written to mapping_summary.tsv")


if __name__ == "__main__":
    main()