#!/usr/bin/env python3
import argparse
import sys
import time
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from zorp import sniffers, parsers


# -------------------- Helpers --------------------

class BasicVariantWithSampleSize(parsers.BasicVariant):
    """Extended BasicVariant that includes sample size"""
    __slots__ = ('n_samples',)
    _fields = ('chrom', 'pos', 'rsid', 'ref', 'alt', 'neg_log_pvalue', 'beta', 'stderr_beta', 'alt_allele_freq', 'n_samples')
    
    def __init__(self, chrom, pos, rsid, ref, alt, neg_log_pvalue, beta, stderr_beta, alt_allele_freq, n_samples=None):
        super().__init__(chrom, pos, rsid, ref, alt, neg_log_pvalue, beta, stderr_beta, alt_allele_freq)
        self.n_samples = n_samples

def create_parser_with_sample_size():
    """Create a parser that reads n_samples from normalized files (column 11)"""
    MISSING_VALUES = {'', '.', 'NA', 'N/A', 'nan', 'NaN', 'NULL', 'null', 'None'}
    
    base_parser = parsers.GenericGwasLineParser(
        chrom_col=1, pos_col=2, ref_col=4, alt_col=5,
        pval_col=7, is_neg_log_pvalue=False,
        beta_col=8, stderr_beta_col=9, allele_freq_col=10,
        rsid_col=3
    )
    
    def wrapper(line):
        """Parse line and add sample size from column 11"""
        variant = base_parser(line)
        
        # Extract n_samples from column 11 (index 10)
        fields = line.strip().split('\t')
        if len(fields) > 10:
            ns_value = fields[10]
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

# -------------------- Main --------------------

def generate_MAGMA_gwas_input_file(gwas_path: str, phenocode: str, bim_ids: set, bim_chr: set, genome_build: str, include_n_samples: bool) -> dict:
    """Filter GWAS variants by presence in reference BIM and outside the MHC region."""
    
    # Choose reader based on whether we need to read n_samples
    if include_n_samples:
        # Use custom parser to read n_samples from normalized files (column 11)
        parser = create_parser_with_sample_size()
        reader = sniffers.guess_gwas_generic(gwas_path, parser=parser, skip_errors=True)
    else:
        # Use standard reader for files without n_samples
        reader = sniffers.guess_gwas_standard(gwas_path).add_filter("neg_log_pvalue").add_filter("pvalue")

    total, retained, chrom_mismatch = 0, 0, 0
    start = time.time()
    output_file = f"{phenocode}_magma.tsv"
    
    print(f"[INFO] Reading GWAS: {gwas_path}")
    print(f"[INFO] Writing to: {output_file}")

    # Open output file and write directly
    with open(output_file, 'w') as out_f:
        for variant in reader:
            # Skip MHC region (based on GWAS)
            chrom = normalize_chr(variant.chrom)
            
            # Check if chromosome matching
            if chrom not in bim_chr:
                chrom_mismatch += 1

            gwas_id = f"{chrom}:{variant.pos}:{variant.ref}:{variant.alt}"
            total += 1
            
            if gwas_id in bim_ids:
                # Write directly to file
                if include_n_samples and hasattr(variant, 'n_samples') and variant.n_samples is not None:
                    out_f.write(f"{gwas_id}\t{variant.pvalue}\t{variant.n_samples}\n")
                else:
                    out_f.write(f"{gwas_id}\t{variant.pvalue}\n")
                retained += 1

    pct = (retained / total * 100) if total else 0
    print(f"[DONE] Processed {total:,} variants; retained {retained:,} ({pct:.2f}%)")    
    print(f"[TIME] Elapsed: {time.time() - start:.2f}s")
    
    # return stats
    return {
        "phenocode": phenocode,
        "input_file": gwas_path,
        "magma_file": output_file,
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
    parser.add_argument("--include-n-samples", action="store_true", help="Include n_samples column in output (if present in input)")

    args = parser.parse_args()

    # Load data
    bim_ids, bim_chr = load_bim_ids(args.ref_bim, args.genome_build)
    all_stats = []
    input_files_dt = pd.read_csv(args.input_files, sep="\t", header=None, names=["file_path", "phenocode"])
    with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        futures = [
            executor.submit(generate_MAGMA_gwas_input_file, row['file_path'], row['phenocode'], bim_ids, bim_chr, args.genome_build, args.include_n_samples)
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