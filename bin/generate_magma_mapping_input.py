#!/usr/bin/env python3
import argparse
import gzip
import os
from collections import defaultdict
import csv
import logging


def extract_csq_fields(vcf_path):
    """
    Extracts the CSQ field names from a VCF file.
    """
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith("##INFO") and "Format:" in line and "CSQ" in line:
                parts = line.strip().split("|")
                parts[0] = parts[0].split("Format: ")[1]
                parts[-1] = parts[-1].strip('">')
                return parts
    raise ValueError("No CSQ format found in VCF header.")

def prepare_config_dicts(config_file):
    # Prepare config dicts in one pass + open headers
    config_dicts = []
    with open(config_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            row = {k.lower(): v for k, v in row.items()}
            mapping_strategy = (row.get("mapping_strategy") or "").lower()
            if mapping_strategy != "positional":
                logging.warning(f"Skipping config row (unsupported strategy): {row}")
                continue

            up = int(row.get("window_up", 0) or 0)
            down = int(row.get("window_down", 0) or 0)
            out_file = f"magma_{mapping_strategy}_{up}_{down}_annotations.genes.annot"
            tmp_file = f"magma_{mapping_strategy}_{up}_{down}_tmp.txt"

            # Write header immediately
            with open(out_file, 'w') as f_out:
                f_out.write(f"# window_up = {up}\n")
                f_out.write(f"# window_down = {down}\n")

            config_dicts.append({
                "window_up": up,
                "window_down": down,
                "out_file": out_file,
                "tmp_file": tmp_file,
                "gene_to_rsids": set()
                })
    if not config_dicts:
        logging.error("No valid positional configs found. Exiting.")
        return None
    return config_dicts

def flush_config(conf):
    """Flush accumulated gene–rsid pairs to file and clear the set."""
    if not conf["gene_to_rsids"]:
        return
    with open(conf["tmp_file"], 'a') as f_out:
        for gene, rsid in conf["gene_to_rsids"]:
            f_out.write(f"{gene}\t{rsid}\n")
    conf["gene_to_rsids"].clear()


def prepare_magma_mapping(vcf_path, config_file, flush_limit=500000):

    without_gene_terms = {
        "regulatory_region_variant",
        "TF_binding_site_variant",
        "intergenic_variant",
        "intron_variant"
    }

    config_dicts = prepare_config_dicts(config_file)
    if not config_dicts:
        return

    csq_fields = extract_csq_fields(vcf_path)

    # Stream VCF once
    i = 0
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            info = parts[7]
            if "CSQ=" not in info:
                continue

            csq_data = info.split("CSQ=")[1].split(";")[0]
            entries = csq_data.split(',')
            i += 1
            for entry in entries:
                values = entry.split('|')
                csq_dict = dict(zip(csq_fields, values))

                if csq_dict.get("BIOTYPE") != "protein_coding" or csq_dict.get("Feature_type") != "Transcript":
                    continue

                consequences = csq_dict.get("Consequence", "").split("&")
                gene = csq_dict.get("Gene")
                ids = csq_dict.get("Existing_variation", "").split("&")
                rsid = next((j for j in ids if j.startswith("rs")), None)

                if not gene or not rsid:
                    continue

                # For each config, check rules
                for conf in config_dicts:
                    if any(c not in without_gene_terms for c in consequences):
                        if any(c not in ("upstream_gene_variant", "downstream_gene_variant") for c in consequences):
                            if (rsid, gene) not in conf["gene_to_rsids"]:
                                conf["gene_to_rsids"].add((gene, rsid))
                            continue

                        distance = csq_dict.get("DISTANCE")
                        try:
                            distance_int = int(distance)
                        except (ValueError, TypeError):
                            continue

                        if ("upstream_gene_variant" in consequences
                            and conf["window_up"]
                            and distance_int <= conf["window_up"]
                            and (rsid, gene) not in conf["gene_to_rsids"]):
                            conf["gene_to_rsids"].add((gene, rsid))
                        elif ("downstream_gene_variant" in consequences
                              and conf["window_down"]
                              and distance_int <= conf["window_down"]
                              and (rsid, gene) not in conf["gene_to_rsids"]):
                            conf["gene_to_rsids"].add((gene, rsid))
                    if len(conf["gene_to_rsids"]) >= flush_limit:
                        flush_config(conf)

    # Final flush into tmp files
    logging.info("Final flush...")
    for conf in config_dicts:
        flush_config(conf)
        
    # Create final mapping file list
    for conf in config_dicts:
        # Read tmp file from disk and create dict of gene -> set(rsids)
        gene_to_rsids = defaultdict(set)
        if not os.path.exists(conf["tmp_file"]):
            logging.warning(f"No variants mapped for config: {conf}. Skipping.")
            continue
        with open(conf["tmp_file"], 'r') as f_tmp:
            for line in f_tmp:
                gene, rsid = line.strip().split('\t')
                gene_to_rsids[gene].add(rsid)
        os.remove(conf["tmp_file"])
        # Write final file
        with open(conf["out_file"], 'a') as f_out:
            for gene in sorted(gene_to_rsids.keys()):
                rsids = sorted(gene_to_rsids[gene])
                f_out.write(f"{gene}\t1:1:2\t{' '.join(rsids)}\n")
        logging.info(f"Wrote MAGMA mapping file: {conf['out_file']}")


def main():
    parser = argparse.ArgumentParser(description="Prepare MAGMA mapping input files from annotated VCF.")
    parser.add_argument("--vcf", required=True, help="Annotated VCF file (bgzipped).")
    parser.add_argument("--config", required=True, help="MAGMA config file (TSV with mapping_strategy, window_up, window_down).")
    args = parser.parse_args()
    
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")


    prepare_magma_mapping(args.vcf, args.config)


if __name__ == "__main__":
    main()