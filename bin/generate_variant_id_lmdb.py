#!/usr/bin/env python3
"""
Build the variant->rsID LMDB used to resolve dbSNP identifiers at query time.

Reads an annotated VCF, extracts the rsID from each variant's VEP "CSQ"
(Existing_variation) field, and stores, per chromosome, a position-keyed mapping
of "ref/alt" -> integer rsID. rsIDs are stored as integers (the "rs" prefix and
"rs" reconstruction happen at lookup time) to keep the database compact.
"""
import argparse
import os
import struct
import lmdb
import msgpack
import pysam


def build_snp_map_lmdb_from_vcf(vcf_path, out_file, num_chroms=28):
    """
    Build the position->{ref/alt: rsID} LMDB from an annotated VCF.

    Creates one LMDB sub-database per chromosome (lazily). For each variant the
    rsID is parsed out of the VEP CSQ annotation; positions with at least one
    resolvable rsID are written with a 4-byte big-endian position key and a
    msgpack-encoded allele->rsID map value.
    """
    env = lmdb.open(out_file, map_size=1024 ** 4, max_dbs=num_chroms, subdir = False)
    db_handles = {}

    vcf = pysam.VariantFile(vcf_path)
    # Locate the rsID ("Existing_variation") column within the VEP CSQ format string.
    csq_header = vcf.header.info["CSQ"].description.split("Format: ")[1].split("|")
    rsid_index = csq_header.index("Existing_variation")

    with env.begin(write=True) as txn:
        for rec in vcf.fetch():
            chrom = rec.chrom

            if chrom not in db_handles:
                db_handles[chrom] = env.open_db(chrom.encode(), txn=txn)

            db = db_handles[chrom]
            pos = rec.pos
            ref = rec.ref
            refalt_to_rsid = {}

            for alt in rec.alts or []:
                rsid = None
                csq_list = rec.info.get("CSQ")
                if csq_list:
                    # A position may carry several CSQ entries; take the first that
                    # yields a well-formed "rs<digits>" identifier ("&"-joined list).
                    for csq_entry in csq_list:
                        fields = csq_entry.split("|")
                        if len(fields) <= rsid_index:
                            continue
                        candidate_field = fields[rsid_index]
                        for candidate in candidate_field.split("&"):
                            if candidate.startswith("rs") and candidate[2:].isdigit():
                                rsid = candidate
                                break
                        if rsid:
                            break

                if rsid is None:
                    continue

                refalt = f"{ref}/{alt}"
                rsid_int = int(rsid.replace("rs", ""))  # store numeric part only
                refalt_to_rsid[refalt] = rsid_int

            if refalt_to_rsid:
                key_bytes = struct.pack(">I", pos)  # 4-byte big-endian int
                value_bytes = msgpack.packb(refalt_to_rsid, use_bin_type=True)
                txn.put(key_bytes, value_bytes, db=db)

    env.sync()
    env.close()


def setup_rsid_mapping_lmdb(vcf_path, out_prefix, num_chroms=28):
    """
    Create LMDB mapping database if it does not already exist.
    """
    build_snp_map_lmdb_from_vcf(vcf_path, out_prefix, num_chroms=num_chroms)

def main():
    """Parse CLI args and build the rsID LMDB."""
    parser = argparse.ArgumentParser(
        description="Build rsID LMDB mapping from annotated VCF"
    )
    parser.add_argument("--vcf", required=True, help="Path to annotated VCF (bgzipped + indexed)")
    parser.add_argument("--out-file", required=True, help="Output directory for LMDB")
    parser.add_argument("--num-chroms", type=int, default=28, help="Max number of chromosome DBs")
    args = parser.parse_args()

    setup_rsid_mapping_lmdb(args.vcf, args.out_file, num_chroms=args.num_chroms)

if __name__ == "__main__":
    main()