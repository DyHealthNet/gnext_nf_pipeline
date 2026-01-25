#!/usr/bin/env python3
import os
import argparse
import sys
import numpy as np
import pysam
from math import ceil
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import shared_memory

# ------------------------------------------------------------------------------
# Globals
# ------------------------------------------------------------------------------
COLS = {
    "neg_log_pvalue": 5,
    "beta": 7,
    "stderr_beta": 8,
    "alt_allele_freq": 9,
}
METRICS = ["neg_log_pvalue", "beta", "stderr_beta", "alt_allele_freq"]

# ------------------------------------------------------------------------------
# Trait loading
# ------------------------------------------------------------------------------
def load_traits_and_paths(norm_files):
    """
    Derive trait names directly from normalized GWAS file basenames.
    E.g. '12345.gz' -> trait '12345'
    """
    traits = [os.path.basename(f).replace(".gz", "") for f in norm_files]
    return traits, norm_files

# ------------------------------------------------------------------------------
# Variant index from annotated VCF
# ------------------------------------------------------------------------------
def build_chr_variant_index(chrom: str, vcf: str):
    rows = []
    with pysam.VariantFile(vcf) as vcf_in:
        for rec in vcf_in.fetch(str(chrom)):
            c = str(rec.contig)
            p = int(rec.pos)
            r = rec.ref
            a = rec.alts[0]
            vid = f"{c}:{p}:{r}:{a}".encode()
            rows.append((c, p, r, a, vid))
    row_index = {vid: i for i, (_, _, _, _, vid) in enumerate(rows)}
    
    print(f"✓ Loaded {len(rows)} variants from VCF for chr {chrom}", file=sys.stderr)
    return rows, row_index

# ------------------------------------------------------------------------------
# Writers
# ------------------------------------------------------------------------------
def open_writers(chrom, trait_names):
    header = "#" + "\t".join(["chrom", "pos", "ref", "alt", "vid"] + trait_names) + "\n"
    paths = {
        "neg_log_pvalue": f"chr_{chrom}_neg_log_pvalue.tsv.bgz",
        "beta": f"chr_{chrom}_beta.tsv.bgz",
        "stderr_beta": f"chr_{chrom}_stderr_beta.tsv.bgz",
        "alt_allele_freq": f"chr_{chrom}_alt_allele_freq.tsv.bgz",
    }
    handles = {m: pysam.BGZFile(paths[m], "w") for m in paths}
    for h in handles.values():
        h.write(header.encode())
    return paths, handles

def close_and_index(paths, handles):
    for m, h in handles.items():
        h.close()
        pysam.tabix_index(
            paths[m],
            seq_col=0,
            start_col=1,
            end_col=1,
            meta_char="#",
            force=True,
        )

# ------------------------------------------------------------------------------
# Worker for trait filling
# ------------------------------------------------------------------------------
def _worker_fill_traits(chrom, trait_idxs, trait_paths, row_index, i0, i1,
                        shm_names, shape, dtype, cols):
    import numpy as _np
    import pysam as _pysam
    import sys as _sys
    col_p, col_b, col_se, col_af = cols

    p_shm = shared_memory.SharedMemory(name=shm_names['p'])
    b_shm = shared_memory.SharedMemory(name=shm_names['b'])
    se_shm = shared_memory.SharedMemory(name=shm_names['se'])
    af_shm = shared_memory.SharedMemory(name=shm_names['af'])

    p_arr = _np.ndarray(shape, dtype=dtype, buffer=p_shm.buf)
    b_arr = _np.ndarray(shape, dtype=dtype, buffer=b_shm.buf)
    se_arr = _np.ndarray(shape, dtype=dtype, buffer=se_shm.buf)
    af_arr = _np.ndarray(shape, dtype=dtype, buffer=af_shm.buf)
    
    for j in trait_idxs:
        fn = trait_paths[j]
        if not (os.path.exists(fn) and os.path.exists(fn + ".tbi")):
            continue
        
        matched = 0
        try:
            with _pysam.TabixFile(fn) as tbx:
                try:
                    it = tbx.fetch(str(chrom))
                except ValueError:
                    continue
                for line in it:
                    f = line.rstrip("\n").split("\t")
                    if len(f) <= max(col_p, col_b, col_se, col_af, 4):
                        continue
                    vid = f"{f[0]}:{f[1]}:{f[3]}:{f[4]}".encode()
                    
                    gi = row_index.get(vid)
                    if gi is None or gi < i0 or gi >= i1:
                        continue
                    
                    matched += 1
                    li = gi - i0
                    try: p_arr[li, j] = float(f[col_p])
                    except Exception: p_arr[li, j] = _np.nan
                    try: b_arr[li, j] = float(f[col_b])
                    except Exception: b_arr[li, j] = _np.nan
                    try: se_arr[li, j] = float(f[col_se])
                    except Exception: se_arr[li, j] = _np.nan
                    try: af_arr[li, j] = float(f[col_af])
                    except Exception: af_arr[li, j] = _np.nan
            
            print(f"  ✓ Trait {os.path.basename(fn)}: {matched} variants matched", file=_sys.stderr)
        except Exception:
            continue

    p_shm.close(); b_shm.close(); se_shm.close(); af_shm.close()

# ------------------------------------------------------------------------------
# Fill block in parallel
# ------------------------------------------------------------------------------
def fill_block_parallel(chrom, i0, i1, trait_paths, row_index,
                        p_block, b_block, se_block, af_block,
                        max_workers):
    ntraits = p_block.shape[1]
    block_rows = i1 - i0
    dtype = np.float32

    def _mk_shm(arr):
        shm = shared_memory.SharedMemory(create=True, size=arr.nbytes)
        np.ndarray(arr.shape, dtype=arr.dtype, buffer=shm.buf)[:] = arr
        return shm

    p_shm = _mk_shm(p_block)
    b_shm = _mk_shm(b_block)
    se_shm = _mk_shm(se_block)
    af_shm = _mk_shm(af_block)
    shm_names = {'p': p_shm.name, 'b': b_shm.name, 'se': se_shm.name, 'af': af_shm.name}

    shards = []
    if max_workers <= 1:
        shards = [list(range(ntraits))]
    else:
        step = ceil(ntraits / max_workers)
        for s in range(max_workers):
            j0 = s * step
            j1 = min((s + 1) * step, ntraits)
            if j0 < j1:
                shards.append(list(range(j0, j1)))

    cols = (COLS["neg_log_pvalue"], COLS["beta"], COLS["stderr_beta"], COLS["alt_allele_freq"])

    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        futs = [
            ex.submit(
                _worker_fill_traits, chrom, shard, trait_paths, row_index, i0, i1,
                shm_names, (block_rows, ntraits), dtype, cols
            )
            for shard in shards
        ]
        for f in as_completed(futs):
            f.result()

    p_arr = np.ndarray((block_rows, ntraits), dtype=dtype, buffer=p_shm.buf)
    b_arr = np.ndarray((block_rows, ntraits), dtype=dtype, buffer=b_shm.buf)
    se_arr = np.ndarray((block_rows, ntraits), dtype=dtype, buffer=se_shm.buf)
    af_arr = np.ndarray((block_rows, ntraits), dtype=dtype, buffer=af_shm.buf)
    p_block[:] = p_arr; b_block[:] = b_arr; se_block[:] = se_arr; af_block[:] = af_arr

    p_shm.close(); p_shm.unlink()
    b_shm.close(); b_shm.unlink()
    se_shm.close(); se_shm.unlink()
    af_shm.close(); af_shm.unlink()

# ------------------------------------------------------------------------------
# Main builder
# ------------------------------------------------------------------------------
def build_chromosome(chrom: str, vcf: str, norm_files, row_block_size: int, trait_workers: int):
    trait_names, trait_paths = load_traits_and_paths(norm_files)
    ntraits = len(trait_names)

    rows, row_index = build_chr_variant_index(chrom, vcf)
    nvar = len(rows)
    if nvar == 0:
        return

    paths, writers = open_writers(chrom, trait_names)

    nblocks = ceil(nvar / row_block_size)
    for b in range(nblocks):
        i0 = b * row_block_size
        i1 = min((b + 1) * row_block_size, nvar)
        block_rows = i1 - i0

        p_block = np.full((block_rows, ntraits), np.nan, dtype=np.float32)
        b_block = np.full((block_rows, ntraits), np.nan, dtype=np.float32)
        se_block = np.full((block_rows, ntraits), np.nan, dtype=np.float32)
        af_block = np.full((block_rows, ntraits), np.nan, dtype=np.float32)

        fill_block_parallel(
            chrom, i0, i1, trait_paths, row_index,
            p_block, b_block, se_block, af_block,
            max_workers=max(1, trait_workers)
        )
        
        non_nan = np.sum(~np.isnan(p_block))
        print(f"✓ Block {b+1}/{nblocks} filled: {non_nan}/{p_block.size} values", file=sys.stderr)

        for li in range(block_rows):
            c, p, r, a, vid = rows[i0 + li]
            head = [c, str(p), r, a, vid.decode()]
            writers["neg_log_pvalue"].write(("\t".join(head + [("." if not np.isfinite(x) else f"{x:.6g}") for x in p_block[li, :]]) + "\n").encode())
            writers["beta"].write(("\t".join(head + [("." if not np.isfinite(x) else f"{x:.6g}") for x in b_block[li, :]]) + "\n").encode())
            writers["stderr_beta"].write(("\t".join(head + [("." if not np.isfinite(x) else f"{x:.6g}") for x in se_block[li, :]]) + "\n").encode())
            writers["alt_allele_freq"].write(("\t".join(head + [("." if not np.isfinite(x) else f"{x:.6g}") for x in af_block[li, :]]) + "\n").encode())

    close_and_index(paths, writers)
    print(f"✓ Chr {chrom} complete: {nvar} variants, {ntraits} traits written", file=sys.stderr)

# ------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Build per-chromosome BGZF+Tabix metric files for all traits")
    parser.add_argument("--chrom", required=True, help="Chromosome name (e.g., 1 or chr1)")
    parser.add_argument("--vcf", required=True, help="Path to annotated VCF file")
    parser.add_argument("--norm-files", nargs="+", help="List of normalized GWAS summary statistics files (bgzipped+tabixed)")
    parser.add_argument("--norm-files-file", help="File containing one normalized file per line")
    parser.add_argument("--row-block-size", type=int, default=64000)
    parser.add_argument("--trait-workers", type=int, default=8)
    args = parser.parse_args()
    
    # Merge both ways of passing input
    norm_files = []
    if args.norm_files_file:
        with open(args.norm_files_file) as f:
            norm_files.extend([line.strip() for line in f if line.strip()])
    if args.norm_files:
        norm_files.extend(args.norm_files)

    build_chromosome(
        chrom=args.chrom,
        vcf=args.vcf,
        norm_files=norm_files,
        row_block_size=args.row_block_size,
        trait_workers=args.trait_workers,
    )

if __name__ == "__main__":
    sys.exit(main())