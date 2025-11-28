process generate_gene_magma_bgz {
    publishDir "${params.out_dir}/magma_results", mode: 'symlink'

    input:
    path(mapped_genes)
    val(magma_files)

    output:
    path "gene_magma_pvalues.tsv.bgz"
    path "lmdb-data.mdb", emit: lmdb_data
    path "lmdb-data.mdb-lock", emit: lmdb_lock

    script:
    
    """
    # Clean up the Nextflow stringified list "[a, b, c]" -> "a b c"
    files=\$(echo ${magma_files} | tr -d '[],')

    # Write manifest.tsv with phenocode<tab>file format
    > manifest.tsv
    for f in \$files; do
        # Extract phenocode by removing _magma.genes.out suffix
        phenocode=\$(basename "\$f" _magma.genes.out)
        echo -e "\${phenocode}\t\${f}" >> manifest.tsv
    done

    generate_gene_magma_bgz.py \
        --gene-file ${mapped_genes} \
        --manifest manifest.tsv \
        --out-bgz-file gene_magma_pvalues.tsv.bgz \
        --out-lmdb-file lmdb-data.mdb
    """
}