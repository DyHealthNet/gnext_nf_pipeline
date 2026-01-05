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

    // Build manifest in Groovy - safe from shell interpretation
    def manifestContent = magma_files.collect { file ->
        def phenocode = file.toString().replaceAll(/_magma\.genes\.out$/, '')
        phenocode = new File(phenocode).name  // Get basename
        "${phenocode}\t${file}"
    }.join('\n')
    
    """
    # Write manifest safely using heredoc
    cat > manifest.tsv << 'END_MANIFEST'
    ${manifestContent}
    END_MANIFEST

    echo "Manifest file created!"

    generate_gene_magma_bgz.py \
        --gene-file ${mapped_genes} \
        --manifest manifest.tsv \
        --out-bgz-file gene_magma_pvalues.tsv.bgz \
        --out-lmdb-file lmdb-data.mdb
    """
}