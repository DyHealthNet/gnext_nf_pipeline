// Collate every trait's MAGMA gene results into one indexed BGZ matrix.
// Writes a phenocode->file manifest, then runs generate_gene_magma_bgz.py to
// produce the gene x trait p-value table plus an LMDB offset index for fast lookup.
// Inputs : mapped_genes table, and the list of per-trait *_magma.genes.out files.
// Outputs: gene_magma_pvalues.tsv.bgz and the LMDB index (data + lock files).
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