// Annotate SNPs to genes with MAGMA (step 1 of the gene-level analysis).
// Runs run_magma_annotation.sh to map reference-panel SNPs to genes using the
// configured up/downstream window, producing the shared .genes.annot file.
// Inputs : reference PLINK .bim and the gene location file.
// Outputs: the magma_annotation_window_*.genes.annot file.
process generate_magma_annotation {
    publishDir "${params.out_dir}/magma_input", mode: 'symlink'

    cpus 1
    memory '64 GB'

    input:
    path magma_reference_plink_bim
    path magma_gene_location

    output:
    path "*.genes.annot"

    script:
    """

    bash run_magma_annotation.sh \
        ${magma_reference_plink_bim} \
        ${magma_gene_location} \
        ${params.window_up} \
        ${params.window_down} \
        magma_annotation_window_${params.window_up}_up_${params.window_down}_down
    """

}