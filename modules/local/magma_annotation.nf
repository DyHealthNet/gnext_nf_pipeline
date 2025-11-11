process generate_magma_annotation {
    publishDir "${params.out_dir}/magma_input", mode: 'symlink'

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
        ${params.magma_window_up} \
        ${params.magma_window_down} \
        magma_annotation_window_${params.magma_window_up}_up_${params.magma_window_down}_down
    """

}