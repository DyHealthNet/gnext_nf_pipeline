process run_magma_gene_test {
    publishDir "${params.out_dir}/magma_results", mode: 'symlink'

    tag { "${magma_input_files.size()}" }

    cpus params.magma_cpus ?: 16
    memory params.magma_memory ?: '32 GB'

    input:
    tuple val(magma_input_files), path(annotation_file), path(magma_reference_plink_bim), path(magma_reference_plink_bed), path(magma_reference_plink_fam), path(magma_gene_location)
    
    output:
        path "*.genes.out", emit: magma_results
        path "*.log", emit: magma_logs

    script:
    // Construct a manifest safely in Groovy first
    def manifest = magma_input_files.collect { p, f, n -> "${p}\t${f}\t${n}" }.join('\n')

    // Reference file
    def magma_reference_plink_prefix = magma_reference_plink_bim.toString().replaceFirst(/\.bim$/, '')
    """
    set -e

    # Write the manifest to a file in the task working directory
    printf "%s\n" '${manifest}' > manifest.tsv

    IFS=\$(printf '\\t')
    while read -r phenocode file nr_samples <&3; do
        echo "Running MAGMA for \$phenocode using \$file"
        bash run_magma_gene_test.sh \
            ${magma_reference_plink_prefix} \
            ${annotation_file} \
            "\$file" \
            ${task.cpus} \
            "\$nr_samples" \
            "\$phenocode"_magma
        parse_magma_gene_output.py --magma-output "\$phenocode"_magma.genes.out --gene-location ${magma_gene_location}
    done 3< manifest.tsv
    """
}