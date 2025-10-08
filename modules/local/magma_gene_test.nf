process run_magma_gene_test {
    publishDir "${params.out_dir}/magma_results", mode: 'symlink'

    tag { "batch of ${magma_input_files.size()}" }

    cpus params.magma_cpus ?: 16
    memory params.magma_memory ?: '32 GB'

    input:
    tuple val(magma_input_files), path(annotation_file)

    output:
    path "*.genes.out"

    script:
    // Construct a manifest safely in Groovy first
    def manifest = magma_input_files.collect { p, f, n -> "${p}\t${f}\t${n}" }.join('\n')

    """
    set -e

    # Write the manifest to a file in the task working directory
    printf "%s\n" '${manifest}' > manifest.tsv

    IFS=\$(printf '\\t')
    while read -r phenocode file nr_samples <&3; do
        echo "Running MAGMA for \$phenocode using \$file"
        bash run_magma_gene_test.sh \
            ${params.magma_reference_plink} \
            ${annotation_file} \
            "\$file" \
            ${task.cpus} \
            "\$nr_samples" \
            "\$phenocode"_magma
    done 3< manifest.tsv
    """
}