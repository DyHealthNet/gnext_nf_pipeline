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
    // Construct a manifest with phenocode, file, and sample size
    def manifestContent = magma_input_files.collect { p, f, n -> "${p}\t${f}\t${n}" }.join('\n')
    
    // Determine if sample sizes come from pheno file (gwas_rows) or from GWAS file column
    def use_n_column = params.n_samples_column ? "true" : "false"

    // Reference file
    def magma_reference_plink_prefix = magma_reference_plink_bim.toString().replaceFirst(/\.bim$/, '')
    """
    set -e

    # Use heredoc to write manifest - safer for special characters
    cat > manifest.tsv << 'END_MANIFEST'
${manifestContent}
END_MANIFEST

    echo "Manifest file created!"

    IFS=\$(printf '\\t')
    while read -r phenocode magma_file nr_samples <&3; do
        echo "Running MAGMA for \$phenocode using preprocessed file \$magma_file"
        
        # Determine sample size parameter for MAGMA
        if [ "${use_n_column}" = "true" ]; then
            # Sample size is in the preprocessed MAGMA file (column 3)
            n_param="ncol=3"
            echo "Using sample size from n_samples column (column 3) in preprocessed file"
        else
            # Sample size comes from pheno file (gwas_rows)
            n_param="N=\${nr_samples}"
            echo "Using sample size from pheno file: \${nr_samples}"
        fi
        
        # Run MAGMA with preprocessed file
        # Preprocessed columns: SNP=1, P=2, optionally n_samples=3
        bash run_magma_gene_test.sh \
            ${magma_reference_plink_prefix} \
            ${annotation_file} \
            "\$magma_file" \
            ${task.cpus} \
            "\$n_param" \
            "\$phenocode"_magma
        
        parse_magma_gene_output.py --magma-output "\$phenocode"_magma.genes.out --gene-location ${magma_gene_location}
    done 3< manifest.tsv
    """
}