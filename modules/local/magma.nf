process generate_magma_data_input {
    publishDir "${params.out_dir}/magma_input", mode: 'symlink'

    tag { "batch ${gz_files.size()} files" }

    input:
    val gz_files
    path lmdb_path


    output:
    path "*.txt"

    script:
    def cmds = gz_files.flatten().collect { gz ->
        def base = gz.baseName
        "export PYTHONPATH=${projectDir} && generate_magma_data_input.py " +
        "--input ${gz.toString()} " +
        "--phenocode ${base} " +
        "--lmdb ${lmdb_path}"
    }
    
    """
    set -e
    echo "===== COMMANDS ====="
    echo "${cmds.join('\n')}"
    echo "===================="
    parallel -j ${task.cpus} ::: ${cmds.collect{"'${it}'"}.join(" ")}
    """
}

process generate_magma_mapping_input {
    publishDir "${params.out_dir}/magma_input", mode: 'symlink'

    input:
    path annotated_vcf_file
    val magma_config_file

    output:
    path "*_annotations.genes.annot"

    script:
    """
    generate_magma_mapping_input.py \\
        --input ${annotated_vcf_file} \\
        --config ${maga_config_file}") \\
    """

}