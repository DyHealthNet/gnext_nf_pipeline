process merge_batch_reference_vcfs {
    publishDir "${params.out_dir}/annotate", mode: 'symlink'

    cpus params.vcf_cpus ?: 16
    memory params.vcf_memory ?: '64 GB'

    input:
    path(batch_vcfs)

    output:
    path "full_variants.vcf.gz", emit: vcf
    path "full_variants.vcf.gz.tbi", emit: vcf_tbi


    script:
    """
    set -e
    # Create manifest.txt from staged input files
    ls *.vcf > manifest.txt
    merge_reference_vcfs.sh manifest.txt full_variants.vcf ${params.vcf_cpus}
    """
}
