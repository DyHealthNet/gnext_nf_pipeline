process generate_batch_reference_vcf {
    publishDir "${params.outdir}", mode: 'symlink'

    cpus params.vcf_cpus ?: 16
    memory params.vcf_memory ?: '64 GB'


    input:
    tuple val(batch_id), path(batch_list)

    output:
    tuple val(batch_id), path("batch_variants_${batch_id}.vcf")

    script:
    """
    set -e
    # Create manifest.txt from staged input files
    ls *.gz > manifest.txt
    batch_reference_vcf.sh manifest.txt batch_variants_${batch_id}.vcf ${params.vcf_cpus}
    """
}
