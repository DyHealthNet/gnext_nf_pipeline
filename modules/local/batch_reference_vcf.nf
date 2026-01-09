process generate_batch_reference_vcf {
    cache 'lenient'

    cpus params.vcf_cpus ?: 16
    memory params.vcf_memory ?: '64 GB'

    input:
    tuple val(batch_id), path(batch_list)

    output:
    path("batch_variants_${batch_id}.vcf")

    script:
    """
    set -e
    # Create manifest.txt from staged input files
    for f in *.gz; do printf '%s\\n' "\$f"; done > manifest.txt
    batch_reference_vcf.sh manifest.txt batch_variants_${batch_id}.vcf ${params.vcf_cpus}
    """
}
