#!/usr/bin/env nextflow

/*
 * merge_vcfs.nf: Nextflow process for merging batch VCFs
 * Usage: Provide a list of batch VCFs as input, outputs a merged VCF
 */

process merge_batch_reference_vcfs {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(batch_vcfs)

    output:
    path "full_variants.vcf.gz", emit: vcf
    path "full_variants.vcf.gz.tbi", emit: vcf_tbi


    script:
    """
    bash bin/merge_reference_vcfs.sh full_variants.vcf ${batch_vcfs.join(' ')}
    """
}
