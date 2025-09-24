process generate_chrom_bgz {
    publishDir "${params.out_dir}/chr_bgz", mode: 'symlink'

    cpus params.chrom_bgz_cpus ?: 8
    maxForks params.chrom_bgz_maxForks ?: 23
    memory params.chrom_bgz_memory ?: '64GB'

    input:
    tuple val(chrom), path(vcf_file), path(vcf_tbi), path(norm_gz_files_file)

    output:
    tuple val(chrom), path("chr_${chrom}_neg_log_pvalue.tsv.bgz"), path("chr_${chrom}_neg_log_pvalue.tsv.bgz.tbi")
    tuple val(chrom), path("chr_${chrom}_beta.tsv.bgz"), path("chr_${chrom}_beta.tsv.bgz.tbi")
    tuple val(chrom), path("chr_${chrom}_stderr_beta.tsv.bgz"), path("chr_${chrom}_stderr_beta.tsv.bgz.tbi")
    tuple val(chrom), path("chr_${chrom}_alt_allele_freq.tsv.bgz"), path("chr_${chrom}_alt_allele_freq.tsv.bgz.tbi")

    script:
    """
    generate_chrom_bgz.py \
      --chrom ${chrom} \
      --vcf ${vcf_file} \
      --norm-files-file ${norm_gz_files_file} \
      --row-block-size 64000 \
      --trait-workers ${task.cpus}
    """
}