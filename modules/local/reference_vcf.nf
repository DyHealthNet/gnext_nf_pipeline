process generate_vcf {
  publishDir "${params.out_dir}/annotate", mode: 'symlink'

  cpus params.vcf_cpus ?: 16
  memory params.vcf_memory ?: '64 GB'


  input:
  val norm_gz_files


  output:
  path "full_variants.vcf.gz", emit: vcf
  path "full_variants.vcf.gz.tbi", emit: vcf_tbi


  script:
  """
  generate_full_variants_vcf.sh \
    full_variants.vcf \
    ${task.cpus} \
    ${norm_gz_files.join(' ')}
  bgzip -f full_variants.vcf
  tabix -p vcf full_variants.vcf.gz
  """
}