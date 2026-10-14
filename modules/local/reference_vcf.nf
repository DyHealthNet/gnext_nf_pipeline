// Build the reference VCF directly from all normalized files in one shot.
// Runs generate_full_variants_vcf.sh over every normalized .gz, then bgzips and
// tabix-indexes the result. (Single-step alternative to the batch + merge path.)
// Inputs : list of all normalized .gz files.
// Outputs: full_variants.vcf.gz (emit: vcf) and its .tbi (emit: vcf_tbi).
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