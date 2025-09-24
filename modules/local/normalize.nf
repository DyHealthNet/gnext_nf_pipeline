process normalize {
  publishDir "${params.out_dir}/normalize", mode: 'symlink'

  cpus params.normalize_cpus ?: 8
  memory params.normalize_memory ?: '32 GB'

  tag { "batch ${pheno_data.size()} files" }

  input:
  val(pheno_data)  // list of [phenocode, file] tuples

  output:
  path "*.gz", emit: gz
  path "*.gz.tbi", emit: tbi

  script:
  def neglog_flag = params.pval_neglog10 ? "--pval-neglog10" : ""

  // Write manifest.tsv directly
  def manifestContent = pheno_data.collect { filename, phenocode ->
        "${filename}\t${phenocode}"
    }.join("\n")

  """
  set -e

  # Write manifest.tsv inside the task dir
  printf "%s\n" '${manifestContent}' > manifest.tsv

  echo "Manifest written:"
  cat manifest.tsv
  
  normalize.py \
    --input-files manifest.tsv \
    --max-workers ${task.cpus} \
    --chr-col ${params.chr_column} \
    --pos-col ${params.pos_column} \
    --ref-col ${params.ref_column} \
    --alt-col ${params.alt_column} \
    --pval-col ${params.pval_column} \
    ${neglog_flag} \
    --beta-col ${params.beta_column} \
    --se-col ${params.se_column} \
    --af-col ${params.af_column}
  """
}