process normalize {
  cache 'lenient'
  publishDir "${params.out_dir}/normalize", mode: 'symlink'

  cpus params.normalize_cpus ?: 8
  memory params.normalize_memory ?: '32 GB'

  tag { "${pheno_data.size()}" }

  when:
    pheno_data.size() > 0

  input:
  val(pheno_data)  // list of [phenocode, file, nr_samples]

  output:
  path "*.gz", emit: gz
  path "*.gz.tbi", emit: tbi

  script:
  def neglog_flag = params.pval_neglog10 ? "--pval-neglog10" : ""
  def sample_size_flag = params.n_samples_column ? "--sample-size-col ${params.n_samples_column}" : ""

  // Write manifest.tsv directly
  def manifestContent = pheno_data.collect { phenocode, filename, nr_samples ->
        "${filename}\t${phenocode}"
    }.join("\n")

  """
  set -e

  # Use heredoc to write manifest - safer for special characters
  cat > manifest.tsv << 'END_MANIFEST'
  ${manifestContent}
  END_MANIFEST

  echo "Manifest file created!"

  
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
    --af-col ${params.af_column} \
    ${sample_size_flag} \
  
  # Get phenocodes from manifest and process each file
  while IFS=\$'\t' read -r filename phenocode; do
    if [ -f "\${phenocode}" ]; then
      echo "Processing \${phenocode} ..."
      sort -k1,1n -k2,2n "\${phenocode}" | bgzip > "\${phenocode}.gz"
      tabix -f -p vcf "\${phenocode}.gz"
      rm "\${phenocode}"
    else
      echo "Warning: Expected file \${phenocode} not found"
    fi
  done < manifest.tsv

  echo "Final output files:"
  ls -lh *.gz *.gz.tbi
  """
}