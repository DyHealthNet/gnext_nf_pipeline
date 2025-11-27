process generate_magma_data_input {
  publishDir "${params.out_dir}/magma_input", mode: 'symlink'

  cpus params.magma_input_cpus ?: 32
  memory params.magma_input_memory ?: '64 GB'

  tag { "${gz_files.size()}" }

  input:
    val(gz_files)
    path magma_reference_plink_bim

  output:
  path "*_magma.tsv", emit: input_files
  path "mapping_summary.tsv", emit: mapping_summary

  script:
  // Check if we need to include n_samples column
  def n_samples_flag = params.n_samples_column ? "--include-n-samples" : ""

  """
  set -e

  # Clean up the Nextflow stringified list "[a, b, c]" -> "a b c"
  files=\$(echo $gz_files | tr -d '[],')
  
  # Write manifest.tsv with <file>\t<basename>
  > manifest.tsv
  for f in \$files; do
    echo -e "\$f\t\$(basename \$f .gz)" >> manifest.tsv
  done

  generate_magma_data_input.py \
    --input-files manifest.tsv \
    --max-workers ${task.cpus} \
    --genome-build ${params.ensemblvep_genome} \
    --ref-bim ${magma_reference_plink_bim} \
    ${n_samples_flag}
  """
}