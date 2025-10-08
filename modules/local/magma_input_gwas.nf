process generate_magma_data_input {
  publishDir "${params.out_dir}/magma_input", mode: 'symlink'

  cpus params.magma_input_cpus ?: 32
  memory params.magma_input_memory ?: '64 GB'

  tag { "batch ${gz_files.size()} files" }

  input:
    val(gz_files)

  output:
  path "*_magma.tsv", emit: gz

  script:  

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
    --ref-bim ${params.magma_reference_plink + ".bim"} \
  """
}