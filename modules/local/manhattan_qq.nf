process generate_manhattan_qq {
  publishDir "${params.out_dir}/manhattan_qq", mode: 'symlink'
  tag { "${gz_files.size()}" }

  
  cpus params.manhattan_qq_cpus ?: 8
  memory params.manhattan_qq_memory ?: '64GB'


  input:
  val gz_files

  output:
  path "*_manhattan.json", emit: manhattan
  path "*_qq.json", emit: qq


  script:

  """
  # Clean up the Nextflow stringified list "[a, b, c]" -> "a b c"
  files=\$(echo $gz_files | tr -d '[],')
  
  # Write manifest.tsv with <file>\t<basename>
  > manifest.tsv
  for f in \$files; do
    echo -e "\$f\t\$(basename \$f .gz)" >> manifest.tsv
  done

  export PYTHONPATH=${projectDir}
  generate_manhattan_qq.py --input-files manifest.tsv \
      --max-workers ${task.cpus} \
      --manhattan-num-unbinned ${params.manhattan_num_unbinned ?: 500} \
      --manhattan-peak-max-count ${params.manhattan_peak_max_count ?: 500} \
      --manhattan-peak-pval-threshold ${params.manhattan_peak_pval_threshold ?: 1e-6} \
      --manhattan-peak-sprawl-dist ${params.manhattan_peak_sprawl_dist ?: 200000}
  """
  
}