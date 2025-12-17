process generate_manhattan_qq {
  cache 'lenient'
  publishDir "${params.out_dir}/manhattan_qq", mode: 'symlink'
  tag { "${gz_files.size()}" }

  
  cpus params.manhattan_qq_cpus ?: 8
  memory params.manhattan_qq_memory ?: '64GB'


  input:
  val gz_files  // list of tuples (phenocode, file)

  output:
  path "*_manhattan.json", emit: manhattan
  path "*_qq.json", emit: qq


  script:
  
  // Write manifest.tsv directly (similar to MAGMA)
  def manifestContent = gz_files.collect { phenocode, file ->
        "${file}\t${phenocode}"
    }.join("\n")

  """
  set -e
  
  # Pre-create zorp assets directory to avoid race condition
  mkdir -p "\${CONDA_PREFIX}/share/.assets/zorp" || true
  
  # Write manifest.tsv inside the task dir
  printf "%s\n" '${manifestContent}' > manifest.tsv

  export PYTHONPATH=${projectDir}
  generate_manhattan_qq.py --input-files manifest.tsv \
      --max-workers ${task.cpus} \
      --manhattan-num-unbinned ${params.manhattan_num_unbinned ?: 500} \
      --manhattan-peak-max-count ${params.manhattan_peak_max_count ?: 500} \
      --manhattan-peak-pval-threshold ${params.manhattan_peak_pval_threshold ?: 1e-6} \
      --manhattan-peak-sprawl-dist ${params.manhattan_peak_sprawl_dist ?: 200000}
  """
  
}