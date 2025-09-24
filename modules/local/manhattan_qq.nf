process manhattan_qq {
  debug true
  
  cpus params.manhattan_qq_cpus ?: 8
  memory params.manhattan_qq_memory ?: '64GB'

  publishDir "${params.out_dir}/manhattan_qq", mode: 'symlink'

  input:
  tuple val(gz_files), path(lmdb_data_path), path(lmdb_lock_path)

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

  mkdir lmdb-database
  mv lmdb-data.mdb lmdb-database/data.mdb
  mv lmdb-data.mdb-lock lmdb-database/lock.mdb


  export PYTHONPATH=${projectDir}
  generate_manhattan_qq.py --input-files manifest.tsv \
      --lmdb lmdb-database/data.mdb \
      --max-workers ${task.cpus} \
      --manhattan-num-unbinned ${params.manhattan_num_unbinned ?: 500} \
      --within-pheno-mask-around-peak ${params.within_pheno_mask_around_peak ?: 500000} \
      --between-pheno-mask-around-peak ${params.between_pheno_mask_around_peak ?: 1000000} \
      --manhattan-peak-max-count ${params.manhattan_peak_max_count ?: 500} \
      --manhattan-peak-pval-threshold ${params.manhattan_peak_pval_threshold ?: 1e-6} \
      --manhattan-peak-sprawl-dist ${params.manhattan_peak_sprawl_dist ?: 200000} \
      --manhattan-peak-variant-counting-pval-threshold ${params.manhattan_peak_variant_counting_pval_threshold ?: 5e-8}
  """
  
}