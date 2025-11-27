process generate_top_hits {
  publishDir "${params.out_dir}/top_hits", mode: 'symlink'
  cpus 1
  
  tag "top_hits"

  input:
  path manhattan_files_file
  path pheno_file
  path lmdb_gene_file

  output:
  path "top_hits.json"

  script:
    """
    generate_top_hits.py \
      --manhattan-files-file ${manhattan_files_file} \
      --phenocode-file ${pheno_file} \
      --out top_hits.json \
      --pval-cutoff ${params.top_hits_pval_cutoff} \
      --max-limit ${params.top_hits_max_limit} \
      --lmdb-gene-file ${lmdb_gene_file}
    """
}