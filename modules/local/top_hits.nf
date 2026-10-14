// Aggregate the global top-hits table across all phenotypes.
// Runs generate_top_hits.py over every Manhattan JSON, keeping the strongest
// variant per peak (below the cutoff) and annotating nearest genes from the LMDB.
// Inputs : file listing Manhattan JSONs, the pheno file, and the variant->gene LMDB.
// Outputs: top_hits.json.
process generate_top_hits {
  cache 'lenient'
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