process generate_variant_gene_lmdb {
  publishDir "${params.out_dir}/lmdb_gene", mode: 'symlink'
  
  cpus 2
  memory '16 GB'

  input:
  path(annotated_vcf)
  path(annotated_vcf_tbi)
  val(chroms)

  output:
  path "lmdb-data.mdb", emit: lmdb_data
  path "lmdb-data.mdb-lock", emit: lmdb_lock
  path "mapped_genes.tsv", emit: genes_metadata

  script:  
  """
  generate_variant_gene_lmdb.py \
    --vcf ${annotated_vcf} \
    --gene-file ${params.gene_location} \
    --out-file lmdb-data.mdb \
    --window-up ${params.window_up} \
    --window-down ${params.window_down} \
    --num-chroms ${chroms.size()}
  """
}