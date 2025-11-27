process generate_variant_id_lmdb {
  publishDir "${params.out_dir}/lmdb_rsid", mode: 'symlink'
  cpus 2
  memory '16 GB'

  input:
  path(annotated_vcf)
  path(annotated_vcf_tbi)
  val chroms

  output:
  path "lmdb-data.mdb", emit: lmdb_data
  path "lmdb-data.mdb-lock", emit: lmdb_lock


  script:
  """
  generate_variant_id_lmdb.py \
    --vcf ${annotated_vcf} \
    --out-file lmdb-data.mdb \
    --num-chroms ${chroms.size()}
  """
  
}