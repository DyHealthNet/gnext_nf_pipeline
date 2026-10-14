// Build the variant->gene LMDB used for nearest-gene annotation.
// Runs generate_variant_gene_lmdb.py over the annotated VCF, mapping each variant
// position to genes within the configured up/downstream window.
// Inputs : annotated vcf + tbi and the chromosome list (sizes the LMDB).
// Outputs: the LMDB (data + lock) and mapped_genes.tsv (emit: genes_metadata).
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