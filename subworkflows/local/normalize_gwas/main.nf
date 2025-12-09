nextflow.enable.dsl=2

include { normalize } from '../../../modules/local/normalize.nf'

workflow NORMALIZE_GWAS {
    take:
    gwas_rows

    main:
    gwas_batches = gwas_rows.collate(params.pheno_batch_size)
    norm_results = normalize(gwas_batches)
    new_norm_gz = norm_results.gz.collect().flatten()

    // Create independent channel copies for new files only
    new_norm_gz.multiMap { item ->
        vcf: item
        json: item  
        bgz: item
        magma: item
    }.set { gz_files }

    emit:
    vcf    = gz_files.vcf
    json   = gz_files.json
    bgz    = gz_files.bgz
    magma  = gz_files.magma
}