nextflow.enable.dsl=2

include { normalize } from '../../../modules/local/normalize.nf'

// NORMALIZE_GWAS: normalize all raw GWAS files (the pipeline's first step).
// Batches the pheno manifest rows, runs the `normalize` process per batch, then
// fans the resulting .gz files out into four identical channels so the downstream
// subworkflows (VCF/annotation, JSONs, BGZ matrices, MAGMA) can each consume them
// independently.
//   take : gwas_rows (pheno manifest rows: [phenocode, file, nr_samples]).
//   emit : vcf / json / bgz / magma - the same normalized .gz files on four channels.
workflow NORMALIZE_GWAS {
    take:
    gwas_rows

    main:
    // Group rows into batches so each `normalize` task processes several files.
    gwas_batches = gwas_rows.collate(params.pheno_batch_size)
    norm_results = normalize(gwas_batches)
    // Flatten the per-batch output lists into a single stream of .gz files.
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