nextflow.enable.dsl=2

include { NORMALIZE_GWAS } from './subworkflows/local/normalize_gwas/main'
include { ANNOTATE_VARIANTS } from './subworkflows/local/annotate_variants/main'
include { GENERATE_JSONS } from './subworkflows/local/generate_jsons/main'
include { GENERATE_BGZ_FILES } from './subworkflows/local/generate_bgz_files/main'
include { MAGMA_ANALYSIS } from './subworkflows/local/run_magma_analysis/main'
include { CHECK_PARAMETERS } from './subworkflows/local/check_parameters/main'
include { SNAPSHOT_PARAMETERS } from './subworkflows/local/snapshot_parameters/main'

workflow {

    // Check for parameter changes since last run
    //CHECK_PARAMETERS()

    // Read phenotype manifest CSV -> deterministic batching for stable caching
   gwas_rows = Channel
        .fromPath(params.pheno_file)
        .splitCsv(header: true)
        .collect()                                // gather all rows first
        .flatMap { rows ->
            rows
                .sort { a, b -> a.phenocode <=> b.phenocode }  // sort by phenocode
                .collect { row -> 
                    def n_int = row.nr_samples.toString().replace('.0', '').toInteger()
                    tuple(row.phenocode, row.filename.toString(), n_int) }
        }

    // Normalize GWAS summary statistics
    NORMALIZE_GWAS(gwas_rows)

    // Proceed with downstream GWAS exploration and gene statistics

    if(params.steps.contains("gwas_exploration")){

        // Annotate variants using Ensembl VEP
        ANNOTATE_VARIANTS(NORMALIZE_GWAS.out.vcf)

        // Generate JSONs for Manhattan & QQ plots, and top hits
        GENERATE_JSONS(NORMALIZE_GWAS.out.json)

        // Generate per-chromosome BGZ files
        GENERATE_BGZ_FILES(
            NORMALIZE_GWAS.out.bgz,
            ANNOTATE_VARIANTS.out.chroms,
            ANNOTATE_VARIANTS.out.ref_vcf,
            ANNOTATE_VARIANTS.out.ref_tbi
        )
    }

    if(params.steps.contains("gene_statistics")){
        MAGMA_ANALYSIS(
            NORMALIZE_GWAS.out.magma,
            gwas_rows
        )
    }

    SNAPSHOT_PARAMETERS() 
}