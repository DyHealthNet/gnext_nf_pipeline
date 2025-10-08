nextflow.enable.dsl=2

include { normalize          } from './modules/local/normalize.nf'
include { manhattan_qq       } from './modules/local/manhattan_qq.nf'
include { generate_vcf       } from './modules/local/reference_vcf.nf'
include { generate_chrom_bgz } from './modules/local/chrom_bgz.nf'
include { generate_top_hits  } from './modules/local/top_hits.nf'
include { extract_chroms     } from './modules/local/extract_chroms.nf'
include { VCF_ANNOTATE_ENSEMBLVEP } from './subworkflows/nf-core/vcf_annotate_ensemblvep/main'
include { generate_variant_id_lmdb } from './modules/local/variant_id_lmdb.nf'
include { DOWNLOAD_CACHE_VEP } from './subworkflows/local/download_cache_vep/main'
include { generate_magma_annotation } from './modules/local/magma_annotation.nf'
include { generate_magma_data_input } from './modules/local/magma_input_gwas.nf'
include { run_magma_gene_test } from './modules/local/magma_gene_test.nf'

workflow {
    // 1. Read phenotype manifest CSV
    // Deterministic batching for stable caching
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


    // 2. Normalize per trait in batches (of size params.pheno_batch_size)
    gwas_batches = gwas_rows.collate(params.pheno_batch_size) 

    norm_results = normalize(gwas_batches)

    // Create independent channel copies using multiMap
    norm_results.gz.flatten().multiMap { item ->
        vcf: item
        manhattan: item  
        bgz: item
        magma: item
    }.set { norm_channels }

    if(params.steps.contains("gwas_exploration")){

        // 3. Generate one full variants VCF
        vcf_files = generate_vcf(norm_channels.vcf.collect())
        vcf_file = vcf_files.vcf
        vcf_tbi  = vcf_files.vcf_tbi

        // 4. Download Ensembl VEP cache if not already present
        if (!params.ensemblvep_cache) {
            ensemblvep_info = Channel.of([ [ id:"${params.ensemblvep_cache_version}_${params.ensemblvep_genome}" ], params.ensemblvep_genome, params.ensemblvep_species, params.ensemblvep_cache_version ])
            DOWNLOAD_CACHE_VEP(ensemblvep_info)
            cache_ch = DOWNLOAD_CACHE_VEP.out.ensemblvep_cache.map{ meta, cache -> [ cache ] }
        } else {
            cache_ch = Channel.value(params.ensemblvep_cache)
        }

        // 4. Annotate full variants VCF using Ensembl VEP
        vep_results = VCF_ANNOTATE_ENSEMBLVEP(
                vcf_file.map { vcf -> tuple([id: "vcf"], vcf, [])}, //ch_vcf
                [[id: 'null'], []], //ch_fasta
                params.ensemblvep_genome,
                params.ensemblvep_species,
                params.ensemblvep_cache_version,
                cache_ch,
                []
        )
        vep_vcf = vep_results.vcf_tbi.map { it[1] }
        vep_tbi = vep_results.vcf_tbi.map { it[2] }

        // 5. Extract unique chromosomes from the VCF
        chroms = extract_chroms(vcf_file).splitText().map{it.trim()}

        // 6. Generate LMDB database for variant ID mapping
        lmdb_path = generate_variant_id_lmdb(vep_vcf, vep_tbi, chroms.splitText().map{it.trim()}.collect())

        // 7. Generate Manhattan & QQ per trait not in batches because of lmdb
        norm_batches = norm_channels.manhattan.collate(params.pheno_batch_size)
        manhattan_qq_inputs = norm_batches.combine(lmdb_path.lmdb_data).combine(lmdb_path.lmdb_lock)
                               .map { items ->
                                   def lmdb_data = items[-2]
                                   def lmdb_lock = items[-1]
                                   def batch = items[0..-3]
                                   tuple(batch, lmdb_data, lmdb_lock)  // Reconstruct the tuple
                               }

        manhattan_qq_results = manhattan_qq(manhattan_qq_inputs)
   
        // 8. Generate per-chromosome GWAS BGZ files (22 chromosome calls)
        norm_gz_files_file = norm_channels.bgz.map{it.toString()}.collectFile(
            name: 'norm_gwas_files.txt',
            storeDir: "${params.out_dir}/meta",
            newLine: true
        )
        bgz_inputs = chroms.combine(vcf_file).combine(vcf_tbi).combine(norm_gz_files_file)

        bgz_results = generate_chrom_bgz(bgz_inputs)

        // 9. Generate global top hits JSON
        manhattan_files = manhattan_qq_results.manhattan.collect()

        top_hits = generate_top_hits(manhattan_files)
    }

    if(params.steps.contains("gene_statistics")){
        // 10. Generate MAGMA annotation file
        magma_annotation = generate_magma_annotation()

        // 11. Prepare MAGMA input files per batch of phenotypes
        magma_data_inputs = norm_channels.magma.collate(params.pheno_batch_size)
        magma_input_results = generate_magma_data_input(magma_data_inputs)

        // 12. MAGMA execution
        // Extract phenocode from MAGMA input files
        magma_input_with_ids = magma_input_results.flatten().map { file_obj -> 
            def phenocode = file_obj.baseName.replaceFirst(/_magma$/, '')
            tuple(phenocode, file_obj)
        }
        // Join with GWAS metadata to get sample sizes
        meta_info = gwas_rows.map { pheno, gwas_file, n -> tuple(pheno, n)}
        magma_joined = magma_input_with_ids.join(meta_info)

        // Make batches of MAGMA input files
        magma_batches = magma_joined.collate(params.pheno_batch_size)

        // Combine with annotation file, then reconstruct tuple shape
        magma_input_final = magma_batches
            .combine(magma_annotation)
            .map { items ->
                def annotation_file = items[-1]
                def batch = items[0..-2]
                tuple(batch, annotation_file)
            }  

        // Run MAGMA gene test
        magma_final_results = run_magma_gene_test(magma_input_final)

    }
}