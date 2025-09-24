// Create independent channel copies using multiMap
    norm_results.gz.flatten().multiMap { item ->
        vcf: item
        manhattan: item  
        bgz: item
        magma: item
    }.set { norm_channels }

    // 3. Generate one full variants VCF
    vcf_files = generate_vcf(norm_channels.vcf.collect())
    vcf_file = vcf_files.vcf
    vcf_tbi  = vcf_files.vcf_tbi

    // 4. Annotate full variants VCF using Ensembl VEP
    vep_results = VCF_ANNOTATE_ENSEMBLVEP(
            vcf_file.map { vcf -> tuple([id: "vcf"], vcf, [])}, //ch_vcf
            [[id: 'null'], []], //ch_fasta
            params.ensemblvep_genome,
            params.ensemblvep_species,
            params.ensemblvep_cache_version,
            params.ensemblvep_cache,
            []
    )
    vep_vcf = vep_results.vcf_tbi.map { it[1] }
    vep_tbi = vep_results.vcf_tbi.map { it[2] }

    // 5. Extract unique chromosomes from the VCF
    chroms = extract_chroms(vcf_file).splitText().map{it.trim()}

    // 6. Generate LMDB database for variant ID mapping
    lmdb_path = generate_variant_id_lmdb(vep_vcf, vep_tbi, chroms.splitText().map{it.trim()}.collect())

    if(params.steps.contains("gwas_exploration")){
        // 7. Generate Manhattan & QQ per trait in batches (6 batches)
        norm_batches = norm_channels.manhattan.collate(params.pheno_batch_size)
        manhattan_qq_inputs = norm_batches.combine(lmdb_path.lmdb_data)
                               .map { items ->
                                   def lmdb = items.last()
                                   def batch = items[0..-2]
                                   tuple(batch, lmdb)  // Reconstruct the tuple
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
        //10. Generate MAGMA input files

        // Generate MAGMA data input files in batches
        norm_batches = norm_channels.magma.collate(params.pheno_batch_size)
        magma_inputs = generate_magma_data_input(norm_batches, lmdb_path.lmdb_data)
        
        // Generate MAGMA configuration files
        magma_config_files = generate_magma_mapping_input(vep_vcf, params.magma_conf_file)

        //11. Run MAGMA with the different configurations
    }