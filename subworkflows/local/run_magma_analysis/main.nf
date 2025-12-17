nextflow.enable.dsl=2

include{generate_magma_annotation} from '../../../modules/local/magma_annotation.nf'
include{generate_magma_data_input} from '../../../modules/local/magma_input_gwas.nf'
include{run_magma_gene_test} from '../../../modules/local/magma_gene_test.nf'
include {generate_gene_magma_bgz} from '../../../modules/local/gene_magma_bgz.nf'

workflow MAGMA_ANALYSIS {
    take:
        norm_gz_files   // Channel with normalized GWAS results
        gwas_rows       // Channel with GWAS metadata
        mapped_genes  // Channel with mapped gene information

    main:
    
    bim_file = file("${params.magma_reference_plink}.bim")
    gene_location_file = file(params.gene_location)
    magma_annotation = generate_magma_annotation(bim_file,gene_location_file)

    // Create channel of (phenocode, file) pairs
    norm_files_ch = norm_gz_files.map { file ->
        def phenocode = file.baseName.replaceFirst(/\.gz$/, '')
        tuple(phenocode, file)
    }

    // Join gwas_rows with norm_files by phenocode to ensure deterministic ordering
    norm_files_ordered = gwas_rows
        .map { pheno, gwas_file, n -> tuple(pheno, gwas_file, n) }
        .join(norm_files_ch)  // Joins on first element (phenocode)
        .map { pheno, gwas_file, n, norm_file -> norm_file }

    // Batch normalized files for MAGMA input generation
    all_norm_files = norm_files_ordered.collate(params.pheno_batch_size)

    // Batch generate MAGMA input files
    magma_input_results = generate_magma_data_input(all_norm_files, bim_file)
    
    // Map preprocessed MAGMA input files to phenocodes
    all_magma_inputs = magma_input_results.input_files.flatten().map { tsv_file ->
        def phenocode = tsv_file.baseName.replaceFirst(/_magma$/, '')
        tuple(phenocode, tsv_file)
    }

    // Use gwas_rows to drive deterministic ordering, join with MAGMA input files
    magma_batches = gwas_rows
        .map { pheno, gwas_file, n -> tuple(pheno, n) }
        .join(all_magma_inputs)  // Joins on phenocode
        .map { pheno, n, tsv_file -> tuple(pheno, tsv_file, n) }
        .collate(params.pheno_batch_size)
    
    magma_input_final = magma_batches
        .combine(magma_annotation)
        .map { items ->
            def annotation_file = items[-1]
            def batch = items[0..-2]
            def reference_plink_bim = file("${params.magma_reference_plink}.bim")
            def reference_plink_bed = file("${params.magma_reference_plink}.bed")
            def reference_plink_fam = file("${params.magma_reference_plink}.fam")
            def gene_location = file(params.gene_location)
            tuple(batch, annotation_file, reference_plink_bim, reference_plink_bed, reference_plink_fam, gene_location)
        }

    magma_final_results = run_magma_gene_test(magma_input_final)

    // Collect all MAGMA results from all batches into a single list
    all_magma_results = magma_final_results.magma_results.collect()
    
    magma_gene_bgz = generate_gene_magma_bgz(mapped_genes, all_magma_results)

    //emit:
    //magma_final_results
}