nextflow.enable.dsl=2

include{generate_magma_annotation} from '../../../modules/local/magma_annotation.nf'
include{generate_magma_data_input} from '../../../modules/local/magma_input_gwas.nf'
include{run_magma_gene_test} from '../../../modules/local/magma_gene_test.nf'

workflow MAGMA_ANALYSIS {
    take:
        norm_gz_files   // Channel with normalized GWAS results
        gwas_rows       // Channel with GWAS metadata

    main:
    
    bim_file = file("${params.magma_reference_plink}.bim")
    gene_location_file = file(params.magma_gene_location)
    magma_annotation = generate_magma_annotation(bim_file,gene_location_file)

    magma_data_inputs = norm_gz_files.collate(params.pheno_batch_size)
    magma_input_results = generate_magma_data_input(magma_data_inputs, bim_file)
    
    all_magma_inputs = magma_input_results.input_files.flatten().map { file_obj -> 
        def phenocode = file_obj.baseName.replaceFirst(/_magma$/, '')
        tuple(phenocode, file_obj)
    }

    // Join with GWAS metadata to get sample sizes
    meta_info = gwas_rows.map { pheno, gwas_file, n -> tuple(pheno, n)}
    magma_with_meta = all_magma_inputs.join(meta_info)

    
    magma_batches = magma_with_meta.collate(params.pheno_batch_size)
    
    magma_input_final = magma_batches
        .combine(magma_annotation)
        .map { items ->
            def annotation_file = items[-1]
            def batch = items[0..-2]
            def reference_plink_bim = file("${params.magma_reference_plink}.bim")
            def reference_plink_bed = file("${params.magma_reference_plink}.bed")
            def reference_plink_fam = file("${params.magma_reference_plink}.fam")
            def gene_location = file(params.magma_gene_location)
            tuple(batch, annotation_file, reference_plink_bim, reference_plink_bed, reference_plink_fam, gene_location)
        }

    magma_final_results = run_magma_gene_test(magma_input_final)

    //emit:
    //magma_final_results
}