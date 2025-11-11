nextflow.enable.dsl=2

include { generate_chrom_bgz } from '../../../modules/local/chrom_bgz.nf'

workflow GENERATE_BGZ_FILES {

    take:
    norm_gz_files
    chroms
    vcf_file
    vcf_tbi

    main:
    //Generate per-chromosome GWAS BGZ files (22 chromosome calls)
    norm_gz_files_file = norm_gz_files.map{it.toString()}.collectFile(
        name: 'norm_gwas_files.txt',
        newLine: true
    )
    bgz_inputs = chroms.combine(vcf_file).combine(vcf_tbi).combine(norm_gz_files_file)

    bgz_results = generate_chrom_bgz(bgz_inputs)

    //emit:
    //bgz_results
}
