nextflow.enable.dsl=2

include { generate_chrom_bgz } from '../../../modules/local/chrom_bgz.nf'

// GENERATE_BGZ_FILES: build the per-chromosome BGZ metric matrices.
// Collects the normalized file paths into a single list file, fans out over
// chromosomes (combining each with the reference VCF + index + file list), and
// runs generate_chrom_bgz per chromosome.
//   take : norm_gz_files, chroms, vcf_file, vcf_tbi.
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