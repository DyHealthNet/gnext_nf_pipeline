nextflow.enable.dsl=2

include { generate_vcf              } from '../../../modules/local/reference_vcf.nf'
include { extract_chroms            } from '../../../modules/local/extract_chroms.nf'
include { generate_variant_id_lmdb  } from '../../../modules/local/variant_id_lmdb.nf'
include { publish_vep_results       } from '../../../modules/local/publish_vep_results.nf'
include { VCF_ANNOTATE_ENSEMBLVEP   } from '../../nf-core/vcf_annotate_ensemblvep/main'
include { DOWNLOAD_CACHE_VEP        } from '../download_cache_vep/main'
include { generate_variant_gene_lmdb } from '../../../modules/local/variant_gene_lmdb.nf'

workflow ANNOTATE_VARIANTS {
    take:
    norm_gz_files

    main:
    // Generate reference VCF file including all unique variants
    vcf_files = generate_vcf(norm_gz_files.collect())
    vcf_file = vcf_files.vcf
    vcf_tbi  = vcf_files.vcf_tbi

    // Download Ensembl VEP cache if not already present
    if (!params.ensemblvep_cache) {
        ensemblvep_info = Channel.of([ [ id:"${params.ensemblvep_cache_version}_${params.ensemblvep_genome}" ], params.ensemblvep_genome, params.ensemblvep_species, params.ensemblvep_cache_version ])
        DOWNLOAD_CACHE_VEP(ensemblvep_info)
        cache_ch = DOWNLOAD_CACHE_VEP.out.ensemblvep_cache.map{ meta, cache -> [ cache ] }
    } else {
        cache_ch = Channel.value(params.ensemblvep_cache)
    }

    // Annotate full variants VCF using Ensembl VEP
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
    tmp = publish_vep_results(vep_vcf, vep_tbi)

    // Extract unique chromosomes from the VCF
    chroms = extract_chroms(vcf_file).splitText().map{it.trim()}

    // Generate LMDB database for variant ID mapping
    lmdb_path = generate_variant_id_lmdb(vep_vcf, vep_tbi, chroms.splitText().map{it.trim()}.collect())

    // Generate LMDB database for variant to gene mapping
    lmdb_gene_path = generate_variant_gene_lmdb(vep_vcf, vep_tbi, chroms.splitText().map{it.trim()}.collect())

    emit:
    ref_vcf = vcf_file
    ref_tbi = vcf_tbi
    //anno_vcf = vep_vcf
    //anno_tbi = vep_tbi
    chroms = chroms
    lmdb_gene_file = lmdb_gene_path.lmdb_data
    //lmdb_data = lmdb_path.lmdb_data
    //lmdb_lock = lmdb_path.lmdb_lock
}