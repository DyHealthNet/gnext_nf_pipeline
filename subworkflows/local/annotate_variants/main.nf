nextflow.enable.dsl=2

include { generate_vcf              } from '../../../modules/local/reference_vcf.nf'
include { extract_chroms            } from '../../../modules/local/extract_chroms.nf'
include { generate_variant_id_lmdb  } from '../../../modules/local/variant_id_lmdb.nf'
include { publish_vep_results       } from '../../../modules/local/publish_vep_results.nf'
include { VCF_ANNOTATE_ENSEMBLVEP   } from '../../nf-core/vcf_annotate_ensemblvep/main'
include { DOWNLOAD_CACHE_VEP        } from '../download_cache_vep/main'
include { generate_variant_gene_lmdb } from '../../../modules/local/variant_gene_lmdb.nf'
include { generate_batch_reference_vcf} from '../../../modules/local/batch_reference_vcf.nf'
include { merge_batch_reference_vcfs} from '../../../modules/local/merge_reference_vcfs.nf'

workflow ANNOTATE_VARIANTS {
    take:
    norm_gz_files
    gwas_rows

    main:

     // Create channel of (phenocode, file) pairs
    norm_files_ch = norm_gz_files.map { file ->
        def phenocode = file.baseName.replaceFirst(/\.gz$/, '')
        tuple(phenocode, file)
    }

    // Join gwas_rows with norm_files by phenocode, then batch
    norm_batches = gwas_rows
    .map { pheno, gwas_file, n -> tuple(pheno, gwas_file, n) }
    .join(norm_files_ch)
    .map { pheno, gwas_file, n, norm_file -> tuple(norm_file, pheno) }
    .collate(params.pheno_batch_size)
    .toList()
    .flatMap { batches -> 
        batches.withIndex().collect { batch, idx -> 
            tuple(idx, batch.collect { it[0] })  // Extracts just the norm_files
        }
    }
    batch_results = generate_batch_reference_vcf(norm_batches)
    
    vcf_files = merge_batch_reference_vcfs(batch_results.collect())
    
    // Generate reference VCF file including all unique variants
    //vcf_files = generate_vcf(norm_gz_files.collect())
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
    mapped_genes = lmdb_gene_path.genes_metadata
    //lmdb_data = lmdb_path.lmdb_data
    //lmdb_lock = lmdb_path.lmdb_lock
}