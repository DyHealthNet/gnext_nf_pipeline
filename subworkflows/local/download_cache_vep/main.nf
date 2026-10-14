
//
// DOWNLOAD CACHE SNPEFF VEP
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Condition is based on params.step and params.tools
// If and extra condition exists, it's specified in comments

include { ENSEMBLVEP_DOWNLOAD } from '../../../modules/nf-core/ensemblvep/download/main'

// DOWNLOAD_CACHE_VEP: fetch the Ensembl VEP cache when none is provided.
// Wraps the nf-core ENSEMBLVEP_DOWNLOAD module; ANNOTATE_VARIANTS calls this only
// when params.ensemblvep_cache is unset.
//   take : ensemblvep_info (meta, genome, species, cache version).
//   emit : ensemblvep_cache ([meta, cache]) and tool versions.
workflow DOWNLOAD_CACHE_VEP {
    take:
    ensemblvep_info

    main:
    versions = Channel.empty()

    ENSEMBLVEP_DOWNLOAD(ensemblvep_info)

    // Gather versions of all tools used
    versions = versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)

    emit:
    ensemblvep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.collect()  // channel: [ meta, cache ]

    versions // channel: [ versions.yml ]
}
