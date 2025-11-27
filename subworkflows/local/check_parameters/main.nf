nextflow.enable.dsl=2

workflow CHECK_PARAMETERS {
    main:
     // Check input parameters
    const_params = [
        'manhattan_num_unbinned'        : params.manhattan_num_unbinned,
        'manhattan_peak_max_count'      : params.manhattan_peak_max_count,
        'manhattan_peak_pval_threshold' : params.manhattan_peak_pval_threshold,
        'manhattan_peak_sprawl_dist'    : params.manhattan_peak_sprawl_dist,
        'magma_reference_plink'         : params.magma_reference_plink,
        'magma_gene_location'           : params.gene_location,
        'window_up'               : params.window_up,
        'window_down'             : params.window_down
    ]

    if (!params.extend) {
        log.info "Not in extend mode; skipping parameter check."
        return
     }
    snapshot_path = file("${params.out_dir}/meta/workflow_params.json")
    if (snapshot_path.exists() && params.extend) {
        def old_params = new groovy.json.JsonSlurper().parseText(snapshot_path.text)
        def diffs = []

        const_params.each { k, v ->
            if (old_params.containsKey(k) && old_params[k] != v)
                diffs << "Parameter '${k}' changed: was '${old_params[k]}', now '${v}'"
        }

        if (diffs) {
            log.warn "Detected changed parameters since last run:"
            diffs.each { log.warn "   - ${it}" }
            exit 1 // Abort on parameter changes
        } else {
            log.info "Parameters match previous run."
        }
    } else {
        log.error "No previous parameter snapshot found; cannot extend."
        exit 1
    }
}