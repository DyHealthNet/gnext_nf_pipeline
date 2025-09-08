// Default parameter input
params.pheno_file = ""
params.gwas_dir = ""
params.out_dir = ""

// Parameters calculated from input
params.gwas_norm_dir = "${params.out_dir}/GWAS_stats_norm"


workflow {
    gwas_files = Channel
        .fromPath(params.pheno_file)
        .splitCsv(header:true)   // parse CSV with header

    normalize(gwas_files)
}

process normalize {
    conda 'envs/gwas_norm.yaml'
    tag "$row.phenocode"

    input:
    val row

    output:
    path "norm/${row.phenocode}.gz"
    path "norm${row.phenocode}.gz.tbi"
    path "norm/logs/${row.phenocode}.log"

    script:
    """
    normalize_one.py 
      --filename ${row.filename} \
      --phenocode ${row.phenocode} \
      --input-dir ${params.gwas_dir} \
      --chr-col ${params.chr_column} \
      --pos-col ${params.pos_column} \
      --ref-col ${params.ref_column} \
      --alt-col ${params.alt_column} \
      --pval-col ${params.pval_column} \
      --pval-neglog10 ${params.pval_neglog10} \
      --beta-col ${params.beta_column} \
      --se-col ${params.se_column} \
      --af-col ${params.af_column}
    """
}