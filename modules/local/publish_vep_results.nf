process publish_vep_results {
    publishDir "${params.out_dir}/annotate", mode: 'symlink'

    input:
    path vcf 
    path tbi

    output:
    path "out/*.vcf.gz"
    path "out/*.tbi"

    script:
    """
    mkdir -p out
    cp "$vcf" "out/"
    cp "$tbi" "out/"
    """
}