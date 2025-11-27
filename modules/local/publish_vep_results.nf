process publish_vep_results {
    publishDir "${params.out_dir}/annotate", mode: 'symlink'

    memory '32.GB'
    cpus 1

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