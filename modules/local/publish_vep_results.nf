// Publish the VEP-annotated VCF (and its index) to the output directory.
// A thin staging step that copies the annotated VCF into an out/ folder so it
// lands under publishDir as the final annotation artifact.
// Inputs : annotated vcf and its tbi.
// Outputs: out/*.vcf.gz and out/*.tbi.
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