// List the distinct chromosomes present in the reference VCF.
// Used to fan out per-chromosome processing (e.g. generate_chrom_bgz).
// Inputs : the (gzipped) reference VCF.
// Outputs: chroms.txt - one unique chromosome name per line.
process extract_chroms {
    cpus 1
    memory '1 GB'
    time '30m'

    input:
    path vcf_file

    output:
    path "chroms.txt"

    script:
    """
    # Skip header, take first column, unique sort
    gunzip -c ${vcf_file} | awk '!/^#/ {print \$1}' | sort -u > chroms.txt
    """
}