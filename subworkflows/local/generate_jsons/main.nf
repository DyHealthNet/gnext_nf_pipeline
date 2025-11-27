nextflow.enable.dsl=2

include {generate_manhattan_qq} from '../../../modules/local/manhattan_qq.nf'
include {generate_top_hits} from '../../../modules/local/top_hits.nf'

workflow GENERATE_JSONS {
    take:
    norm_gz_files
    lmdb_gene_file
    
    main:
    norm_batches = norm_gz_files.collate(params.pheno_batch_size)
    manhattan_qq_results = generate_manhattan_qq(norm_batches)
    all_manhattan_files = manhattan_qq_results.manhattan

    manhattan_file = all_manhattan_files.flatten().map{it.toString()}.collectFile(
        name: 'manhattan_files.txt',
        newLine: true
    )
    
    top_hits = generate_top_hits(manhattan_file, params.pheno_file, lmdb_gene_file)

    //emit:
    //manhattan = all_manhattan_files
    //top_hits = top_hits
}