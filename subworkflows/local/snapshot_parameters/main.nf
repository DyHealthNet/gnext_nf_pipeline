nextflow.enable.dsl=2

workflow SNAPSHOT_PARAMETERS {
    main:
    // Store workflow configurations to meta folder
    def paramsFile = "${params.out_dir}/meta/workflow_params.json"
    // Write parameters as pretty-printed JSON and create directory
    new File(paramsFile).parentFile.mkdirs()
    new File(paramsFile).withWriter {writer -> 
    writer << groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params))

    // Store processed phenotype manifest to meta folder -> copy params.pheno_file to out_dir/meta
    def phenoManifestFile = "${params.out_dir}/meta/phenotype_manifest.csv"
    new File(phenoManifestFile).withWriter {phenoWriter ->
        phenoWriter << new File(params.pheno_file).text
    }
    }   
}