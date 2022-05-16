// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'
params.options = [:]
options        = initOptions(params.options)

process MAXQUANT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda (params.enable_conda ? "bioconda::maxquant=2.0.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/https://depot.galaxyproject.org/singularity/maxquant:2.0.1.0--py39hdfd78af_1"
    } else {
        container "quay.io/biocontainers/quay.io/biocontainers/maxquant:2.0.1.0--py39hdfd78af_1"
    }

    input:
       path mqparameters
       path rawfile 
       path fastafile 
         

    output:

        path "combined/txt/*"   , emit: mq_out
        path "*.version.txt"                    , emit: version


    script:
    def software = getSoftwareName(task.process)

    """ 
    sed -i "s|PLACEHOLDER|\$PWD/|g" "${mqparameters}"
    mkdir temp
    maxquant ${mqparameters}
    maxquant --version > $maxquant.version.txt
    """
}
