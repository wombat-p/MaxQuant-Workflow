// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'


params.options = [:]
options        = initOptions(params.options)

process SDRFPIPELINES {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }


    conda (params.enable_conda ? "bioconda::sdrf-pipelines=0.0.14--py_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/sdrf-pipelines:0.0.14"
    } else {
        container "quay.io/biocontainers/sdrf-pipelines:0.0.14--py_0"
    }

    input:
    path sdrf
    path fasta
    

    output:
    path "mqpar.xml"         , emit: maxquantpar
    path "exp_design.tsv"         , emit: tsv
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
   parse_sdrf \\
    convert-maxquant \\
    -s "${sdrf}" \\
    -f "PLACEHOLDER${fasta}" \\
#    -m ${params.match} \\
    -r PLACEHOLDER \\
    $options.args \\
#    -pef ${params.peptidefdr} \\
#    -prf ${params.proteinfdr} \\
    -t PLACEHOLDERtemp \\
    -o2 exp_design.tsv \\
    -n ${task.cpus} 
     
    echo $VERSION > ${software}.version.txt

    """
}
