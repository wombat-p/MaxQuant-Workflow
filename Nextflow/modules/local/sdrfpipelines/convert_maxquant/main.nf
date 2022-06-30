process CONVERT_MAXQUANT {
    publishDir "${params.outdir}/sdrfmerge"
    label 'process_medium'
    conda (params.enable_conda ? "bioconda::sdrf-pipelines=0.0.21--py_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker://wombatp/maxquant-pipeline:dev"
    } else {
//        container "quay.io/biocontainers/sdrf-pipelines:0.0.21--py_0"
        container "wombatp/maxquant-pipeline:dev"
    }

    input:
    path sdrf
    path fasta
    

    output:
    path "mqpar.xml"         , emit: maxquantpar
    path "exp_design.tsv"         , emit: exp_design
    path "*.version.txt"          , emit: version

    script:
    """
    parse_sdrf \\
    convert-maxquant \\
    -s "${sdrf}" \\
    -f "PLACEHOLDER${fasta}" \\
    -r PLACEHOLDER \\
    -t PLACEHOLDERtemp \\
    -o2 exp_design.tsv \\
    -n ${task.cpus} 
    echo "Preliminary" > sdrf_merge.version.txt

    parse_sdrf \\
    convert-normalyzerde \\
    -s "${sdrf}" \\
    -o exp_design.tsv
    """
}
