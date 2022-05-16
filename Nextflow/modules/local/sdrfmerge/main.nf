// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'


params.options = [:]
options = initOptions(params.options)
def VERSION = '0.0.2'


process SDRFMERGE { 
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process)) }

    conda (params.enable_conda ? "bioconda::sdrf-pipelines=0.0.21" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/sdrf-pipelines:0.0.21--py_0"
    } else {
        container "wombatp/maxquant-pipeline:dev"
    }


    input:
      path sdrf
      path parameters
      path map
     

    output:
      path "sdrf_local.tsv"         , emit: sdrf_local
    // TODO nf-core: List additional required output channels/values here
      path "*.version.txt"          , emit: version


    script:
      def software = getSoftwareName(task.process)
  
    """
    if [ "$sdrf" -ne "sdrf.tsv" ]
    then
	cp "${sdrf}" sdrf.tsv
    fi
    cp "${parameters}" params.yml
    cp "${map}" params2sdrf.yml
    python $projectDir/scripts/sdrf-mapper.py
    echo "preliminary version" > sdrf-merge.version.txt
    """
}
