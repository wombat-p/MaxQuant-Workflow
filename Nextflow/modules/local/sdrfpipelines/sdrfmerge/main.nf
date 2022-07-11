process SDRFMERGE { 
    publishDir "${params.outdir}/sdrfmerge"
    label 'process_medium'
    conda (params.enable_conda ? "bioconda::sdrf-pipelines=0.0.21" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker://wombatp/maxquant-pipeline:dev"
    } else {
        container "wombatp/maxquant-pipeline:v0.14"
    }


    input:
      path sdrf
      path parameters
      path map    

    output:
      path "sdrf_local.tsv"         , emit: sdrf_local


    script:
    """
    if [[ "$sdrf" != "sdrf.tsv" ]]
    then
	cp "${sdrf}" sdrf.tsv
    fi
    if [[ "$parameters" != "params.yml" ]] 
    then
        cp "${parameters}" params.yml
    fi
    if [[ "$map" != "params2sdrf.yml" ]]
    then
        cp "${map}" params2sdrf.yml
    fi
    # TODO change to package when available
    python $projectDir/scripts/add_data_analysis_param.py
    python $projectDir/scripts/sdrf2params.py
    """
}
