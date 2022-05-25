process CALCBENCHMARKS {
    label 'process_low'
    label 'process_single_thread'
    publishDir "${params.outdir}/", mode:'copy'
    conda (params.enable_conda ? "conda-forge::notyetavailable" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/notyetavailable"
    } else {
        container "wombatp/maxquant-pipeline:dev"
    }


    input:
        val foo
        path exp_design_file
        path std_prot_file
        path std_pep_file
        path fasta_file
 
  output:
   path "params.json",   emit: parameters
   path "benchmarks.json",  emit:  benchmarks
  
  script:
  """
  echo '$foo' > params.json
  cp "${fasta_file}" database.fasta
  Rscript $baseDir/scripts/CalcBenchmarks.R

  """
}


