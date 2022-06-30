process NORMALYZERDE {
    label 'process_intermediate'
    label 'process_single_thread'
    publishDir "${params.outdir}/normalyzerde", mode:'copy'
    conda (params.enable_conda ? "bioconda:bionconductor-normalyzerde=1.14.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker://computationalproteomics/normalyzerde:1.14.0"
    } else {
        container "computationalproteomics/normalyzerde:1.14.0"
    }

input:
    path maxquant
    path exp_file
    val normalization
  
    output:
	path "exp_design_calcb.tsv" , emit: exp_design
	path "NormalyzerProteins/*"   , emit:  normalyzer_proteins
	path "NormalyzerPeptides/*"   , emit:  normalyzer_peptides
        path "stand_prot_quant_merged.csv"    , emit: std_prots
        path "stand_pep_quant_merged.csv"    , emit: std_peps

    when:
    params.run_statistics

    script:
    """
    cp "${exp_file}" Normalyzer_design.tsv
    cp "proteinGroups.txt" protein_file.txt
    cp "peptides.txt" peptide_file.txt
    Rscript $baseDir/scripts/runNormalyzer.R --comps="${params.comps}" --method="${normalization}"
    """

}
