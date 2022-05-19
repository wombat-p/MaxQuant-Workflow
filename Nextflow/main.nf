#!/usr/bin/env nextflow
/*
========================================================================================
                         MaxQuant-Pipeline
========================================================================================
(NOT YET A nf-core!)
 #### Homepage / Documentation
----------------------------------------------------------------------------------------
*/

def helpMessage() {
log.info nfcoreHeader()
log.info"""
Usage: 

The typical command for running the pipeline is as follwos:
nextflow run main.nf --raws "path/to/raws" --sdrf "path/to/*.tsv" --fasta "path/to/*.fasta" -profile docker

Mandatory arguments: 
--raws      raw files of input data (must be surrounded with quotes)
--sdrf      Path to sdrf.tsv file (must be surrounded with quotes, and in tsv format)
--fasta     Path tp the fasta file(must be surrounded with quotes, and in .fasta)



-profile    Configuration profile to use. Can use multiple (comma separated)
            Available: docker


Other options:
    --peptidefdr    posterior error probability calculation based on target-decoy search, default=0.01 
    --proteinfdr    protein score = product of peptide PEPs (one for each sequence)', default=0.01
    --match         "True" via matching between runs to boosts number of identifications, default = True
    --outdir        The output directory where the results will be saved
    --email         Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
    --experiment_design               Text-file containing 2 columns: first with raw file names and second with names for 
experimental conditions
    --run_statistics                  Either true or false
    --normalyzerMethod Method for normalization of samples: can be log2, CycLoess, median, mean, Quantile, GI, RLR or VSN
    --comparisons   Comparisons for statistical tests. Needs to be in the format "Cond2-Cond1,Cond3-Cond1" assuming that CondN are the conditions in the experimental design file. 
Run statistics


    -name           Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.


""".stripIndent()
}

/*
* Read and setup the variables
*/ 
params.raws = params.raws ?: { log.error "No raw data provided. Make sure you have used the '--raws' option."; exit 1 }()
params.sdrf = params.sdrf ?: {log.error "No SDRF file provided. Make sure you have used the '--sdrf' option."; exit 1}()
params.fasta = params.fasta ?: { log.error "No fasta file provided. Make sure you have used the '--fasta' option."; exit 1}()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()
 
/*
* Define the default paramaters. 
* The parameters are found in the xml file and is change within the file.
*/

// SDRF parameters
params.peptidefdr = 0.01
params.proteinfdr = 0.01

params.match = "True"

// NormalyzerDE parameters
// Which groups to compare
params.run_statistics = true
params.comparisons = ''
params.normalyzerMethod = "log2"

/*
 * SET UP CONFIGURATION VARIABLES
 */


// Configurable variables
params.name = false
params.email = false
params.plaintext_email = false


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
    custom_runName = workflow.runName
}

/*
* Generate the channels for the raw files
*/ 
input_raw_maxquant = Channel.fromPath(params.raws)

/*
* Generate the channels for the sdrf files
*/ 
input_sdrf_parse = Channel.fromPath (params.sdrf)

/*
* Generate the channels for the fasta file
*/ 
input_fasta = Channel.fromPath (params.fasta)
input_fasta.into { input_fasta_parse; input_fasta_maxquant }

/*
* Generate channel for experimental design file
*/
input_exp_design = Channel.fromPath(params.experiment_design)


/* 
* STEP 1 - Generate the parameter and experimental design for maxquant through sdrf
*/
process run_sdrf {
    label 'process_high'
    // this is run as high to pass the (maximal) number of threads (which will be 4 as default without label)

    publishDir "${params.outdir}/params", mode: 'copy'
    
    input: 
    file sdrf_file from input_sdrf_parse
    file fasta_file from input_fasta_parse
    
    output:
    file "mqpar.xml" into mq_parameters
    file "exp_design.tsv" into exp_design_parsed

    script: 
    """
    parse_sdrf convert-maxquant -s "${sdrf_file}" -f "PLACEHOLDER${fasta_file}" -m ${params.match} -r PLACEHOLDER -pef ${params.peptidefdr} -prf ${params.proteinfdr} -t PLACEHOLDERtemp -o2 exp_design.tsv -n ${task.cpus} 
    """
}

/*
* STEP 2 - Run maxQuant with parameters from mqpar
*/
process run_maxquant {
    label 'process_high'

    publishDir "${params.outdir}/maxq", mode: 'copy'
    
    input:
    file rawfile from input_raw_maxquant.collect()
    file fastafile from input_fasta_maxquant
    file mqparameters from mq_parameters

    output:
    file "combined/txt/proteinGroups.txt" into protein_groups
    file "combined/txt/peptides.txt" into peptides
    file "${mqparameters}" into mq_par_final
    
    script:
    """
    sed -i "s|PLACEHOLDER|\$PWD/|g" "${mqparameters}"
    mkdir temp
    maxquant ${mqparameters}
    cp -R "\$PWD/combined/txt" "${params.outdir}"
    """
}

/*
* STEP 3 - Run NormalyzerDE
*/
process run_normalyzerde {
    label 'process_low'
    label 'process_single_thread'

    publishDir "${params.outdir}/normalyzerde", mode:'copy'

    input:
    file exp_file_parsed from exp_design_parsed
    file exp_file_input from input_exp_design
    file protein_file from protein_groups
    file peptide_file from peptides

    output:
	file "NormalyzerProteins/NormalyzerProteins_stats.tsv" into normalyzer_prot_stats
	file "NormalyzerProteins/${params.normalyzerMethod}-normalized.txt" into normalyzer_prot_normalized
        file "Normalyzer_design.tsv" into normalyzer_design
        file "NormalyzerProteins/*.pdf" into normalyzer_prot_reports
        file "NormalyzerPeptides/NormalyzerPeptides_stats.tsv" into normalyzer_peptide_stats
        file "NormalyzerPeptides/${params.normalyzerMethod}-normalized.txt" into normalyzer_peptide_normalized
        file "NormalyzerPeptides/*.pdf" into normalyzer_peptide_reports 

    when:
    params.run_statistics

    script:
    """
    cp "${exp_file_parsed}" exp_file.tsv
    cp "${exp_file_input}" exp_file2.tsv 
    cp "${protein_file}" protein_file.txt
    cp "${peptide_file}" peptide_file.txt
    Rscript $baseDir/scripts/runNormalyzer.R --comps="${params.comparisons}" --method="${params.normalyzerMethod}"
    """
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the files in the following folder --> $params.outdir\n" : "Oops .. something went wrong" )
}
