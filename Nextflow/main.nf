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
params.raws = params.raws ?: { log.error "No read data provided. Make sure you have used the '--raws' option."; exit 1 }()
params.sdrf = params.sdrf ?: {log.error "No read data provided. Make sure you have used the '--sdrf' option."; exit 1}()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()
params.fasta = params.fasta ?: { log.error "No fasta file provided. Make sure you have used the '--fasta' option."; exit 1}()
 
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
Channel 
    .fromPath (params.raws).into {input_raw; input_raw2}

/*
* Generate the channels for the sdrf files
*/ 
Channel 
    .fromPath (params.sdrf).into {input_sdrf; input_sdrf2}

/*
* Generate the channels for the fasta file
*/ 
Channel
    .fromPath (params.fasta).into {input_fasta; input_fasta2}

/*
* Generate channel for experimental design file
*/
Channel
     .fromPath(params.experiment_design).into {input_exp_design}

/* 
* STEP 1 - Generate the parameter and experimental design for maxquant through sdrf
*/
process run_sdrf {
    publishDir "${params.outdir}"
    input: 
        path sdrf_file from input_sdrf
        path fasta_file from input_fasta
    
    output:
        file "mqpar.xml" into outParameters
        file "exp_design.tsv" into outExpDesign

    script: 
    """
    parse_sdrf convert-maxquant -s "${sdrf_file}" -f "PLACEHOLDER${fasta_file}" -m ${params.match} -r PLACEHOLDER -pef ${params.peptidefdr} -prf ${params.proteinfdr} -t PLACEHOLDERtemp -o2 exp_design.tsv -n ${task.cpus} 
    """
}

/*
* STEP 2 - Run maxQuant with parameters from mqpar
*/
process run_maxquant {
   
    publishDir "${params.outdir}"

    input:
        file rawfile from input_raw2.collect()
	file fastafile from input_fasta2
        path mqparameters from outParameters

    output:
        file "combined/txt/proteinGroups.txt" into input_proteinGroups	


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

    publishDir "${params.outdir}"

    input:
        path sdrf_file from input_sdrf2
        path exp_file from outExpDesign
        path protein_file from input_proteinGroups
        file exp_file2 from input_exp_design

    output:
	file "Normalyzer/Normalyzer_stats.tsv" into normalyzer_out
	file "Normalyzer/${params.normalyzerMethod}-normalized.txt" into normalyzer_out2

    when:
      params.run_statistics

    script:
    """
     cp "${exp_file}" exp_file.tsv
     cp "${exp_file2}" exp_file2.tsv 
     cp "${protein_file}" protein_file.txt
     Rscript $baseDir/runNormalyzer.R --comps="${params.comparisons}" --method="${params.normalyzerMethod}"
     cp -R Normalyzer "${params.outdir}"
    """   
}


workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the files in the following folder --> $params.outdir\n" : "Oops .. something went wrong" )
}
