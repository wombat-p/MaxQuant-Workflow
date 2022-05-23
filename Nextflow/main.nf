#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/wombatp
========================================================================================
 nf-core/wombatp Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/wombatp
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2


def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info nfcoreHeader()
    log.info"""

    Usage:

    There are two ways to run this workflow:
    a) Using an sdrf file that contains paths to the raw spectrum files
    The typical command for running the pipeline is as follows:
    nextflow run nf-core/wombatp --sdrf "path/to/sdrffile" --fasta "path/to/sdrfile" --parameters "path/to/param_yml_file" -profile docker

    b) Using a file with the experimental design and corresponding raw files
    The typical command for running the pipeline is as follows:
    nextflow run nf-core/wombatp --raws "path/to/rawfiles" --exp_design "path/to/exp_design_file" --fasta "path/to/sdrfile" --parameters "path/to/param_yml_file" -profile docker

    


    Mandatory arguments (either "--raws" or "--sdrf" containing valid file paths:
      --raws [file]                   Path to folder with raw data files or wildcard (must be surrounded with single quotes)
      --sdrf			      Path to sdrf file for experimental design and some of the parameters. Requires valid (online) paths to raw files.
      --fasta			      Path to fasta database file

    Optional argument (recommended):
      --parameters		      Path to yaml file containing the workflow parameters  (see example in data folder)
      --exp_design		      Path to tsv file with experimental design (see data folder for an example). This parameter is only needed when not specifying an sdrf file

    For defining environment (mandatory):  
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated)
                                      Available: conda, docker, singularity, test, awsbatch, <institute> and more

    Other options:
      --outdir [file]                 The output directory where the results will be saved
      --publish_dir_mode [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Add any reference files that are needed
// 

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}



// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
//ch_sdrfmapping = file("$projectDir/assets/param2sdrf.yml", checkIfExists: true)
ch_sdrfmapping = file("https://raw.githubusercontent.com/bigbio/proteomics-metadata-standard/master/sdrf-proteomics/assets/param2sdrf.yml", checkIfExists: true)
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)




/*
================================================================================
                                 Create channels
================================================================================
*/

/*
 * Create a channel for input read files
 */


params.parameters = "no_params"
params.exp_design = "no_exp_design"
params.sdrf = "no_sdrf"
params.raws = "no_raws"
params.fasta = "no_fasta"

input_raws = Channel
     .fromPath(params.raws)

if (params.raws == "no_raws" & !file(params.sdrf).exists()) {
	exit 1, "ERROR: neither raw files found nor sdrf file given. This could also be due to an empty folder - no input raw files supplied" } 


if(!file(params.fasta).exists()) {
	exit 1, "ERROR: no fasta database specified or file not existing "}
input_fasta = Channel.from(file(params.fasta, checkIfExists: true))

input_sdrf = Channel
     .fromPath(params.sdrf)

input_parameters = Channel.fromPath(params.parameters)
input_exp_design = Channel.fromPath(params.exp_design)

// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Raw files']            = params.raws
summary['Fasta database']       = params.fasta
summary['SDRF file']		= params.sdrf
summary['Parameter file']	= params.parameters
summary['Experimental design file']	= params.exp_design
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-wombatp-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/wombatp Workflow Summary'
    section_href: 'https://github.com/nf-core/wombatp'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }


/*
================================================================================
                                 Modules and workflow
================================================================================
*/


params.collected_options = [:]
//include { GET_SOFTWARE_VERSIONS }         from './modules/local/getsoftwareversions/main'   addParams(options: [publish_files : ['csv':'']])
include { PREPARE_FILES }                 from './modules/local/prepare_files/main'   addParams(options: params.collected_options)
include { CONVERT_MAXQUANT }                 from './modules/local/sdrfpipelines/convert_maxquant/main'   addParams(options: params.collected_options)
include { SDRFMERGE }                 from './modules/local/sdrfpipelines/sdrfmerge/main'         addParams(options: params.collected_options)
include { MAXQUANT_LFQ }                      from './modules/nf-core/modules/maxquant/lfq/main'              addParams(options: [:] )
include { NORMALYZERDE }                  from './modules/local/normalyzerde/main'          addParams(options: params.collected_options)


workflow maxquantpipeline {
    ch_software_versions = Channel.empty()

    PREPARE_FILES (input_sdrf, input_parameters, input_exp_design, input_raws.collect(), ch_sdrfmapping)
    SDRFMERGE (PREPARE_FILES.out.sdrf_local, PREPARE_FILES.out.params, ch_sdrfmapping)
    CONVERT_MAXQUANT (SDRFMERGE.out.sdrf_local, input_fasta)
 //   input = [ [ id:'run' ], // meta map
 //              input_fasta, CONVERT_MAXQUANT.out.maxquantpar ]
    MAXQUANT_LFQ ( input_fasta , CONVERT_MAXQUANT.out.maxquantpar, PREPARE_FILES.out.raws.collect() )
//    MAXQUANT_LFQ (input, input_raws.collect())
    NORMALYZERDE (MAXQUANT_LFQ.out.maxquant_txt, PREPARE_FILES.out.exp_design, CONVERT_MAXQUANT.out.exp_design)
 //   ch_software_versions = ch_software_versions.mix(NORMALIZERDE.out.version.first().ifEmpty(null))
//    GET_SOFTWARE_VERSIONS (
//        ch_software_versions.map { it }.collect()
//    )
}


////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {
maxquantpipeline()
}



/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs
    file images from ch_output_docs_images

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/wombatp] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/wombatp] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.max_multiqc_email_size)
    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/wombatp] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/wombatp] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/wombatp] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/wombatp] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/wombatp]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/wombatp]${c_red} Pipeline completed with errors${c_reset}-"
    }

}


def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/wombatp v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}

