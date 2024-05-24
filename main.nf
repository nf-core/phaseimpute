#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/phaseimpute
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/phaseimpute
    Website: https://nf-co.re/phaseimpute
    Slack  : https://nfcore.slack.com/channels/phaseimpute
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PHASEIMPUTE             } from './workflows/phaseimpute'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_phaseimpute_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_phaseimpute_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_phaseimpute_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_PHASEIMPUTE {

    take:
    ch_input       // channel: samplesheet read in from --input
    ch_input_truth // channel: samplesheet read in from --input-truth
    ch_fasta       // channel: reference genome FASTA file with index
    ch_panel       // channel: reference panel variants file
    ch_regions     // channel: regions to use [[chr, region], region]
    ch_depth       // channel: depth of coverage file [[depth], depth]
    ch_map         // channel: map file for imputation
    ch_posfile     // channel: samplesheet read in from --posfile
    ch_chunks      // channel: samplesheet read in from --chunks
    ch_versions    // channel: versions of software used

    main:

    //
    // Initialise input channels
    //

    ch_input_impute         = Channel.empty()
    ch_input_simulate       = Channel.empty()
    ch_input_validate       = Channel.empty()

    if (params.steps.split(',').contains("simulate") || params.steps.split(',').contains("all")) {
        ch_input_simulate = ch_input
    } else if (params.steps.split(',').contains("impute")) {
        ch_input_impute   = ch_input
    } else if (params.steps.split(',').contains("validate")) {
        ch_input_validate = ch_input
    }

    //
    // WORKFLOW: Run pipeline
    //
    PHASEIMPUTE (
        ch_input_impute,
        ch_input_simulate,
        ch_input_validate,
        ch_input_truth,
        ch_fasta,
        ch_panel,
        ch_regions,
        ch_depth,
        ch_map,
        ch_posfile,
        ch_chunks,
        ch_versions
    )


    emit:
    multiqc_report = PHASEIMPUTE.out.multiqc_report // channel: /path/to/multiqc_report.html

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_PHASEIMPUTE (
        PIPELINE_INITIALISATION.out.input,
        PIPELINE_INITIALISATION.out.input_truth,
        PIPELINE_INITIALISATION.out.fasta,
        PIPELINE_INITIALISATION.out.panel,
        PIPELINE_INITIALISATION.out.regions,
        PIPELINE_INITIALISATION.out.depth,
        PIPELINE_INITIALISATION.out.map,
        PIPELINE_INITIALISATION.out.posfile,
        PIPELINE_INITIALISATION.out.chunks,
        PIPELINE_INITIALISATION.out.versions
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_PHASEIMPUTE.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
