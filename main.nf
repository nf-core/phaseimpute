#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/phaseimpute
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/phaseimpute
    Website: https://nf-co.re/phaseimpute
    Slack  : https://nfcore.slack.com/channels/phaseimpute
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PHASEIMPUTE                } from './workflows/phaseimpute'
include { CHRCHECK as CHRCHECK_INPUT } from './workflows/chrcheck'
include { CHRCHECK as CHRCHECK_TRUTH } from './workflows/chrcheck'
include { CHRCHECK as CHRCHECK_PANEL } from './workflows/chrcheck'
include { PIPELINE_INITIALISATION    } from './subworkflows/local/utils_nfcore_phaseimpute_pipeline'
include { PIPELINE_COMPLETION        } from './subworkflows/local/utils_nfcore_phaseimpute_pipeline'
include { PREPARE_GENOME             } from './subworkflows/local/prepare_genome'
include { getGenomeAttribute         } from 'plugin/nf-core-utils'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta             = getGenomeAttribute('fasta')
params.fasta_fai         = getGenomeAttribute('fasta_fai')
params.fasta_gzi         = getGenomeAttribute('fasta_gzi')

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
    steps            // array: List of steps to perform
    tools            // array: List of tools to use
    ch_input         // channel: samplesheet read in from --input
    ch_input_truth   // channel: samplesheet read in from --input-truth
    ch_fasta         // channel: reference genome FASTA file with index
    ch_panel         // channel: reference panel variants file
    ch_regions       // channel: regions to use [[chr, region], region]
    ch_depth         // channel: depth of coverage file [[depth], depth]
    ch_map           // channel: map file for imputation
    ch_posfile       // channel: samplesheet read in from --posfile
    ch_chunks        // channel: samplesheet read in from --chunks
    sheets_given     // map: input sheets given
    rename_chr       // parameter: rename chromosome prefix
    max_chr_names    // parameter: max number of chr to show in message
    params_simulate  // map: parameters use for simulation step [depth: float, genotype: path]
    params_panelprep // map: parameters use for panelprep step  [normalize: boolean, remove_samples: string, compute_freq: boolean, phase: boolean, chunk_model: string ]
    params_impute    // map: parameters use for imputation step [batch_size: integer, k_val: integer, n_gen: integer, buffer: integer]
    params_validate  // map: parameters use for validation step [bins: string, min_val_gl: float, min_val_dp: integer]
    params_multiqc   // map: parameters use for multiqc report  [config: path, logo: path, methods_description: string]
    seed             // integer
    outdir           // path

    main:

    //
    // Initialise input channels
    //

    ch_input_impute         = channel.empty()
    ch_input_simulate       = channel.empty()
    ch_input_validate       = channel.empty()

    //  Check input files for contigs names consistency
    lst_chr = ch_regions.map {meta, _region -> meta.chr }
        .unique()
        .collect()
        .toList()

    CHRCHECK_INPUT(ch_input.combine(lst_chr), rename_chr, max_chr_names)
    ch_input = CHRCHECK_INPUT.out.output

    CHRCHECK_TRUTH(ch_input_truth.combine(lst_chr), rename_chr, max_chr_names)
    ch_input_truth = CHRCHECK_TRUTH.out.output

    CHRCHECK_PANEL(
        ch_panel.map{ meta, file, index -> [meta, file, index, [meta.chr]]},
        rename_chr, max_chr_names
    )
    ch_panel = CHRCHECK_PANEL.out.output

    if (steps.contains("simulate") || steps.contains("all")) {
        ch_input_simulate = ch_input
    } else if (steps.contains("impute")) {
        ch_input_impute   = ch_input
    } else if (steps.contains("validate")) {
        ch_input_validate = ch_input
    }

    //
    // WORKFLOW: Run pipeline
    //
    PHASEIMPUTE (
        steps,
        tools,
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
        sheets_given,
        params_simulate,
        params_panelprep,
        params_impute,
        params_validate,
        params_multiqc,
        seed,
        outdir
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
    def steps = params.steps.split(',') as List
    def tools = params.tools ? params.tools.split(',') as List : []

    PREPARE_GENOME(
        params.genome ?: file(params.fasta, checkIfExists:true).getBaseName(),
        params.fasta,
        params.fasta_fai,
        params.fasta_gzi
    )

    def sheets_given = [
        input_target : params.input,
        input_truth  : params.input_truth,
        input_region : params.input_region,
        input_panel  : params.panel,
        input_posfile: params.posfile,
        input_chunks : params.chunks,
        input_map    : params.map,
    ]

    def params_simulate = [
        depth   : params.depth,
        genotype: params.genotype
    ]

    def params_panelprep = [
        normalize     : params.normalize,
        remove_samples: params.remove_samples,
        compute_freq  : params.compute_freq,
        phase         : params.phase,
        chunk_model   : params.chunk_model
    ]

    def params_impute = [
        batch_size: params.batch_size,
        k_val     : params.k_val,
        n_gen     :params.n_gen,
        buffer    :params.buffer,
    ]

    def params_validate = [
        bins      : params.bins,
        min_val_gl: params.min_val_gl,
        min_val_dp: params.min_val_gl
    ]

    def params_multiqc = [
        config             : params.multiqc_config,
        logo               : params.multiqc_logo,
        methods_description: params.multiqc_methods_description
    ]

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.help,
        params.help_full,
        params.show_hidden,
        PREPARE_GENOME.out.ch_fasta_fai_gzi,
        sheets_given,
        steps,
        tools,
        params.max_chr_names,
        params_simulate,
        params_panelprep,
        params_impute
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_PHASEIMPUTE (
        steps,
        tools,
        PIPELINE_INITIALISATION.out.ch_input_target,
        PIPELINE_INITIALISATION.out.ch_input_truth,
        PIPELINE_INITIALISATION.out.ch_fasta_index,
        PIPELINE_INITIALISATION.out.ch_panel,
        PIPELINE_INITIALISATION.out.ch_regions,
        PIPELINE_INITIALISATION.out.ch_depth,
        PIPELINE_INITIALISATION.out.ch_map,
        PIPELINE_INITIALISATION.out.ch_posfile,
        PIPELINE_INITIALISATION.out.ch_chunks,
        sheets_given,
        params.rename_chr,
        params.max_chr_names,
        params_simulate,
        params_panelprep,
        params_impute,
        params_validate,
        params_multiqc,
        params.seed,
        params.outdir,
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
