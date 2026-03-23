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
    steps
    tools
    ch_input       // channel: samplesheet read in from --input
    ch_input_truth // channel: samplesheet read in from --input-truth
    ch_fasta       // channel: reference genome FASTA file with index
    ch_panel       // channel: reference panel variants file
    ch_regions     // channel: regions to use [[chr, region], region]
    ch_depth       // channel: depth of coverage file [[depth], depth]
    ch_map         // channel: map file for imputation
    ch_posfile     // channel: samplesheet read in from --posfile
    ch_chunks      // channel: samplesheet read in from --chunks
    chunk_model    // parameter: chunk model
    rename_chr     // parameter: rename chromosome prefix
    max_chr_names  // parameter: max number of chr to show in message
    input_region
    input_truth
    posfile
    chunks
    panel
    depth
    normalize
    remove_samples
    compute_freq
    phase
    batch_size
    k_val
    ngen
    buffer
    bins
    min_val_gl
    min_val_dp
    seed
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

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
        chunk_model,
        input_region,
        input_truth,
        posfile,
        chunks,
        panel,
        depth,
        normalize,
        remove_samples,
        compute_freq,
        phase,
        batch_size,
        k_val,
        ngen,
        buffer,
        bins,
        min_val_gl,
        min_val_dp,
        seed,
        multiqc_config,
        multiqc_logo,
        multiqc_methods_description,
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

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
        PREPARE_GENOME.out.ch_fasta_fai_gzi,
        params.input,
        params.input_truth,
        params.input_region,
        params.panel,
        params.posfile,
        params.chunks,
        params.map,
        steps,
        tools,
        params.batch_size,
        params.max_chr_names,
        params.depth,
        params.genotype,
        params.remove_samples,
        params.normalize,
        params.chunk_model
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_PHASEIMPUTE (
        steps,
        tools,
        PIPELINE_INITIALISATION.out.input,
        PIPELINE_INITIALISATION.out.input_truth,
        PIPELINE_INITIALISATION.out.fasta,
        PIPELINE_INITIALISATION.out.panel,
        PIPELINE_INITIALISATION.out.regions,
        PIPELINE_INITIALISATION.out.depth,
        PIPELINE_INITIALISATION.out.gmap,
        PIPELINE_INITIALISATION.out.posfile,
        PIPELINE_INITIALISATION.out.chunks,
        params.chunk_model,
        params.rename_chr,
        params.max_chr_names,
        params.input_region,
        params.input_truth,
        params.posfile,
        params.chunks,
        params.panel,
        params.depth,
        params.normalize,
        params.remove_samples,
        params.compute_freq,
        params.phase,
        params.batch_size,
        params.k_val,
        params.ngen,
        params.buffer,
        params.bins,
        params.min_val_gl,
        params.min_val_dp,
        params.seed,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir
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
