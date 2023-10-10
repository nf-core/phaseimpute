/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowPhaseimpute.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { BAM_DOWNSAMPLE                     } from '../subworkflows/local/bam_downsample.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SAMTOOLS_FAIDX              } from '../modules/nf-core/samtools/faidx/main'
include { BAM_REGION                  } from '../subworkflows/local/bam_region'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INITIALIZE PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
//
map      = params.map      ? Channel.fromPath(params.map).collect()     : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT PLUGINS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Add plugins
//
include { fromSamplesheet } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow PHASEIMPUTE {

    // To gather all QC reports for MultiQC
    reports = Channel.empty()

    // To gather used softwares versions for MultiQC
    ch_versions = Channel.empty()

    // Create single fasta channel
    ch_fasta = Channel.of([[genome: params.genome]])
        .combine(Channel.fromPath(params.fasta).collect())

    // Gather regions to use and create the meta map
    if (params.input_region_string == "all") {
        SAMTOOLS_FAIDX(ch_fasta, [[],[]])
        ch_region = SAMTOOLS_FAIDX.output.fai
            .splitCsv(header: ["chr", "size", "offset", "lidebase", "linewidth", "qualoffset"], sep: "\t")
            .map{ meta, row -> [meta + ["chr": row.chr], row.chr + ":0-" + row.size]}
            .map{ metaC, region -> [metaC + ["region": region], region]}
    } else {
        ch_region = Channel.fromSamplesheet("input_region_file")
            .map{ chr, start, end -> [["chr": chr], chr + ":" + start + "-" + end]}
            .map{ metaC, region -> [metaC + ["region": region], region]}
    }

    //
    // Simulate data if asked
    //
    if (params.step == 'simulate') {
        //
        // Read in samplesheet, validate and stage input_simulate files
        //
        ch_input_sim = Channel.fromSamplesheet("input")

        // Split the bam into the region specified
        ch_bam_region = BAM_REGION(ch_input_sim, ch_region, fasta)

        // Initialize channel to impute
        ch_bam_to_impute = Channel.empty()

        if (params.depth) {
            // Create channel from depth parameter
            ch_depth = Channel.fromList(params.depth)

            // Downsample input to desired depth
            BAM_DOWNSAMPLE(ch_bam_region, ch_depth, ch_fasta)
            ch_versions = ch_versions.mix(BAM_DOWNSAMPLE.out.versions.first())

            ch_bam_to_impute = ch_bam_to_impute
                .combine(BAM_DOWNSAMPLE.out.bam_emul)
        }

        if (params.genotype) {
            // Create channel from samplesheet giving the chips snp position
            ch_chip_snp = Channel.fromSamplesheet("input_chip_snp")
            BAM_TO_GENOTYPE(ch_input_sim, ch_chip_snp, ch_fasta)

            ch_bam_to_impute = ch_bam_to_impute
                .combine(BAM_TO_GENOTYPE.out.bam_emul)
        }
    }

    //
    // Prepare panel
    //
    if (params.step == 'panelprep') {
        ch_panel = Channel.fromSamplesheet("input")

        // Remove if necessary "chr"
        if (params.panel_rename = true) {
            ch_panel = VCF_CHR_RENAME(ch_panel, "./assets/chr_rename.txt")
        }

        GET_PANEL(ch_panel)
        ch_versions = ch_versions.mix(GET_PANEL.out.versions.first())

        // Register all panel preparation to csv
    }

    if (params.step.contains("impute")) {
        // Read from panel preparation csv

        if (params.tools.contains("glimpse1")){
            print("Impute with Glimpse1")
            // Glimpse1 subworkflow
        }
        if (params.tools.contains("glimpse2")){
            print("Impute with Glimpse2")
            // Glimpse2 subworkflow
        }
    }

    //
    // MODULE : dump softwares versions
    //
    //CUSTOM_DUMPSOFTWAREVERSIONS (
    //    ch_versions.unique().collectFile(name: 'collated_versions.yml')
    //)

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowPhaseimpute.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowPhaseimpute.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    //ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    //MULTIQC (
    //    ch_multiqc_files.collect(),
    //    ch_multiqc_config.toList(),
    //    ch_multiqc_custom_config.toList(),
    //    ch_multiqc_logo.toList()
    //)
    //multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log) //, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
