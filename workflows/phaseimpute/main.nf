/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../../modules/nf-core/multiqc/main'
include { paramsSummaryMap            } from 'plugin/nf-validation'
include { paramsSummaryMultiqc        } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML      } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../../subworkflows/local/utils_nfcore_phaseimpute_pipeline'

include { BAM_REGION                  } from '../../subworkflows/local/bam_region'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { BAM_DOWNSAMPLE              } from '../../subworkflows/local/bam_downsample'
include { COMPUTE_GL as GL_TRUTH      } from '../../subworkflows/local/compute_gl'
include { COMPUTE_GL as GL_INPUT      } from '../../subworkflows/local/compute_gl'
include { VCF_CHR_RENAME              } from '../../subworkflows/local/vcf_chr_rename'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TEST PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Assert that the different parameters are set correctly for each step
//

if (params.step.contains("impute")) {
    assert params.tools, "You must specify at least one imputation tool to use"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PHASEIMPUTE {

    take:
    ch_input       // channel: samplesheet read in from --input
    ch_fasta       // channel: fasta file
    ch_panel       // channel: panel file
    ch_region      // channel: region to use [meta, region]
    ch_map         // channel: genetic map

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    /*
    //
    // Simulate data if asked
    //
    if (params.step == 'simulate') {
        //
        // Read in samplesheet, validate and stage input_simulate files
        //
        ch_sim_input = Channel.fromSamplesheet("input")

        // Output channel of simulate process
        ch_sim_output = Channel.empty()

        // Split the bam into the region specified
        ch_bam_region = BAM_REGION(ch_input_sim, ch_region, fasta)

        // Initialize channel to impute
        ch_bam_to_impute = Channel.empty()

        if (params.depth) {
            // Create channel from depth parameter
            ch_depth = Channel.fromList(params.depth)

            // Downsample input to desired depth
            BAM_DOWNSAMPLE(ch_sim_input, ch_region, ch_depth, ch_fasta)
            ch_versions = ch_versions.mix(BAM_DOWNSAMPLE.out.versions.first())

            ch_sim_output = ch_sim_output.mix(BAM_DOWNSAMPLE.out.bam_emul)
        }

        if (params.genotype) {
            // Create channel from samplesheet giving the chips snp position
            ch_chip_snp = Channel.fromSamplesheet("input_chip_snp")
            BAM_TO_GENOTYPE(ch_sim_input, ch_region, ch_chip_snp, ch_fasta)
            ch_sim_output = ch_sim_output.mix(BAM_TO_GENOTYPE.out.bam_emul)
        }
    }
    */
    //
    // Prepare panel
    //
    if (params.step == 'impute') {
        // Remove if necessary "chr"
        if (params.panel_rename = true) {
            ch_panel = VCF_CHR_RENAME(ch_panel, "./assets/chr_rename.txt")
        }

        GET_PANEL(ch_panel)
        ch_versions = ch_versions.mix(GET_PANEL.out.versions.first())

        // Register all panel preparation to csv
        ch_panel_sites  = GET_PANEL.out.panel_sites
        ch_panel_tsv    = GET_PANEL.out.panel_tsv
        ch_panel_phased = GET_PANEL.out.panel_phased

        // Output channel of input process
        ch_impute_output = Channel.empty()

        if (params.tools.contains("glimpse1")) {
            print("Impute with Glimpse1")
            // Glimpse1 subworkflow
            GL_INPUT(
                ch_input,
                ch_region,
                ch_panel_sites,
                ch_panel_tsv
            )
            
            impute_input = GL_INPUT.out.vcf
                | combine(Channel.of([[]]))
                | map{meta, vcf, index, sample -> [meta, vcf, index, sample, meta.region]}

            VCF_IMPUTE_GLIMPSE(impute_input,
                ch_panel.phased,
                ch_map)
            
            ch_impute_output = ch_impute_output.mix(VCF_IMPUTE_GLIMPSE.out.merged_variants)
        }
        if (params.tools.contains("glimpse2")) {
            print("Impute with Glimpse2")
            error "Glimpse2 not yet implemented"
            // Glimpse2 subworkflow
        }
        if (params.tools.contains("quilt")) {
            print("Impute with quilt")
            error "Quilt not yet implemented"
            // Quilt subworkflow
        }
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
