/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../../modules/nf-core/multiqc'
include { paramsSummaryMap            } from 'plugin/nf-validation'
include { paramsSummaryMultiqc        } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML      } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../../subworkflows/local/utils_nfcore_phaseimpute_pipeline'
include { getAllFilesExtension        } from '../../subworkflows/local/utils_nfcore_phaseimpute_pipeline'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

// Simulate subworkflows
include { BAM_REGION                                 } from '../../subworkflows/local/bam_region'
include { BAM_DOWNSAMPLE                             } from '../../subworkflows/local/bam_downsample'
include { CHANNEL_SIMULATE_CREATE_CSV                } from '../../subworkflows/local/channel_simulate_create_csv'

// Panelprep subworkflows
include { VCF_CHR_CHECK                              } from '../../subworkflows/local/vcf_chr_check'
include { VCF_NORMALIZE_BCFTOOLS                     } from '../../subworkflows/local/vcf_normalize_bcftools'
include { VCF_SITES_EXTRACT_BCFTOOLS                 } from '../../subworkflows/local/vcf_sites_extract_bcftools'
include { VCF_PHASE_PANEL                            } from '../../subworkflows/local/vcf_phase_panel'
include { CHUNK_PREPARE_CHANNEL                      } from '../../subworkflows/local/chunk_prepare_channel'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_PANEL   } from '../../subworkflows/local/vcf_concatenate_bcftools'
include { CHANNEL_POSFILE_CREATE_CSV                 } from '../../subworkflows/local/channel_posfile_create_csv'
include { CHANNEL_CHUNKS_CREATE_CSV                  } from '../../subworkflows/local/channel_chunks_create_csv'
include { CHANNEL_PANEL_CREATE_CSV                   } from '../../subworkflows/local/channel_panel_create_csv'

// Imputation subworkflows
include { CHANNEL_IMPUTE_CREATE_CSV                   } from '../../subworkflows/local/channel_impute_create_csv'

// GLIMPSE1 subworkflows
include { VCF_IMPUTE_GLIMPSE1                        } from '../../subworkflows/local/vcf_impute_glimpse1'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_GLIMPSE1} from '../../subworkflows/local/vcf_concatenate_bcftools'

// GLIMPSE2 subworkflows
include { VCF_IMPUTE_GLIMPSE2                        } from '../../subworkflows/local/vcf_impute_glimpse2'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_GLIMPSE2} from '../../subworkflows/local/vcf_concatenate_bcftools'

// QUILT subworkflows
include { VCF_CHUNK_GLIMPSE                          } from '../../subworkflows/local/vcf_chunk_glimpse'
include { BAM_IMPUTE_QUILT                           } from '../../subworkflows/local/bam_impute_quilt'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_QUILT   } from '../../subworkflows/local/vcf_concatenate_bcftools'

// STITCH subworkflows
include { PREPARE_INPUT_STITCH                       } from '../../subworkflows/local/prepare_input_stitch'
include { BAM_IMPUTE_STITCH                          } from '../../subworkflows/local/bam_impute_stitch'
include { VCF_SAMPLES_BCFTOOLS                       } from '../../subworkflows/local/vcf_samples_bcftools'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_STITCH  } from '../../subworkflows/local/vcf_concatenate_bcftools'

// Concordance subworkflows
include { BAM_GL_BCFTOOLS as GL_TRUTH                } from '../../subworkflows/local/bam_gl_bcftools'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_TRUTH   } from '../../subworkflows/local/vcf_concatenate_bcftools'
include { VCF_CONCORDANCE_GLIMPSE2                   } from '../../subworkflows/local/vcf_concordance_glimpse2'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PHASEIMPUTE {

    take:
    ch_input_impute         // channel: input file    [ [id], file, index ]
    ch_input_sim            // channel: input file    [ [id], file, index ]
    ch_input_validate       // channel: input file    [ [id], file, index ]
    ch_input_validate_truth // channel: truth file    [ [id], file, index ]
    ch_fasta                // channel: fasta file    [ [genome], fasta, fai ]
    ch_panel                // channel: panel file    [ [id, chr], vcf, index ]
    ch_region               // channel: region to use [ [chr, region], region]
    ch_depth                // channel: depth select  [ [depth], depth ]
    ch_map                  // channel: genetic map   [ [chr], map]
    ch_posfile              // channel: posfile       [ [chr], txt]
    ch_chunks               // channel: chunks        [ [chr], txt]
    ch_versions             // channel: versions of software used

    main:

    ch_multiqc_files = Channel.empty()

    //
    // Simulate data if asked
    //
    if (params.steps.split(',').contains("simulate") || params.steps.split(',').contains("all")) {
        // Output channel of simulate process
        ch_sim_output = Channel.empty()

        // Test if the input are all bam files
        getAllFilesExtension(ch_input_sim)
            .map{ if (it != "bam") {
                error "All input files must be in BAM format to perform simulation"
            } }

        // Split the bam into the region specified
        BAM_REGION(ch_input_sim, ch_region, ch_fasta)
        ch_versions = ch_versions.mix(BAM_REGION.out.versions)

        // Initialize channel to impute
        ch_bam_to_impute = Channel.empty()

        if (params.depth) {
            // Downsample input to desired depth
            BAM_DOWNSAMPLE(
                BAM_REGION.out.bam_region,
                ch_depth,
                ch_fasta
            )
            ch_versions             = ch_versions.mix(BAM_DOWNSAMPLE.out.versions)
            ch_multiqc_files        = ch_multiqc_files.mix(BAM_DOWNSAMPLE.out.coverage.map{ [it[1]] })
            ch_input_impute         = BAM_DOWNSAMPLE.out.bam_emul
            ch_input_validate_truth = ch_input_sim
        }

        if (params.genotype) {
            error "Genotype simulation not yet implemented"
        }

        // Create CSV from simulate step
        CHANNEL_SIMULATE_CREATE_CSV(ch_input_impute, params.outdir)
    }

    //
    // Prepare panel
    //
    if (params.steps.split(',').contains("panelprep") || params.steps.split(',').contains("validate") || params.steps.split(',').contains("all")) {
        // Check chr prefix and remove if necessary
        VCF_CHR_CHECK(ch_panel, ch_fasta)
        ch_versions = ch_versions.mix(VCF_CHR_CHECK.out.versions)

        // Normalize indels in panel
        VCF_NORMALIZE_BCFTOOLS(VCF_CHR_CHECK.out.vcf, ch_fasta)
        ch_versions = ch_versions.mix(VCF_NORMALIZE_BCFTOOLS.out.versions)

        // Extract sites from normalized vcf
        VCF_SITES_EXTRACT_BCFTOOLS(VCF_NORMALIZE_BCFTOOLS.out.vcf_tbi)
        ch_versions = ch_versions.mix(VCF_SITES_EXTRACT_BCFTOOLS.out.versions)

        // If required, phase panel (currently not working, a test should be added)
        // Phase panel with tool of choice (e.g. SHAPEIT5)
        VCF_PHASE_PANEL(VCF_NORMALIZE_BCFTOOLS.out.vcf_tbi)
        ch_versions = ch_versions.mix(VCF_PHASE_PANEL.out.versions)

        // Generate posfile channels
        ch_posfile_glimpse = VCF_SITES_EXTRACT_BCFTOOLS.out.panel_sites
            .join(VCF_SITES_EXTRACT_BCFTOOLS.out.panel_tsv_glimpse)
            .map{ metaPC, sites, s_index, tsv, t_index -> [metaPC, sites, tsv]}

        ch_posfile_stitch = VCF_SITES_EXTRACT_BCFTOOLS.out.panel_tsv_stitch

        CONCAT_PANEL(VCF_SITES_EXTRACT_BCFTOOLS.out.panel_sites)
        ch_versions    = ch_versions.mix(CONCAT_PANEL.out.versions)

        ch_panel_sites = CONCAT_PANEL.out.vcf_tbi_join
        ch_panel_phased = VCF_PHASE_PANEL.out.vcf_tbi

        // Create chunks from reference VCF
        VCF_CHUNK_GLIMPSE(ch_panel_phased, ch_map)
        ch_versions = ch_versions.mix(VCF_CHUNK_GLIMPSE.out.versions)

        // Create CSVs from panelprep step
        CHANNEL_POSFILE_CREATE_CSV(VCF_SITES_EXTRACT_BCFTOOLS.out.panel_tsv_stitch, params.outdir)
        CHANNEL_CHUNKS_CREATE_CSV(VCF_CHUNK_GLIMPSE.out.chunks, params.outdir)
        CHANNEL_PANEL_CREATE_CSV(ch_panel_phased,
                VCF_NORMALIZE_BCFTOOLS.out.hap_legend,
                params.outdir)

    }

    if (params.steps.split(',').contains("impute") || params.steps.split(',').contains("all")) {
            if (params.tools.split(',').contains("glimpse1")) {
                print("Impute with GLIMPSE1")

                if (params.chunks) {
                    ch_chunks_glimpse1 = CHUNK_PREPARE_CHANNEL(ch_chunks, "glimpse").out.chunks
                } else if (params.panel && params.steps.split(',').find { it in ["all", "panelprep"] } && !params.chunks) {
                    ch_chunks_glimpse1 = VCF_CHUNK_GLIMPSE.out.chunks_glimpse1
                }

                // if (params.posfile) {
                //     ch_posfile_glimpse = // User supplied posfile
                // } else if (params.panel && params.steps.split(',').find { it in ["all", "panelprep"] } && !params.posfile) {
                //     ch_posfile_glimpse = POSFILE_PREPARE_CHANNEL(ch_posfile, "glimpse").out.posfile
                // }

                // if (params.panel && !params.steps.split(',').find { it in ["all", "panelprep"] }) {
                //     ch_panel_phased = // User supplied phased panel
                // }

                VCF_IMPUTE_GLIMPSE1(
                    ch_input_impute,
                    ch_posfile_glimpse,
                    ch_panel_phased,
                    ch_chunks_glimpse1,
                    ch_fasta
                )
                ch_versions = ch_versions.mix(VCF_IMPUTE_GLIMPSE1.out.versions)
                ch_multiqc_files = ch_multiqc_files.mix(VCF_IMPUTE_GLIMPSE1.out.multiqc_files)

                // Concatenate by chromosomes
                CONCAT_GLIMPSE1(VCF_IMPUTE_GLIMPSE1.out.vcf_tbi)
                ch_versions = ch_versions.mix(CONCAT_GLIMPSE1.out.versions)

                // Add results to input validate
                ch_input_validate = ch_input_validate.mix(CONCAT_GLIMPSE1.out.vcf_tbi_join)

            }
            if (params.tools.split(',').contains("glimpse2")) {
                print("Impute with GLIMPSE2")

                // Use previous chunks if --steps panelprep
                if (params.panel && params.steps.split(',').find { it in ["all", "panelprep"] } && !params.chunks) {
                    ch_chunks = VCF_CHUNK_GLIMPSE.out.chunks_glimpse1 // Chunks from glimpse2 are wrong
                } else if (params.chunks) {
                    ch_chunks = CHUNK_PREPARE_CHANNEL(ch_chunks, "glimpse").out.chunks
                }

                // if (params.posfile) {
                //     ch_posfile_glimpse = // User supplied posfile
                // } else if (params.panel && params.steps.split(',').find { it in ["all", "panelprep"] } && !params.posfile) {
                //     ch_posfile_glimpse = POSFILE_PREPARE_CHANNEL(ch_posfile, "glimpse").out.posfile
                // }

                // if (params.panel && !params.steps.split(',').find { it in ["all", "panelprep"] }) {
                //     ch_panel_phased = // User supplied phased panel
                // }

                // Run imputation
                VCF_IMPUTE_GLIMPSE2(
                    ch_input_impute,
                    ch_panel_phased,
                    ch_chunks,
                    ch_fasta
                )
                ch_versions = ch_versions.mix(VCF_IMPUTE_GLIMPSE2.out.versions)
                // Concatenate by chromosomes
                CONCAT_GLIMPSE2(VCF_IMPUTE_GLIMPSE2.out.vcf_tbi)
                ch_versions = ch_versions.mix(CONCAT_GLIMPSE2.out.versions)

                // Add results to input validate
                ch_input_validate = ch_input_validate.mix(CONCAT_GLIMPSE2.out.vcf_tbi_join)
            }
            if (params.tools.split(',').contains("stitch")) {
                print("Impute with STITCH")

                if (params.posfile) {
                    ch_posfile_stitch = ch_posfile
                }

                // Prepare inputs
                PREPARE_INPUT_STITCH(ch_input_impute, ch_posfile_stitch, ch_region)
                ch_versions = ch_versions.mix(PREPARE_INPUT_STITCH.out.versions)

                // Impute with STITCH
                BAM_IMPUTE_STITCH (
                    PREPARE_INPUT_STITCH.out.stitch_parameters,
                    PREPARE_INPUT_STITCH.out.stitch_samples,
                    ch_fasta
                )
                ch_versions = ch_versions.mix(BAM_IMPUTE_STITCH.out.versions)

                // Concatenate by chromosomes
                CONCAT_STITCH(BAM_IMPUTE_STITCH.out.vcf_tbi)
                ch_versions = ch_versions.mix(CONCAT_STITCH.out.versions)

                // Separate by samples
                VCF_SAMPLES_BCFTOOLS(CONCAT_STITCH.out.vcf_tbi_join)
                ch_versions = ch_versions.mix(VCF_SAMPLES_BCFTOOLS.out.versions)

                // Add results to input validate
                ch_input_validate = ch_input_validate.mix(VCF_SAMPLES_BCFTOOLS.out.vcf_tbi_join)

            }
            if (params.tools.split(',').contains("quilt")) {
                print("Impute with QUILT")

                // Use previous chunks if --steps panelprep
                if (params.panel && params.steps.split(',').find { it in ["all", "panelprep"] } && !params.chunks) {
                    ch_chunks_quilt = VCF_CHUNK_GLIMPSE.out.chunks_quilt
                // Use provided chunks if --chunks
                } else if (params.chunks) {
                    ch_chunks_quilt = CHUNK_PREPARE_CHANNEL(ch_chunks, "quilt").out.chunks
                }

                // Impute BAMs with QUILT
                BAM_IMPUTE_QUILT(
                    ch_input_impute,
                    VCF_NORMALIZE_BCFTOOLS.out.hap_legend,
                    VCF_CHUNK_GLIMPSE.out.chunks_quilt
                )
                ch_versions = ch_versions.mix(BAM_IMPUTE_QUILT.out.versions)

                // Concatenate by chromosomes
                CONCAT_QUILT(BAM_IMPUTE_QUILT.out.vcf_tbi)
                ch_versions = ch_versions.mix(CONCAT_QUILT.out.versions)

                // Add results to input validate
                ch_input_validate = ch_input_validate.mix(CONCAT_QUILT.out.vcf_tbi_join)
            }
            // Create CSV from imputation step
            CHANNEL_IMPUTE_CREATE_CSV(ch_input_validate, params.outdir)
        }

    if (params.steps.split(',').contains("validate") || params.steps.split(',').contains("all")) {
        ch_truth_vcf = Channel.empty()
        // Get extension of input files
        truth_ext = getAllFilesExtension(ch_input_validate_truth)

        // Channels for branching
        ch_truth = ch_input_validate_truth
            .combine(truth_ext)
            .branch {
                bam: it[3] == 'bam'
                vcf: it[3] =~ 'vcf|bcf'
            }

        GL_TRUTH(
            ch_truth.bam.map { [it[0], it[1], it[2]] },
            ch_posfile_glimpse,
            ch_fasta
        )
        ch_multiqc_files = ch_multiqc_files.mix(GL_TRUTH.out.multiqc_files)
        ch_versions      = ch_versions.mix(GL_TRUTH.out.versions)

        // Mix the original vcf and the computed vcf
        ch_truth_vcf = ch_truth.vcf
            .map { [it[0], it[1], it[2]] }
            .mix(GL_TRUTH.out.vcf)

        // Concatenate truth vcf by chromosomes
        CONCAT_TRUTH(ch_truth_vcf)
        ch_versions = ch_versions.mix(CONCAT_TRUTH.out.versions)

        // Compute concordance analysis
        VCF_CONCORDANCE_GLIMPSE2(
            ch_input_validate,
            CONCAT_TRUTH.out.vcf_tbi_join,
            ch_panel_sites,
            ch_region
        )
        ch_multiqc_files = ch_multiqc_files.mix(VCF_CONCORDANCE_GLIMPSE2.out.multiqc_files)
        ch_versions      = ch_versions.mix(VCF_CONCORDANCE_GLIMPSE2.out.versions)
    }

    if (params.steps.split(',').contains("refine")) {
        error "refine steps not yet implemented"
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
