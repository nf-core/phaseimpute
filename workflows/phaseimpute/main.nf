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
include { getAllFilesExtension        } from '../../subworkflows/local/utils_nfcore_phaseimpute_pipeline'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

// Simulate subworkflows
include { BAM_REGION                                 } from '../../subworkflows/local/bam_region'
include { BAM_DOWNSAMPLE                             } from '../../subworkflows/local/bam_downsample'

// Panelprep subworkflows
include { VCF_CHR_CHECK                              } from '../../subworkflows/local/vcf_chr_check'
include { VCF_NORMALIZE_BCFTOOLS                     } from '../../subworkflows/local/vcf_normalize_bcftools/vcf_normalize_bcftools'
include { PANEL_PREPARE_CHANNELS                     } from '../../subworkflows/local/panel_prepare_channels'
include { VCF_SITES_EXTRACT_BCFTOOLS                 } from '../../subworkflows/local/vcf_sites_extract_bcftools'
include { VCF_PHASE_PANEL                            } from '../../subworkflows/local/vcf_phase_panel'

// GLIMPSE subworkflows
include { VCF_IMPUTE_GLIMPSE as VCF_IMPUTE_GLIMPSE1  } from '../../subworkflows/nf-core/vcf_impute_glimpse'
include { COMPUTE_GL as GL_TRUTH                     } from '../../subworkflows/local/compute_gl'
include { COMPUTE_GL as GL_INPUT                     } from '../../subworkflows/local/compute_gl'

// QUILT subworkflows
include { MAKE_CHUNKS                                } from '../../subworkflows/local/make_chunks/make_chunks'
include { IMPUTE_QUILT                               } from '../../subworkflows/local/impute_quilt/impute_quilt'

// STITCH subworkflows
include { PREPARE_INPUT_STITCH                       } from '../../subworkflows/local/prepare_input_stitch/prepare_input_stitch'
include { BAM_IMPUTE_STITCH                          } from '../../subworkflows/local/bam_impute_stitch/bam_impute_stitch'

// CONCAT subworkflows
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_IMPUT   } from '../../subworkflows/local/vcf_concatenate_bcftools'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_TRUTH   } from '../../subworkflows/local/vcf_concatenate_bcftools'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_PANEL   } from '../../subworkflows/local/vcf_concatenate_bcftools'

// Concordance subworkflows
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
    ch_panel                // channel: panel file    [ [id, chr], chr, vcf, index ]
    ch_region               // channel: region to use [ [chr, region], region]
    ch_depth                // channel: depth select  [ [depth], depth ]
    ch_map                  // channel: genetic map   [ [chr], map]
    ch_versions             // channel: versions of software used

    main:

    ch_multiqc_files = Channel.empty()

    //
    // Simulate data if asked
    //
    if (params.step == 'simulate' || params.step == 'all') {
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
            ch_input_validate_truth = BAM_REGION.out.bam_region
        }

        if (params.genotype) {
            error "Genotype simulation not yet implemented"
        }
    }

    //
    // Prepare panel
    //
    if (params.step == 'impute' || params.step == 'panel_prep' || params.step == 'validate' || params.step == 'all') {
        // Check chr prefix and remove if necessary
        VCF_CHR_CHECK(ch_panel, ch_fasta)
        ch_versions = ch_versions.mix(VCF_CHR_CHECK.out.versions)

        // Normalize indels in panel
        VCF_NORMALIZE_BCFTOOLS(VCF_CHR_CHECK.out.vcf, ch_fasta)
        ch_versions = ch_versions.mix(VCF_NORMALIZE_BCFTOOLS.out.versions)

        // Extract sites from normalized vcf
        VCF_SITES_EXTRACT_BCFTOOLS(VCF_NORMALIZE_BCFTOOLS.out.vcf_tbi)
        ch_versions = ch_versions.mix(VCF_SITES_EXTRACT_BCFTOOLS.out.versions)

        // Phase panel
        VCF_PHASE_PANEL(VCF_SITES_EXTRACT_BCFTOOLS.out.vcf_tbi,
                        VCF_SITES_EXTRACT_BCFTOOLS.out.vcf_tbi,
                        VCF_SITES_EXTRACT_BCFTOOLS.out.panel_sites,
                        VCF_SITES_EXTRACT_BCFTOOLS.out.panel_tsv)
        ch_versions = ch_versions.mix(VCF_PHASE_PANEL.out.versions)

        // Generate channels (to be simplified)
        ch_panel_sites_tsv = VCF_PHASE_PANEL.out.panel
                        .map{ metaPC, norm, n_index, sites, s_index, tsv, t_index, phased, p_index
                        -> [metaPC, sites, tsv]
                        }
        CONCAT_PANEL(VCF_PHASE_PANEL.out.panel
                        .map{ metaPC, norm, n_index, sites, s_index, tsv, t_index, phased, p_index
                        -> [[id:metaPC.panel], sites, s_index]
                        }
        )
        ch_panel_sites = CONCAT_PANEL.out.vcf_tbi_join
        ch_versions    = ch_versions.mix(CONCAT_PANEL.out.versions)

        ch_panel_phased = VCF_PHASE_PANEL.out.panel
            .map{ metaPC, norm, n_index, sites, s_index, tsv, t_index, phased, p_index
                -> [metaPC, phased, p_index]
            }

        if (params.step == 'impute' || params.step == 'all') {
            // Output channel of input process
            ch_impute_output = Channel.empty()
            if (params.tools.split(',').contains("glimpse1")) {
                println "Impute with Glimpse1"
                // Glimpse1 subworkflow
                GL_INPUT( // Compute GL for input data once per panel
                    ch_input_impute,
                    ch_panel_sites_tsv,
                    ch_fasta
                )
                ch_multiqc_files = ch_multiqc_files.mix(GL_INPUT.out.multiqc_files)
                ch_versions = ch_versions.mix(GL_INPUT.out.versions)

                impute_input = GL_INPUT.out.vcf // [metaIPC, vcf, index]
                    .map {metaIPC, vcf, index -> [metaIPC.subMap("panel", "chr"), metaIPC, vcf, index] }
                    .combine(ch_panel_phased, by: 0)
                    .combine(Channel.of([[]]))
                    .map { metaPC, metaIPC, vcf, index, panel, p_index, sample ->
                        [metaPC.subMap("chr"), metaIPC, vcf, index, panel, p_index, sample]}
                    .combine(ch_region
                        .map {metaCR, region -> [metaCR.subMap("chr"), metaCR, region]},
                        by: 0)
                    .combine(ch_map, by: 0)
                    .map{
                        metaC, metaIPC, vcf, index, panel, p_index, sample, metaCR, region, map
                        -> [metaIPC+metaCR.subMap("Region"), vcf, index, sample, region, panel, p_index, map]
                    } //[ metaIPCR, vcf, csi, sample, region, ref, ref_index, map ]

                VCF_IMPUTE_GLIMPSE1(impute_input)
                output_glimpse1 = VCF_IMPUTE_GLIMPSE1.out.merged_variants
                    .combine(VCF_IMPUTE_GLIMPSE1.out.merged_variants_index, by: 0)
                    .map{ metaIPCR, vcf, csi -> [metaIPCR + [tools: "Glimpse1"], vcf, csi] }
                ch_multiqc_files = ch_multiqc_files.mix(VCF_IMPUTE_GLIMPSE1.out.chunk_chr.map{ [it[1]]})
                ch_versions      = ch_versions.mix(VCF_IMPUTE_GLIMPSE1.out.versions)

                // Add to output channel
                ch_impute_output = ch_impute_output.mix(output_glimpse1)
            }
            if (params.tools.split(',').contains("glimpse2")) {
                error "Glimpse2 not yet implemented"
                // Glimpse2 subworkflow
            }

            if (params.tools.split(',').contains("stitch")) {
                print("Impute with STITCH")

                ch_posfile = []
                // Obtain the user's posfile if provided or calculate it from ref panel file
                if (params.posfile) { // Untested
                    ch_posfile = Channel.of([id:'posfile'], file(params.posfile), checkIfExists:true)
                } else if (params.panel) {
                    // It should do all the panelprep functions if a panel is provided
                    // Currently: the panelprep functions are run by default
                    ch_posfile = VCF_SITES_EXTRACT_BCFTOOLS.out.panel_sites
                } else {
                    error "No posfile or reference panel was included"
                }
                // Prepare inputs
                PREPARE_INPUT_STITCH(ch_posfile, ch_fasta, ch_input_impute)
                ch_versions    = ch_versions.mix(PREPARE_INPUT_STITCH.out.versions)

                // Impute with STITCH
                BAM_IMPUTE_STITCH ( PREPARE_INPUT_STITCH.out.stitch_parameters,
                                    PREPARE_INPUT_STITCH.out.stitch_samples,
                                    ch_fasta )
                ch_versions    = ch_versions.mix(BAM_IMPUTE_STITCH.out.versions)

                // Output channel to concat
                ch_impute_output = ch_impute_output.mix(BAM_IMPUTE_STITCH.out.vcf_tbi)

            }

            if (params.tools.split(',').contains("quilt")) {
                print("Impute with QUILT")

                // Quilt subworkflow

                    // Create chunks from reference VCF
                    MAKE_CHUNKS(ch_panel, ch_fasta)

                    // Impute BAMs with QUILT
                    IMPUTE_QUILT(MAKE_CHUNKS.out.ch_hap_legend, ch_input_impute, MAKE_CHUNKS.out.ch_chunks)
                    ch_versions = ch_versions.mix(IMPUTE_QUILT.out.versions)

                    // Add to output channel
                    ch_impute_output = ch_impute_output.mix(IMPUTE_QUILT.out.vcf_tbi)
            }
            // Concatenate by chromosomes
            CONCAT_IMPUT(ch_impute_output)
            ch_versions       = ch_versions.mix(CONCAT_IMPUT.out.versions)
            ch_input_validate = ch_input_validate.mix(CONCAT_IMPUT.out.vcf_tbi_join)
        }

    }

    if (params.step == 'validate' || params.step == 'all') {
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
            ch_panel_sites_tsv,
            ch_fasta
        )
        ch_multiqc_files = ch_multiqc_files.mix(GL_TRUTH.out.multiqc_files)
        ch_versions      = ch_versions.mix(GL_TRUTH.out.versions)

        // Mix the original vcf and the computed vcf
        ch_truth_vcf = ch_truth.vcf
            .map { [it[0], it[1], it[2]] }
            .mix(GL_TRUTH.out.vcf)

        // Concatenate by chromosomes
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
        ch_versions = ch_versions.mix(VCF_CONCORDANCE_GLIMPSE2.out.versions)
    }

    if (params.step == 'refine') {
        error "refine step not yet implemented"
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
