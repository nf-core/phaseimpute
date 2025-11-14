/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../../modules/nf-core/multiqc'
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { paramsSummaryMultiqc        } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML      } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../../subworkflows/local/utils_nfcore_phaseimpute_pipeline'
include { getFilesSameExt             } from '../../subworkflows/local/utils_nfcore_phaseimpute_pipeline'
include { getFileExtension            } from '../../subworkflows/local/utils_nfcore_phaseimpute_pipeline'
include { exportCsv                   } from '../../subworkflows/local/utils_nfcore_phaseimpute_pipeline'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

// Simulate subworkflows
include { BAM_EXTRACT_REGION_SAMTOOLS                } from '../../subworkflows/local/bam_extract_region_samtools'
include { BAM_DOWNSAMPLE_SAMTOOLS                    } from '../../subworkflows/local/bam_downsample_samtools'
include { SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_INP } from '../../modules/nf-core/samtools/coverage'
include { SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_DWN } from '../../modules/nf-core/samtools/coverage'
include { GAWK as FILTER_CHR_INP                     } from '../../modules/nf-core/gawk'
include { GAWK as FILTER_CHR_DWN                     } from '../../modules/nf-core/gawk'

// Panelprep subworkflows
include { VCF_NORMALIZE_BCFTOOLS                     } from '../../subworkflows/local/vcf_normalize_bcftools'
include { VCF_SITES_EXTRACT_BCFTOOLS                 } from '../../subworkflows/local/vcf_sites_extract_bcftools'
include { VCF_PHASE_SHAPEIT5                         } from '../../subworkflows/local/vcf_phase_shapeit5'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_PANEL   } from '../../subworkflows/local/vcf_concatenate_bcftools'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_PANEL     } from '../../modules/nf-core/bcftools/stats'
include { chunkPrepareChannel                        } from './function.nf'

// Imputation
include { LISTTOFILE                                 } from '../../modules/local/listtofile'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_IMPUTED   } from '../../modules/nf-core/bcftools/query'
include { GAWK as GAWK_IMPUTED                       } from '../../modules/nf-core/gawk'
include { VCF_SPLIT_BCFTOOLS as SPLIT_IMPUTED        } from '../../subworkflows/local/vcf_split_bcftools'

// GLIMPSE1 subworkflows
include { BAM_GL_BCFTOOLS as GL_GLIMPSE1             } from '../../subworkflows/local/bam_gl_bcftools'
include { VCF_IMPUTE_GLIMPSE1                        } from '../../subworkflows/local/vcf_impute_glimpse1'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_GLIMPSE1} from '../../subworkflows/local/vcf_concatenate_bcftools'

// GLIMPSE2 subworkflows
include { BAM_VCF_IMPUTE_GLIMPSE2                    } from '../../subworkflows/local/bam_vcf_impute_glimpse2'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_GLIMPSE2} from '../../subworkflows/local/vcf_concatenate_bcftools'

// QUILT subworkflows
include { VCF_CHUNK_GLIMPSE                          } from '../../subworkflows/local/vcf_chunk_glimpse'
include { BAM_IMPUTE_QUILT                           } from '../../subworkflows/local/bam_impute_quilt'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_QUILT   } from '../../subworkflows/local/vcf_concatenate_bcftools'

// STITCH subworkflows
include { BAM_IMPUTE_STITCH                          } from '../../subworkflows/local/bam_impute_stitch'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_STITCH  } from '../../subworkflows/local/vcf_concatenate_bcftools'

// BEAGLE5 subworkflows
include { VCF_IMPUTE_BEAGLE5                         } from '../../subworkflows/local/vcf_impute_beagle5'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_BEAGLE5 } from '../../subworkflows/local/vcf_concatenate_bcftools'

// MINIMAC4 subworkflows
include { VCF_IMPUTE_MINIMAC4                        } from '../../subworkflows/local/vcf_impute_minimac4'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_MINIMAC4} from '../../subworkflows/local/vcf_concatenate_bcftools'

// Imputation stats
include { BCFTOOLS_STATS as BCFTOOLS_STATS_TOOLS     } from '../../modules/nf-core/bcftools/stats'

// Concordance subworkflows
include { BAM_GL_BCFTOOLS as GL_TRUTH                } from '../../subworkflows/local/bam_gl_bcftools'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_TRUTH     } from '../../modules/nf-core/bcftools/query'
include { GAWK as GAWK_TRUTH                         } from '../../modules/nf-core/gawk'
include { VCF_SPLIT_BCFTOOLS as SPLIT_TRUTH          } from '../../subworkflows/local/vcf_split_bcftools'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_TRUTH     } from '../../modules/nf-core/bcftools/stats'
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
    ch_input_truth          // channel: truth file    [ [id], file, index ]
    ch_fasta                // channel: fasta file    [ [genome], fasta, fai ]
    ch_panel                // channel: panel file    [ [id, chr], vcf, index ]
    ch_region               // channel: region to use [ [chr, region], region]
    ch_depth                // channel: depth select  [ [depth], depth ]
    ch_map                  // channel: genetic map   [ [chr], map]
    ch_posfile              // channel: posfile       [ [id, chr], vcf, index, hap, legend]
    ch_chunks               // channel: chunks        [ [chr], txt]
    chunk_model             // parameter: chunk model
    ch_versions             // channel: versions of software used

    main:

    ch_multiqc_files = Channel.empty()

    //
    // Simulate data if asked
    //
    if (params.steps.split(',').contains("simulate") || params.steps.split(',').contains("all")) {
        // Test if the input are all bam files
        getFilesSameExt(ch_input_sim)
            .map{ if (it != "bam" & it != "cram") {
                error "All input files must be in the same format, either BAM or CRAM, to perform simulation: ${it}"
            } }

        if (params.input_region) {
            // Split the bam into the regions specified
            BAM_EXTRACT_REGION_SAMTOOLS(ch_input_sim, ch_region, ch_fasta)
            ch_versions  = ch_versions.mix(BAM_EXTRACT_REGION_SAMTOOLS.out.versions)
            ch_input_sim = BAM_EXTRACT_REGION_SAMTOOLS.out.bam_region
        }

        // Use input for simulation as truth for validation step
        // if no truth is provided
        if (!params.input_truth) {
            ch_input_truth = ch_input_sim
        }

        // Program to filter chromosomes
        filter_chr_program = ch_region
            .collect{ meta, _region -> meta.chr }
            .map { chr ->
                "BEGIN { FS=\"\t\";\nsplit(\"" + chr.join(" ") + '", chr, " ");\n' +
                'for (i in chr) {\nchr_map[chr[i]] = 1;\n}\n}\n' +
                'NR == 1 || (\$1 in chr_map){\nprint \$0;\n}'
            }
            .collectFile(name:"program.txt")
            .collect()

        // Compute coverage of input files
        SAMTOOLS_COVERAGE_INP(ch_input_sim, ch_fasta)
        ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE_INP.out.versions.first())

        FILTER_CHR_INP(
            SAMTOOLS_COVERAGE_INP.out.coverage,
            filter_chr_program,
            false
        )
        ch_versions = ch_versions.mix(FILTER_CHR_INP.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(FILTER_CHR_INP.out.output.map{ it[1] })

        if (params.depth) {
            // Downsample input to desired depth
            BAM_DOWNSAMPLE_SAMTOOLS(ch_input_sim, ch_depth, ch_fasta)
            ch_versions     = ch_versions.mix(BAM_DOWNSAMPLE_SAMTOOLS.out.versions)
            ch_input_impute = BAM_DOWNSAMPLE_SAMTOOLS.out.bam_emul

            // Compute coverage of input files
            SAMTOOLS_COVERAGE_DWN(BAM_DOWNSAMPLE_SAMTOOLS.out.bam_emul, ch_fasta)
            ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE_DWN.out.versions.first())

            FILTER_CHR_DWN(
                SAMTOOLS_COVERAGE_DWN.out.coverage,
                filter_chr_program,
                false
            )
            ch_versions = ch_versions.mix(FILTER_CHR_DWN.out.versions.first())
            ch_multiqc_files = ch_multiqc_files.mix(FILTER_CHR_DWN.out.output.map{ it[1] })
        }

        if (params.genotype) {
            error "Genotype simulation not yet implemented"
        }

        // Create CSV from simulate step
        exportCsv(
            ch_input_impute.map{ meta, file, index ->
                [meta, [2:"simulation/samples", 3:"simulation/samples"], file, index]
            },
            ["id"], "sample,file,index",
            "simulate.csv", "simulation/csv"
        )
    }

    //
    // Prepare panel
    //
    if (params.steps.split(',').contains("panelprep") || params.steps.split(',').contains("all")) {
        // Normalize indels in panel
        VCF_NORMALIZE_BCFTOOLS(ch_panel, ch_fasta)
        ch_panel_phased = VCF_NORMALIZE_BCFTOOLS.out.vcf_tbi
        ch_versions = ch_versions.mix(VCF_NORMALIZE_BCFTOOLS.out.versions)

        // Extract sites from normalized vcf
        VCF_SITES_EXTRACT_BCFTOOLS(ch_panel_phased, ch_fasta)
        ch_versions = ch_versions.mix(VCF_SITES_EXTRACT_BCFTOOLS.out.versions)

        // Generate all necessary channels
        ch_posfile          = VCF_SITES_EXTRACT_BCFTOOLS.out.posfile

        // Phase panel with Shapeit5
        if (params.phase == true) {
            VCF_PHASE_SHAPEIT5(
                VCF_NORMALIZE_BCFTOOLS.out.vcf_tbi.combine(Channel.of([[]])),
                ch_region,
                [[],[],[]],
                [[],[],[]],
                ch_map,
                chunk_model
            )
            ch_panel_phased = VCF_PHASE_SHAPEIT5.out.vcf_tbi
            ch_versions = ch_versions.mix(VCF_PHASE_SHAPEIT5.out.versions)
        }

        // Create chunks from reference VCF
        VCF_CHUNK_GLIMPSE(ch_panel_phased, ch_map, chunk_model)
        ch_versions = ch_versions.mix(VCF_CHUNK_GLIMPSE.out.versions)

        // Assign chunks channels
        ch_chunks_glimpse1  = VCF_CHUNK_GLIMPSE.out.chunks_glimpse1
        ch_chunks_glimpse2  = VCF_CHUNK_GLIMPSE.out.chunks_glimpse2
        ch_chunks_quilt     = VCF_CHUNK_GLIMPSE.out.chunks_quilt

        // Create CSVs from panelprep step
        // Phased panel
        exportCsv(
            ch_panel_phased.map{ meta, vcf, index ->
                [meta, [2:"prep_panel/panel", 3:"prep_panel/panel"], vcf, index]
            },
            ["id", "chr"], "panel,chr,vcf,index",
            "panel.csv", "prep_panel/csv"
        )
        // Posfile
        exportCsv(
            ch_posfile.map{ meta, vcf, index, hap, legend ->
                [meta, [2:"prep_panel/sites", 3:"prep_panel/sites", 4:"prep_panel/haplegend", 5:"prep_panel/haplegend"], vcf, index, hap, legend]
            },
            ["id", "chr"], "panel,chr,vcf,index,hap,legend",
            "posfile.csv", "prep_panel/csv"
        )
        // Chunks
        exportCsv(
            VCF_CHUNK_GLIMPSE.out.chunks.map{ meta, file ->
                [meta, [2:"prep_panel/chunks/glimpse1"], file]
            },
            ["id", "chr"], "panel,chr,file",
            "chunks_glimpse1.csv", "prep_panel/csv"
        )
    }

    //
    // Impute target files
    //
    if (params.steps.split(',').contains("impute") || params.steps.split(',').contains("all")) {
        // Split input files into BAMs and VCFs
        ch_input_type = ch_input_impute
            .branch {
                bam: it[1] =~ 'bam|cram'
                vcf: it[1] =~ '(vcf|bcf)(.gz)*'
                other: true
            }

        // Check if input files are only BAM/CRAM or VCF/BCF
        ch_input_type.other
            .map{ error "Input files must be either BAM/CRAM or VCF/BCF" }

        // Group BAMs by batch size
        def nb_batch = -1
        ch_input_bams = ch_input_type.bam
            .toSortedList { it1, it2 -> it1[0]["id"] <=> it2[0]["id"] }
            .map { list -> list.collate(params.batch_size)
                .collect{ nb_batch += 1; [[id: "all", batch: nb_batch], it] } }
            .map { list -> [list.collect{ it[0] }, list.collect{ it[1] }] }
            .transpose()
            .map { metaI, filestuples-> [
                metaI + [metas: filestuples.collect{it[0].findAll{it.key != "batch"}}],
                filestuples.collect{it[1]}, filestuples.collect{it[2]}
            ] }

        LISTTOFILE(
            ch_input_bams.map{ meta, file, _index -> [
                meta, file, meta.metas.collect { it.id }
            ] }
        )

        ch_input_bams_withlist = ch_input_bams
            .join(LISTTOFILE.out.txt)

        // Use panel from parameters if provided
        if (params.panel && !params.steps.split(',').find { it in ["all", "panelprep"] }) {
            ch_panel_phased = ch_panel
        }

        if (params.tools.split(',').contains("glimpse1")) {
            log.info("Impute with GLIMPSE1")

            // Use chunks from parameters if provided or use previous chunks from panelprep
            if (params.chunks) {
                ch_chunks_glimpse1 = chunkPrepareChannel(ch_chunks, "glimpse")
            }

            // Glimpse1 subworkflow
            // Compute GL from BAM files and merge them
            GL_GLIMPSE1(
                ch_input_type.bam,
                ch_posfile.map{ [it[0], it[4]] },
                ch_fasta
            )
            ch_multiqc_files = ch_multiqc_files.mix(GL_GLIMPSE1.out.multiqc_files)
            ch_versions = ch_versions.mix(GL_GLIMPSE1.out.versions)

            // Combine vcf and processed bam
            ch_input_glimpse1 = ch_input_type.vcf
                .mix(GL_GLIMPSE1.out.vcf_tbi)

            // Run imputation
            VCF_IMPUTE_GLIMPSE1(
                ch_input_glimpse1,
                ch_panel_phased,
                ch_chunks_glimpse1
            )
            ch_versions = ch_versions.mix(VCF_IMPUTE_GLIMPSE1.out.versions)

            // Concatenate by chromosomes
            CONCAT_GLIMPSE1(VCF_IMPUTE_GLIMPSE1.out.vcf_tbi)
            ch_versions = ch_versions.mix(CONCAT_GLIMPSE1.out.versions)

            // Add results to input validate
            ch_input_validate = ch_input_validate.mix(CONCAT_GLIMPSE1.out.vcf_tbi)

        }

        if (params.tools.split(',').contains("glimpse2")) {
            log.info("Impute with GLIMPSE2")

            if (params.chunks) {
                ch_chunks_glimpse2 = chunkPrepareChannel(ch_chunks, "glimpse")
            }

            // Run imputation
            BAM_VCF_IMPUTE_GLIMPSE2(
                ch_input_bams_withlist
                    .map{ [it[0], it[1], it[2], it[3]] }
                    .mix(ch_input_type.vcf.combine(Channel.of([[]]))),
                ch_panel_phased,
                ch_chunks_glimpse2,
                ch_fasta
            )
            ch_versions = ch_versions.mix(BAM_VCF_IMPUTE_GLIMPSE2.out.versions)
            // Concatenate by chromosomes
            CONCAT_GLIMPSE2(BAM_VCF_IMPUTE_GLIMPSE2.out.vcf_tbi)
            ch_versions = ch_versions.mix(CONCAT_GLIMPSE2.out.versions)

            // Add results to input validate
            ch_input_validate = ch_input_validate.mix(CONCAT_GLIMPSE2.out.vcf_tbi)
        }

        if (params.tools.split(',').contains("stitch")) {
            log.info("Impute with STITCH")

            // Impute with STITCH
            BAM_IMPUTE_STITCH (
                ch_input_bams_withlist.map{ [it[0], it[1], it[2], it[4], it[5]] },
                ch_posfile.map{ [it[0], it[4]] },
                ch_region,
                ch_fasta
            )
            ch_versions = ch_versions.mix(BAM_IMPUTE_STITCH.out.versions)

            // Concatenate by chromosomes
            CONCAT_STITCH(BAM_IMPUTE_STITCH.out.vcf_tbi)
            ch_versions = ch_versions.mix(CONCAT_STITCH.out.versions)

            // Add results to input validate
            ch_input_validate = ch_input_validate.mix(CONCAT_STITCH.out.vcf_tbi)

        }

        if (params.tools.split(',').contains("quilt")) {
            log.info("Impute with QUILT")

            // Use provided chunks if --chunks
            if (params.chunks) {
                ch_chunks_quilt = chunkPrepareChannel(ch_chunks, "quilt")
            }

            // Impute BAMs with QUILT
            BAM_IMPUTE_QUILT(
                ch_input_bams_withlist.map{ [it[0], it[1], it[2], it[4], it[5]] },
                ch_posfile.map{ [it[0], it[3], it[4]] },
                ch_chunks_quilt,
                ch_fasta.map{ [it[0], it[1]] }
            )
            ch_versions = ch_versions.mix(BAM_IMPUTE_QUILT.out.versions)

            // Concatenate by chromosomes
            CONCAT_QUILT(BAM_IMPUTE_QUILT.out.vcf_tbi)
            ch_versions = ch_versions.mix(CONCAT_QUILT.out.versions)

            // Add results to input validate
            ch_input_validate = ch_input_validate.mix(CONCAT_QUILT.out.vcf_tbi)
        }

        if (params.tools.split(',').contains("beagle5")) {
            // Create input channel combining VCF with regions
            ch_input_beagle5 = ch_input_type.vcf
                .combine(ch_region)
                .map { meta_vcf, vcf, index, meta_region, _region ->
                    [meta_vcf + meta_region, vcf, index]
                }

            // Impute with BEAGLE5
            VCF_IMPUTE_BEAGLE5(
                ch_input_beagle5,
                ch_panel_phased,
                ch_map
            )
            ch_versions = ch_versions.mix(VCF_IMPUTE_BEAGLE5.out.versions)

            // Concatenate by chromosomes
            CONCAT_BEAGLE5(VCF_IMPUTE_BEAGLE5.out.vcf_index)
            ch_versions = ch_versions.mix(CONCAT_BEAGLE5.out.versions)

            // Add results to input validate
            ch_input_validate = ch_input_validate.mix(CONCAT_BEAGLE5.out.vcf_tbi)
        }

        if (params.tools.split(',').contains("minimac4")) {
            log.info("Impute with MINIMAC4")

            // Create input channel combining VCF with regions
            ch_input_minimac4 = ch_input_type.vcf
                .combine(ch_region)
                .map { meta_vcf, vcf, index, meta_region, region ->
                    [meta_vcf + meta_region, vcf, index]
                }

            // Run imputation with MINIMAC4
            VCF_IMPUTE_MINIMAC4(
                ch_input_minimac4,
                ch_panel_phased,
                ch_map,
                ch_posfile
            )
            ch_versions = ch_versions.mix(VCF_IMPUTE_MINIMAC4.out.versions)

            // Concatenate by chromosomes
            CONCAT_MINIMAC4(VCF_IMPUTE_MINIMAC4.out.vcf_index)
            ch_versions = ch_versions.mix(CONCAT_MINIMAC4.out.versions)

            // Add results to input validate
            ch_input_validate = ch_input_validate.mix(CONCAT_MINIMAC4.out.vcf_tbi)
        }

        // Prepare renaming file
        BCFTOOLS_QUERY_IMPUTED(ch_input_validate, [], [], [])
        GAWK_IMPUTED(BCFTOOLS_QUERY_IMPUTED.out.output, [], false)
        ch_split_imputed = ch_input_validate.join(GAWK_IMPUTED.out.output)

        // Split result by samples
        SPLIT_IMPUTED(ch_split_imputed)
        ch_versions = ch_versions.mix(SPLIT_IMPUTED.out.versions)
        ch_input_validate = SPLIT_IMPUTED.out.vcf_tbi

        // Compute stats on imputed files
        BCFTOOLS_STATS_TOOLS(
            ch_input_validate,
            [[],[]],
            [[],[]],
            [[],[]],
            [[],[]],
            ch_fasta.map{ [it[0], it[1]] })
        ch_versions = ch_versions.mix(BCFTOOLS_STATS_TOOLS.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_STATS_TOOLS.out.stats.map{ [it[1]] })

        // Export all files to csv
        exportCsv(
            ch_input_validate.map{ meta, file, index ->
                [meta, [2:"imputation/${meta.tools}/samples", 3:"imputation/${meta.tools}/samples"], file, index]
            },
            ["id", "tools"], "sample,tools,file,index",
            "impute.csv", "imputation/csv"
        )
    }

    if (params.steps.split(',').contains("validate") || params.steps.split(',').contains("all")) {
        // Concatenate all sites into a single VCF (for GLIMPSE concordance)
        CONCAT_PANEL(ch_posfile.map{ [it[0], it[1], it[2]] })
        ch_versions    = ch_versions.mix(CONCAT_PANEL.out.versions)
        ch_panel_sites = CONCAT_PANEL.out.vcf_tbi

        // Compute stats on panel
        BCFTOOLS_STATS_PANEL(
            ch_panel_sites,
            [[],[]],
            [[],[]],
            [[],[]],
            [[],[]],
            ch_fasta.map{ [it[0], it[1]] })
        ch_versions = ch_versions.mix(BCFTOOLS_STATS_PANEL.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_STATS_PANEL.out.stats.map{ [it[1]] })

        ch_truth_vcf = Channel.empty()

        // Channels for branching
        ch_truth = ch_input_truth
            .map { [it[0], it[1], it[2], getFileExtension(it[1])] }
            .branch {
                bam: it[3] =~ 'bam|cram'
                vcf: it[3] =~ '(vcf|bcf)(.gz)*'
                other: true
            }

        ch_truth.other
            .map{ error "Input files must be either BAM/CRAM or VCF/BCF" }

        GL_TRUTH(
            ch_truth.bam.map { [it[0], it[1], it[2]] },
            ch_posfile.map{ [it[0], it[4]] },
            ch_fasta
        )
        ch_versions = ch_versions.mix(GL_TRUTH.out.versions)

        // Mix the original vcf and the computed vcf
        ch_truth_vcf = ch_truth.vcf
            .map { [it[0], it[1], it[2]] }
            .mix(GL_TRUTH.out.vcf_tbi)

        // Concatenate truth vcf by chromosomes
        CONCAT_TRUTH(ch_truth_vcf)
        ch_versions = ch_versions.mix(CONCAT_TRUTH.out.versions)

        // Prepare renaming file
        BCFTOOLS_QUERY_TRUTH(CONCAT_TRUTH.out.vcf_tbi, [], [], [])
        GAWK_TRUTH(BCFTOOLS_QUERY_TRUTH.out.output, [], false)
        ch_split_truth = CONCAT_TRUTH.out.vcf_tbi.join(GAWK_TRUTH.out.output)

        // Split truth vcf by samples
        SPLIT_TRUTH(ch_split_truth)
        ch_versions = ch_versions.mix(SPLIT_TRUTH.out.versions)

        // Compute stats on truth files
        BCFTOOLS_STATS_TRUTH(
            SPLIT_TRUTH.out.vcf_tbi,
            [[],[]],
            [[],[]],
            [[],[]],
            [[],[]],
            ch_fasta.map{ [it[0], it[1]] }
        )
        ch_versions = ch_versions.mix(BCFTOOLS_STATS_TRUTH.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_STATS_TRUTH.out.stats.map{ [it[1]] })

        // Compute concordance analysis
        VCF_CONCORDANCE_GLIMPSE2(
            ch_input_validate,
            SPLIT_TRUTH.out.vcf_tbi,
            ch_panel_sites,
            ch_region
        )
        ch_multiqc_files = ch_multiqc_files.mix(VCF_CONCORDANCE_GLIMPSE2.out.multiqc_files)
        ch_versions      = ch_versions.mix(VCF_CONCORDANCE_GLIMPSE2.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true, newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))
    ch_multiqc_replace_names              = params.multiqc_replace_names ? Channel.fromPath(params.multiqc_replace_names, checkIfExists: true) : Channel.empty()
    ch_multiqc_sample_names               = params.multiqc_sample_names ? Channel.fromPath(params.multiqc_sample_names, checkIfExists: true) : Channel.empty()

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        ch_multiqc_replace_names.toList(),
        ch_multiqc_sample_names.toList()
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
