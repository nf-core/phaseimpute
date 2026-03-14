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
include { BAM_SUBSAMPLEDEPTH_SAMTOOLS                } from '../../subworkflows/nf-core/bam_subsampledepth_samtools'
include { SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_INP } from '../../modules/nf-core/samtools/coverage'
include { SAMTOOLS_COVERAGE as SAMTOOLS_COVERAGE_DWN } from '../../modules/nf-core/samtools/coverage'
include { GAWK as FILTER_CHR_INP                     } from '../../modules/nf-core/gawk'
include { GAWK as FILTER_CHR_DWN                     } from '../../modules/nf-core/gawk'

// Panelprep subworkflows
include { VCF_NORMALIZE_BCFTOOLS                     } from '../../subworkflows/local/vcf_normalize_bcftools'
include { VCF_SITES_EXTRACT_BCFTOOLS                 } from '../../subworkflows/local/vcf_sites_extract_bcftools'
include { VCF_PHASE_SHAPEIT5                         } from '../../subworkflows/nf-core/vcf_phase_shapeit5'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_PANEL   } from '../../subworkflows/local/vcf_concatenate_bcftools'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_PANEL     } from '../../modules/nf-core/bcftools/stats'
include { VCF_CHUNK_GLIMPSE                          } from '../../subworkflows/local/vcf_chunk_glimpse'
include { chunkPrepareChannel                        } from './function.nf'

// Imputation
include { LISTTOFILE                                 } from '../../modules/local/listtofile'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_IMPUTED   } from '../../modules/nf-core/bcftools/query'
include { GAWK as GAWK_IMPUTED                       } from '../../modules/nf-core/gawk'
include { VCF_SPLIT_BCFTOOLS as SPLIT_IMPUTED        } from '../../subworkflows/local/vcf_split_bcftools'

// GLIMPSE1 subworkflows
include { BAM_VARIANT_CALLING_MPILEUP_BCFTOOLS as GL_GLIMPSE1 } from '../../subworkflows/nf-core/bam_variant_calling_mpileup_bcftools'
include { VCF_IMPUTE_GLIMPSE                                  } from '../../subworkflows/nf-core/vcf_impute_glimpse'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_GLIMPSE1         } from '../../subworkflows/local/vcf_concatenate_bcftools'

// GLIMPSE2 subworkflows
include { BAM_VCF_IMPUTE_GLIMPSE2                    } from '../../subworkflows/nf-core/bam_vcf_impute_glimpse2'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_GLIMPSE2} from '../../subworkflows/local/vcf_concatenate_bcftools'

// QUILT subworkflows
include { GAWK as GAWK_POSFILE_QUILT                 } from '../../modules/nf-core/gawk'
include { TABIX_BGZIP as BGZIP_POSFILE_QUILT         } from '../../modules/nf-core/tabix/bgzip'
include { BAM_IMPUTE_QUILT                           } from '../../subworkflows/nf-core/bam_impute_quilt'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_QUILT   } from '../../subworkflows/local/vcf_concatenate_bcftools'

// STITCH subworkflows
include { GAWK as GAWK_POSFILE_STITCH                } from '../../modules/nf-core/gawk'
include { TABIX_BGZIP as BGZIP_POSFILE_STITCH        } from '../../modules/nf-core/tabix/bgzip'
include { BAM_IMPUTE_STITCH                          } from '../../subworkflows/nf-core/bam_impute_stitch'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_STITCH  } from '../../subworkflows/local/vcf_concatenate_bcftools'

// BEAGLE5 subworkflows
include { VCF_IMPUTE_BEAGLE5                         } from '../../subworkflows/nf-core/vcf_impute_beagle5'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_BEAGLE5 } from '../../subworkflows/local/vcf_concatenate_bcftools'

// MINIMAC4 subworkflows
include { VCF_IMPUTE_MINIMAC4                        } from '../../subworkflows/nf-core/vcf_impute_minimac4'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_MINIMAC4} from '../../subworkflows/local/vcf_concatenate_bcftools'

// Imputation stats
include { BCFTOOLS_STATS as BCFTOOLS_STATS_TOOLS     } from '../../modules/nf-core/bcftools/stats'

// Concordance subworkflows
include { BAM_VARIANT_CALLING_MPILEUP_BCFTOOLS as GL_TRUTH } from '../../subworkflows/nf-core/bam_variant_calling_mpileup_bcftools'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_TRUTH           } from '../../modules/nf-core/bcftools/query'
include { GAWK as GAWK_TRUTH                               } from '../../modules/nf-core/gawk'
include { VCF_SPLIT_BCFTOOLS as SPLIT_TRUTH                } from '../../subworkflows/local/vcf_split_bcftools'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_TRUTH           } from '../../modules/nf-core/bcftools/stats'
include { VCF_CONCATENATE_BCFTOOLS as CONCAT_TRUTH         } from '../../subworkflows/local/vcf_concatenate_bcftools'
include { VCF_CONCORDANCE_GLIMPSE2                         } from '../../subworkflows/local/vcf_concordance_glimpse2'


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
    ch_posfile              // channel: posfile       [ [id, chr], vcf, index, hap, legend, posfile]
    ch_chunks               // channel: chunks        [ [chr], txt]
    chunk_model             // parameter: chunk model

    main:

    ch_multiqc_files = channel.empty()

    //
    // Simulate data if asked
    //
    if (params.steps.split(',').contains("simulate") || params.steps.split(',').contains("all")) {
        // Test if the input are all bam files
        getFilesSameExt(ch_input_sim)
            .map{ ext -> if (ext != "bam" & ext != "cram") {
                error "All input files must be in the same format, either BAM or CRAM, to perform simulation: ${ext}"
            } }

        if (params.input_region) {
            // Split the bam into the regions specified
            BAM_EXTRACT_REGION_SAMTOOLS(ch_input_sim, ch_region, ch_fasta)
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

        FILTER_CHR_INP(
            SAMTOOLS_COVERAGE_INP.out.coverage,
            filter_chr_program,
            false
        )
        ch_multiqc_files = ch_multiqc_files.mix(FILTER_CHR_INP.out.output.map{ _meta, file -> file })

        if (params.depth) {
            // Downsample input to desired depth
            BAM_SUBSAMPLEDEPTH_SAMTOOLS(ch_input_sim, ch_depth, ch_fasta)
            ch_input_impute = BAM_SUBSAMPLEDEPTH_SAMTOOLS.out.bam_subsampled
                .map{ meta, bam, index ->
                    def keysToKeep = meta.keySet() - ['subsample_fraction']
                    def newMeta = meta.subMap(keysToKeep)
                    [ newMeta, bam, index ]
                }

            // Compute coverage of input files
            SAMTOOLS_COVERAGE_DWN(ch_input_impute, ch_fasta)

            FILTER_CHR_DWN(
                SAMTOOLS_COVERAGE_DWN.out.coverage,
                filter_chr_program,
                false
            )
            ch_multiqc_files = ch_multiqc_files.mix(FILTER_CHR_DWN.out.output.map{ _meta, file -> file })
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

        // Extract sites from normalized vcf
        VCF_SITES_EXTRACT_BCFTOOLS(ch_panel_phased, ch_fasta)

        // Generate all necessary channels
        if (!params.posfile){
            ch_posfile  = VCF_SITES_EXTRACT_BCFTOOLS.out.posfile
        }

        // Use glimpse 1 for chunks if not provided
        if (!params.chunks){
            // Create chunks from reference VCF
            VCF_CHUNK_GLIMPSE(
                VCF_NORMALIZE_BCFTOOLS.out.vcf_tbi,
                ch_map,
                chunk_model
            )
            ch_chunks  = VCF_CHUNK_GLIMPSE.out.chunks

            // Chunks
            exportCsv(
                ch_chunks
                .map{ meta, file ->
                    [meta, [2:"prep_panel/chunks/glimpse1"], file]
                },
                ["panel_id", "chr"], "panel,chr,file",
                "chunks_glimpse1.csv", "prep_panel/csv"
            )
        }

        // Phase panel with Shapeit5
        if (params.phase == true) {
            // Use chunks from parameters and use region with buffer region
            ch_chunks_phase = chunkPrepareChannel(ch_chunks, ch_region, "glimpse1")

            VCF_PHASE_SHAPEIT5(
                VCF_NORMALIZE_BCFTOOLS.out.vcf_tbi.combine(channel.of([[], []])), // No pedigree, no region
                ch_chunks_phase.map{ meta, _regionin, regionout -> [meta, regionout]},
                VCF_NORMALIZE_BCFTOOLS.out.vcf_tbi.map{ meta, _vcf, _index -> [meta, [], []]}, // No ref
                VCF_NORMALIZE_BCFTOOLS.out.vcf_tbi.map{ meta, vcf, index -> [meta, [], []]}, // No scaffold
                ch_map,
                false,
                chunk_model
            )
            ch_panel_phased = VCF_PHASE_SHAPEIT5.out.vcf_index
        }

        // Create CSVs from panelprep step
        // Phased panel
        exportCsv(
            ch_panel_phased.map{ meta, vcf, index ->
                [meta, [2:"prep_panel/panel", 3:"prep_panel/panel"], vcf, index]
            },
            ["panel_id", "chr"], "panel,chr,vcf,index",
            "panel.csv", "prep_panel/csv"
        )
        // Posfile
        exportCsv(
            ch_posfile.map{ meta, vcf, index, hap, legend, posfile ->
                [
                    meta, [2:"prep_panel/sites", 3:"prep_panel/sites", 4:"prep_panel/haplegend", 5:"prep_panel/haplegend", 6:"prep_panel/posfile"],
                    vcf, index, hap, legend, posfile
                ]
            },
            ["panel_id", "chr"], "panel,chr,vcf,index,hap,legend,posfile",
            "posfile.csv", "prep_panel/csv"
        )
    }

    //
    // Impute target files
    //
    if (params.steps.split(',').contains("impute") || params.steps.split(',').contains("all")) {
        // Split input files into BAMs and VCFs
        ch_input_type = ch_input_impute
            .branch { _meta, file, _index ->
                bam: file =~ 'bam|cram'
                vcf: file =~ '(vcf|bcf)(.gz)*'
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
                .collect{ it -> nb_batch += 1; [
                    [id: "all_samples", batch: nb_batch], it]
                }
            }
            .map { list -> [
                list.collect{ it -> it[0] },
                list.collect{ it -> it[1] }
            ] }
            .transpose()
            .map { metaI, filestuples-> [
                metaI + [metas: filestuples.collect{ meta, _file, _index -> meta.findAll{keys -> keys.key != "batch"}}],
                filestuples.collect{_meta, file, _index -> file},
                filestuples.collect{_meta, _file, index -> index}
            ] }

        LISTTOFILE(
            ch_input_bams.map{ meta, file, _index -> [
                meta, file, meta.metas.collect {meta_i -> meta_i.id }
            ] }
        )

        ch_input_bams_withlist = ch_input_bams
            .join(LISTTOFILE.out.txt)

        // Use panel from parameters if provided
        if (params.panel && !params.steps.split(',').find { step -> step in ["all", "panelprep"] }) {
            ch_panel_phased = ch_panel
        }

        if (params.tools.split(',').contains("glimpse1")) {
            log.info("Impute with GLIMPSE1")

            // Use chunks from parameters if provided or use previous chunks from panelprep
            ch_chunks_glimpse1 = chunkPrepareChannel(ch_chunks, ch_region, "glimpse1")

            // Glimpse1 subworkflow
            // Compute GL from BAM files and merge them
            GL_GLIMPSE1(
                ch_input_type.bam,
                ch_posfile.map{
                    meta, _site, _site_index, _hap, _legend, posfile -> [
                        meta, posfile
                    ]
                },
                ch_fasta,
                "id",
                "all_samples",
                [ "panel_id", "id", "batch", "tools" ],
                false,
                true
            )
            ch_multiqc_files = ch_multiqc_files.mix(GL_GLIMPSE1.out.multiqc_files)

            // Combine vcf and processed bam
            ch_input_glimpse1 = ch_input_type.vcf
                .mix(GL_GLIMPSE1.out.vcf_index)
                .map{
                    meta, vcf, index -> [
                        meta, vcf, index, [] // Ignore infos for the moment
                    ]
                }

            // Run imputation
            VCF_IMPUTE_GLIMPSE(
                ch_input_glimpse1,
                ch_panel_phased.map{
                    meta, file, index ->
                    [meta, file, index, []] // Region ignored as chunks are provided
                },
                ch_chunks_glimpse1,
                ch_map,
                false // Do not compute chunks
            )

            // Concatenate by chromosomes
            CONCAT_GLIMPSE1(VCF_IMPUTE_GLIMPSE.out.vcf_index.map{
                meta, vcf, index -> [meta + [tools:"glimpse1"], vcf, index]
            })

            // Add results to input validate
            ch_input_validate = ch_input_validate.mix(CONCAT_GLIMPSE1.out.vcf_index)

        }

        if (params.tools.split(',').contains("glimpse2")) {
            log.info("Impute with GLIMPSE2")

            ch_chunks_glimpse2 = chunkPrepareChannel(ch_chunks, ch_region, "glimpse1")

            // Run imputation
            BAM_VCF_IMPUTE_GLIMPSE2(
                ch_input_bams_withlist.map{
                    meta, file, index, bampath_id, _bampath_noid, _bamnames->
                    [meta, file, index, bampath_id, []]
                },
                ch_panel_phased.map{
                    meta, file, index ->
                    [meta, file, index, []] // Region ignored as chunks are provided
                },
                ch_chunks_glimpse2,
                ch_map,
                ch_fasta,
                false, "", false
            )

            // Concatenate by chromosomes
            CONCAT_GLIMPSE2(BAM_VCF_IMPUTE_GLIMPSE2.out.vcf_index.map{
                meta, vcf, index -> [meta + [tools:"glimpse2"], vcf, index]
            })

            // Add results to input validate
            ch_input_validate = ch_input_validate.mix(CONCAT_GLIMPSE2.out.vcf_index)
        }

        if (params.tools.split(',').contains("stitch")) {
            log.info("Impute with STITCH")

            ch_chunks_stitch = chunkPrepareChannel(ch_chunks, ch_region, "quilt")

            // Transform posfile to tabulated format
            GAWK_POSFILE_STITCH(
                ch_posfile.map{
                    meta, _site, _site_index, _hap, _legend, posfile -> [
                        meta, posfile
                    ]
                }, [], false)

            BGZIP_POSFILE_STITCH(GAWK_POSFILE_STITCH.out.output)

            // Impute with STITCH
            BAM_IMPUTE_STITCH (
                ch_input_bams_withlist.map{
                    meta, file, index, _bampath_id, bampath_noid, bamnames -> [
                        meta, file, index, bampath_noid, bamnames
                    ]
                },
                BGZIP_POSFILE_STITCH.out.output,
                ch_chunks_stitch,
                ch_map,
                ch_fasta,
                params.k_val,
                params.ngen,
                params.seed
            )

            // Concatenate by chromosomes
            CONCAT_STITCH(BAM_IMPUTE_STITCH.out.vcf_index.map{
                meta, vcf, index -> [meta + [tools:"stitch"], vcf, index]
            })

            // Add results to input validate
            ch_input_validate = ch_input_validate.mix(CONCAT_STITCH.out.vcf_index)

        }

        if (params.tools.split(',').contains("quilt")) {
            log.info("Impute with QUILT")

            // Use provided chunks if --chunks or whole chromosome
            ch_chunks_quilt = chunkPrepareChannel(ch_chunks, ch_region, "quilt")

            // Transform posfile to tabulated format
            GAWK_POSFILE_QUILT(
                ch_posfile.map{
                    meta, _site, _site_index, _hap, _legend, posfile -> [
                        meta, posfile
                    ]
                }, [], false
            )

            BGZIP_POSFILE_QUILT(GAWK_POSFILE_QUILT.out.output)

            ch_posfile_quilt = ch_posfile
                .map{
                    meta, _site, _site_index, hap, legend, _posfile -> [
                        meta, hap, legend
                    ]
                }.join(BGZIP_POSFILE_QUILT.out.output)

            // Impute BAMs with QUILT
            BAM_IMPUTE_QUILT(
                ch_input_bams_withlist.map{
                    meta, file, index, _bampath_id, bampath_noid, bamnames -> [
                        meta, file, index, bampath_noid, bamnames
                    ]
                },
                ch_posfile_quilt,
                ch_chunks_quilt,
                ch_map,
                ch_fasta,
                params.ngen,
                params.buffer
            )

            // Concatenate by chromosomes
            CONCAT_QUILT(BAM_IMPUTE_QUILT.out.vcf_index
                .map{
                    meta, vcf, index -> [meta + [tools:"quilt"], vcf, index]
                }
            )

            // Add results to input validate
            ch_input_validate = ch_input_validate.mix(CONCAT_QUILT.out.vcf_index)
        }

        if (params.tools.split(',').contains("beagle5")) {
            log.info("Impute with BEAGLE5")
            ch_chunks_beagle5 = chunkPrepareChannel(ch_chunks, ch_region, "glimpse1")
                .map{ meta, _regionin, regionout -> [meta, regionout]}

            // Impute with BEAGLE5
            VCF_IMPUTE_BEAGLE5(
                ch_input_type.vcf,
                ch_panel_phased,
                ch_chunks_beagle5,
                ch_map
            )

            // Concatenate by chromosomes
            CONCAT_BEAGLE5(VCF_IMPUTE_BEAGLE5.out.vcf_index.map{
                meta, vcf, index -> [meta + [tools:"beagle5"], vcf, index]
            })

            // Add results to input validate
            ch_input_validate = ch_input_validate.mix(CONCAT_BEAGLE5.out.vcf_index)
        }

        if (params.tools.split(',').contains("minimac4")) {
            log.info("Impute with MINIMAC4")

            ch_chunks_minimac4 = chunkPrepareChannel(ch_chunks, ch_region, "glimpse1")
                .map{ meta, _regionin, regionout -> [meta, regionout]}

            // Create input channel combining VCF with regions
            ch_input_minimac4 = ch_input_type.vcf
                .combine(ch_region)
                .map { meta_vcf, vcf, index, meta_region, _region ->
                    [meta_vcf + meta_region, vcf, index]
                }

            // Run imputation with MINIMAC4
            VCF_IMPUTE_MINIMAC4(
                ch_input_minimac4,
                ch_panel_phased,
                ch_posfile.map{
                    meta, site, site_index, _hap, _legend, _posfile -> [
                        meta, site, site_index
                    ]
                },
                ch_chunks_minimac4,
                ch_map
            )

            // Concatenate by chromosomes
            CONCAT_MINIMAC4(VCF_IMPUTE_MINIMAC4.out.vcf_index.map{
                meta, vcf, index -> [meta + [tools:"minimac4"], vcf, index]
            })

            // Add results to input validate
            ch_input_validate = ch_input_validate.mix(CONCAT_MINIMAC4.out.vcf_index)
        }

        // Prepare renaming file
        BCFTOOLS_QUERY_IMPUTED(ch_input_validate, [], [], [])
        GAWK_IMPUTED(BCFTOOLS_QUERY_IMPUTED.out.output, [], false)
        ch_split_imputed = ch_input_validate.join(GAWK_IMPUTED.out.output)

        // Split result by samples
        SPLIT_IMPUTED(ch_split_imputed)
        ch_input_validate = SPLIT_IMPUTED.out.vcf_tbi

        // Compute stats on imputed files
        BCFTOOLS_STATS_TOOLS(
            ch_input_validate,
            [[],[]],
            [[],[]],
            [[],[]],
            [[],[]],
            ch_fasta.map{ meta, fasta, _index -> [ meta, fasta ] }
        )
        ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_STATS_TOOLS.out.stats.map{ _meta, file -> [ file ] })

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
        CONCAT_PANEL(ch_posfile.map{
            meta, site, site_index, _hap, _legend, _posfile -> [
                meta, site, site_index
            ]
        })
        ch_panel_sites = CONCAT_PANEL.out.vcf_index

        // Compute stats on panel
        BCFTOOLS_STATS_PANEL(
            ch_panel_sites,
            [[],[]],
            [[],[]],
            [[],[]],
            [[],[]],
            ch_fasta.map{ meta, fasta, _index -> [meta, fasta] }
        )
        ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_STATS_PANEL.out.stats.map{ _meta, file -> [ file ] })

        ch_truth_vcf = channel.empty()

        // Channels for branching
        ch_truth = ch_input_truth
            .map { meta, file, index -> [meta, file, index, getFileExtension(file)] }
            .branch { _meta, _file, _index, ext ->
                bam: ext =~ 'bam|cram'
                vcf: ext =~ '(vcf|bcf)(.gz)*'
                other: true
            }

        ch_truth.other
            .subscribe { error "Input files must be either BAM/CRAM or VCF/BCF" }

        GL_TRUTH(
            ch_truth.bam.map { meta, file, index, _ext -> [meta, file, index] },
            ch_posfile.map{
                meta, _site, _site_index, _hap, _legend, posfile -> [
                    meta, posfile
                ]
            },
            ch_fasta,
            "id",
            "all_samples",
            [ "panel_id", "id" ],
            false,
            true
        )

        // Mix the original vcf and the computed vcf
        ch_truth_vcf = ch_truth.vcf
            .map { meta, file, index, _ext -> [meta, file, index] }
            .mix(GL_TRUTH.out.vcf_index)

        // Prepare renaming file
        BCFTOOLS_QUERY_TRUTH(ch_truth_vcf, [], [], [])
        GAWK_TRUTH(BCFTOOLS_QUERY_TRUTH.out.output, [], false)
        ch_split_truth = ch_truth_vcf.join(GAWK_TRUTH.out.output)

        // Split truth vcf by samples
        SPLIT_TRUTH(ch_split_truth)

        // Compute stats on truth files
        BCFTOOLS_STATS_TRUTH(
            SPLIT_TRUTH.out.vcf_tbi,
            [[],[]],
            [[],[]],
            [[],[]],
            [[],[]],
            ch_fasta.map{ meta, fasta, _index -> [meta, fasta] }
        )
        ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_STATS_TRUTH.out.stats.map{ _meta, file -> [ file ] })

        // Compute concordance analysis
        VCF_CONCORDANCE_GLIMPSE2(
            ch_input_validate,
            SPLIT_TRUTH.out.vcf_tbi,
            ch_panel_sites,
            ch_region
        )
        ch_multiqc_files = ch_multiqc_files.mix(VCF_CONCORDANCE_GLIMPSE2.out.multiqc_files)
    }

    //
    // Collate and save software versions
    //

    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'phaseimpute_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true) : channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? channel.fromPath(params.multiqc_logo, checkIfExists: true) : channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))
    ch_multiqc_replace_names              = params.multiqc_replace_names ? channel.fromPath(params.multiqc_replace_names, checkIfExists: true) : channel.empty()
    ch_multiqc_sample_names               = params.multiqc_sample_names ? channel.fromPath(params.multiqc_sample_names, checkIfExists: true) : channel.empty()

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
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
