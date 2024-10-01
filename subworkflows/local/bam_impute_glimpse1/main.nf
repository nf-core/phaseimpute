
include { BAM_GL_BCFTOOLS                    } from '../bam_gl_bcftools'
include { GLIMPSE_PHASE                      } from '../../../modules/nf-core/glimpse/phase'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_1 } from '../../../modules/nf-core/bcftools/index'
include { GLIMPSE_LIGATE                     } from '../../../modules/nf-core/glimpse/ligate'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2 } from '../../../modules/nf-core/bcftools/index'

workflow BAM_IMPUTE_GLIMPSE1 {

    take:
    ch_input        // channel (mandatory): [ [id], bam, bai ]
    ch_posfile      // channel (mandatory): [ [panel, chr], legend ]
    ch_panel        // channel (mandatory): [ [panel, chr], vcf, tbi ]
    ch_chunks       // channel  (optional): [ [chr], region1, region2 ]
    ch_fasta        // channel (mandatory): [ [genome], fa, fai ]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Channels for branching
    ch_input = ch_input
        .branch {
            bam: it[1] =~ 'bam|cram'
            vcf: it[1] =~ '(vcf|bcf)(.gz)*'
            other: true
        }
    
    ch_input.other
        .map{ error "Input files must be either BAM/CRAM or VCF/BCF" }

    // Glimpse1 subworkflow
    BAM_GL_BCFTOOLS( // Compute GL for input data once per panel by chromosome
        ch_input.bam,
        ch_posfile,
        ch_fasta
    )
    ch_multiqc_files = ch_multiqc_files.mix(BAM_GL_BCFTOOLS.out.multiqc_files)
    ch_versions = ch_versions.mix(BAM_GL_BCFTOOLS.out.versions)

    ch_impute = ch_input.vcf
        .combine(ch_posfile).view()
        .map{ metaI, vcf, index, metaPC, legend ->
            [metaI + ["panel": metaPC.id, "chr": metaPC.chr], vcf, index]
        }
        .mix(BAM_GL_BCFTOOLS.out.vcf_tbi)

    samples_file = Channel.of([[]]).collect()
    gmap_file    = Channel.of([[]]).collect()

    // Combine chunks with panel
    ch_chunks_panel = ch_chunks
        .combine(ch_panel, by:0)
        .map{ metaPC, regionin, regionout, panel, index ->
            [["panel": metaPC.id, "chr": metaPC.chr], regionin, regionout, panel, index]
        }

    // Join input and chunks reference
    ch_phase_input = ch_impute
        .map{ metaIPC, vcf, index -> [metaIPC.subMap("panel", "chr"), metaIPC, vcf, index] }
        .combine(samples_file).view()
        .combine(ch_chunks_panel, by: 0)
        .combine(gmap_file)
        .map{ metaPC, metaIPC, bam, bai, samples, regionin, regionout, panel, panel_index, gmap ->
            [metaIPC + ["chunk": regionout],
            bam, bai, samples, regionin, regionout, panel, panel_index, gmap]
        }.view()

    GLIMPSE_PHASE ( ch_phase_input ) // [meta, vcf, index, sample, regionin, regionout, ref, ref_index, map]
    ch_versions = ch_versions.mix(GLIMPSE_PHASE.out.versions )

    BCFTOOLS_INDEX_1 ( GLIMPSE_PHASE.out.phased_variants )
    ch_versions = ch_versions.mix( BCFTOOLS_INDEX_1.out.versions )

    // Ligate all phased files in one and index it
    ligate_input = GLIMPSE_PHASE.out.phased_variants
        .join( BCFTOOLS_INDEX_1.out.csi )
        .map{ metaIPCR, vcf, index -> [metaIPCR.subMap("id", "panel", "chr"), vcf, index] }
        .groupTuple()

    GLIMPSE_LIGATE ( ligate_input )
    ch_versions = ch_versions.mix(GLIMPSE_LIGATE.out.versions )

    BCFTOOLS_INDEX_2 ( GLIMPSE_LIGATE.out.merged_variants )
    ch_versions = ch_versions.mix( BCFTOOLS_INDEX_2.out.versions )

    // Join imputed and index files
    ch_imputed_vcf_tbi = GLIMPSE_LIGATE.out.merged_variants
        .join(BCFTOOLS_INDEX_2.out.tbi)
        .map{ metaIPC, vcf, index -> [metaIPC + [tools: "glimpse1"], vcf, index] }

    emit:
    vcf_tbi             = ch_imputed_vcf_tbi    // channel: [ [id, panel, chr, tool], vcf, tbi ]
    versions            = ch_versions           // channel: [ versions.yml ]
    multiqc_files       = ch_multiqc_files      // channel: [ multiqc_files.yml ]
}
