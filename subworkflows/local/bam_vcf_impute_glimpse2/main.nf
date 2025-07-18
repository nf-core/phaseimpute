include { GLIMPSE2_PHASE                     } from '../../../modules/nf-core/glimpse2/phase'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_1 } from '../../../modules/nf-core/bcftools/index'
include { GLIMPSE2_LIGATE                    } from '../../../modules/nf-core/glimpse2/ligate'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2 } from '../../../modules/nf-core/bcftools/index'

workflow BAM_VCF_IMPUTE_GLIMPSE2 {

    take:
    ch_input        // channel (mandatory): [ [id], [file], [index], bamlist ]
    ch_panel        // channel (mandatory): [ [panel, chr], vcf, index ]
    ch_chunks       // channel  (optional): [ [panel, chr], region1, region2 ]
    ch_fasta        // channel (mandatory): [ [genome], fa, fai ]

    main:

    ch_versions = Channel.empty()

    // Impute with Glimpse2 without using binary files
    samples_file = Channel.of([[]]).collect()
    gmap_file    = Channel.of([[]]).collect()

    // Create input channel to impute with Glimpse2

    // Join chunks and panel
    ch_chunks_panel = ch_chunks
        .combine(ch_panel, by:0)
        .map{ metaPC, regionin, regionout, panel, index ->
            [["panel": metaPC.id, "chr": metaPC.chr], regionin, regionout, panel, index]
        }

    // Join input and chunks reference
    ch_phase_input = ch_input
        .combine(samples_file)
        .combine(ch_chunks_panel)
        .combine(gmap_file)
        .map{ metaI, bam, bai, bamlist, samples, metaPC, regionin, regionout, panel, panel_index, gmap ->
            [metaI + metaPC + ["chunk": regionout],
            bam, bai, bamlist, samples, regionin, regionout, panel, panel_index, gmap]
        }

    // Impute with Glimpse2
    GLIMPSE2_PHASE(ch_phase_input, ch_fasta)
    ch_versions = ch_versions.mix(GLIMPSE2_PHASE.out.versions.first())

    // Index phased file
    BCFTOOLS_INDEX_1(GLIMPSE2_PHASE.out.phased_variants)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_1.out.versions.first())

    // Ligate all phased files in one and index it
    ligate_input = GLIMPSE2_PHASE.out.phased_variants
        .join( BCFTOOLS_INDEX_1.out.csi )
        .map{ metaIPCR, vcf, index -> [metaIPCR.subMap("id", "panel", "chr", "batch"), vcf, index] }
        .groupTuple()

    GLIMPSE2_LIGATE(ligate_input)
    ch_versions = ch_versions.mix(GLIMPSE2_LIGATE.out.versions.first())

    BCFTOOLS_INDEX_2(GLIMPSE2_LIGATE.out.merged_variants)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_2.out.versions.first())

    // Join imputed and index files
    ch_imputed_vcf_tbi = GLIMPSE2_LIGATE.out.merged_variants
        .join(BCFTOOLS_INDEX_2.out.tbi)
        .map{ metaIPC, vcf, index -> [metaIPC + [tools: "glimpse2"], vcf, index] }

    emit:
    vcf_tbi             = ch_imputed_vcf_tbi    // channel: [ [id, panel, chr, tool], vcf, tbi ]
    versions            = ch_versions           // channel: [ versions.yml ]
}
