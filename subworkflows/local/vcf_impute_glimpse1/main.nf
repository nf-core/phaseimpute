include { GLIMPSE_PHASE                      } from '../../../modules/nf-core/glimpse/phase'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_1 } from '../../../modules/nf-core/bcftools/index'
include { GLIMPSE_LIGATE                     } from '../../../modules/nf-core/glimpse/ligate'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_2 } from '../../../modules/nf-core/bcftools/index'

workflow VCF_IMPUTE_GLIMPSE1 {

    take:
    ch_input        // channel (mandatory): [ [id], vcf, tbi ]
    ch_panel        // channel (mandatory): [ [panel_id, chr], vcf, tbi ]
    ch_chunks       // channel  (optional): [ [panel_id, chr], region1, region2 ]

    main:

    ch_versions = Channel.empty()

    samples_file = Channel.of([[]]).collect()
    gmap_file    = Channel.of([[]]).collect()

    // Combine chunks with panel
    ch_chunks_panel = ch_chunks
        .combine(ch_panel, by:0)
        .map{ metaPC, regionin, regionout, panel, index ->
            [["panel_id": metaPC.panel_id, "chr": metaPC.chr], regionin, regionout, panel, index]
        }

    // Join input and chunks reference
    ch_phase_input = ch_input
        .map{ metaIPC, vcf, index -> [metaIPC.subMap("panel_id", "chr"), metaIPC, vcf, index] }
        .combine(samples_file)
        .combine(ch_chunks_panel, by: 0)
        .combine(gmap_file)
        .map{ _metaPC, metaIPC, bam, bai, samples, regionin, regionout, panel, panel_index, gmap ->
            [metaIPC + ["chunk": regionout],
            bam, bai, samples, regionin, regionout, panel, panel_index, gmap]
        }

    GLIMPSE_PHASE ( ch_phase_input ) // [meta, vcf, index, sample, regionin, regionout, ref, ref_index, map]
    ch_versions = ch_versions.mix(GLIMPSE_PHASE.out.versions )

    BCFTOOLS_INDEX_1 ( GLIMPSE_PHASE.out.phased_variants )
    ch_versions = ch_versions.mix( BCFTOOLS_INDEX_1.out.versions )

    // Ligate all phased files in one and index it
    ligate_input = GLIMPSE_PHASE.out.phased_variants
        .join( BCFTOOLS_INDEX_1.out.csi )
        .map{ metaIPCR, vcf, index -> [metaIPCR.subMap("id", "panel_id", "chr", "batch"), vcf, index] }
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
}
