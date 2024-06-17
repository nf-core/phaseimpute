include { SAMTOOLS_VIEW                      } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_1 } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_MERGE                     } from '../../../modules/nf-core/samtools/merge'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_2 } from '../../../modules/nf-core/samtools/index'

workflow BAM_REGION {

    take:
    ch_bam    // channel: [ [id], bam, bai ]
    ch_region // channel: [ [chr, region], val(chr:start-end) ]
    ch_fasta  // channel: [ [genome], fasta, fai ]
    main:

    ch_versions = Channel.empty()

    // Add fasta and region to bam channel
    ch_input_region = ch_bam
        .combine(ch_region)
        .map{ metaI, bam, index, metaCR, region ->
            [ metaI + metaCR, bam, index, region, [] ]
        }

    // Extract region of interest
    SAMTOOLS_VIEW(
        ch_input_region,
        [[], []],
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    // Index region of interest
    SAMTOOLS_INDEX_1(SAMTOOLS_VIEW.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_1.out.versions.first())

    ch_bam_region = SAMTOOLS_VIEW.out.bam
        .combine(SAMTOOLS_INDEX_1.out.bai, by: 0)

    SAMTOOLS_MERGE(
        ch_bam_region
            .map{
                metaICR, bam, index -> [metaICR.subMap("id") + [chr: "all"], bam, index]
            }
            .groupTuple(),
        ch_fasta
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

    SAMTOOLS_INDEX_2(SAMTOOLS_MERGE.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_2.out.versions.first())

    ch_bam_region_all = SAMTOOLS_MERGE.out.bam
        .combine(SAMTOOLS_INDEX_2.out.bai, by:0)

    emit:
        bam_region = ch_bam_region_all // channel: [ [id, chr], bam, index ]
        versions   = ch_versions       // channel: [ versions.yml ]
}
