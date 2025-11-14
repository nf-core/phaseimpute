include { SAMTOOLS_VIEW  } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_MERGE } from '../../../modules/nf-core/samtools/merge'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index'

workflow BAM_EXTRACT_REGION_SAMTOOLS {

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
        [[], [], []],
        [],
        "csi"
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    ch_bam_region = SAMTOOLS_VIEW.out.bam
        .join(SAMTOOLS_VIEW.out.csi)

    SAMTOOLS_MERGE(
        ch_bam_region
            .map{
                metaICR, bam, index -> [metaICR.subMap("id", "batch") + [chr: "all"], bam, index]
            }
            .groupTuple(sort: true),
        ch_fasta
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

    SAMTOOLS_INDEX(SAMTOOLS_MERGE.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_bam_region_all = SAMTOOLS_MERGE.out.bam
        .join(SAMTOOLS_INDEX.out.bai)

    emit:
        bam_region = ch_bam_region_all // channel: [ [id, chr], bam, index ]
        versions   = ch_versions       // channel: [ versions.yml ]
}
