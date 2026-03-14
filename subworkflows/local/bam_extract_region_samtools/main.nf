include { SAMTOOLS_VIEW  } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_MERGE } from '../../../modules/nf-core/samtools/merge'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index'

workflow BAM_EXTRACT_REGION_SAMTOOLS {

    take:
    ch_bam    // channel: [ [id], bam, bai ]
    ch_region // channel: [ [chr, region], val(chr:start-end) ]
    ch_fasta  // channel: [ [genome], fasta, fai ]

    main:

    // Add fasta and region to bam channel
    ch_input_region = ch_bam
        .combine(ch_region)
        .map{ metaI, bam, index, metaCR, region ->
            [ metaI + metaCR + ["region_selected": region], bam, index ]
        }

    // Extract region of interest
    SAMTOOLS_VIEW(
        ch_input_region,
        [[], [], []],
        [],
        "csi"
    )

    ch_bam_region = SAMTOOLS_VIEW.out.bam
        .join(SAMTOOLS_VIEW.out.csi)

    SAMTOOLS_MERGE(
        ch_bam_region
            .map{
                metaICR, bam, index ->
                def meta_keys = metaICR.keySet() - ['chr', 'region', 'region_selected']
                [metaICR.subMap(meta_keys) + [chr: "all"], bam, index]
            }
            .groupTuple(sort: true),
        ch_fasta.map{meta, fasta, fai -> [meta, fasta, fai, []]}
    )

    SAMTOOLS_INDEX(SAMTOOLS_MERGE.out.bam)

    ch_bam_region_all = SAMTOOLS_MERGE.out.bam
        .join(SAMTOOLS_INDEX.out.bai)

    emit:
        bam_region = ch_bam_region_all // channel: [ [id, chr], bam, index ]
}
