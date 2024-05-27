include { SAMTOOLS_INDEX  } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_VIEW   } from '../../../modules/nf-core/samtools/view'

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
        ch_fasta.map{ metaG, fasta, fai -> [metaG, fasta] },
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    // Index region of interest
    SAMTOOLS_INDEX(SAMTOOLS_VIEW.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_bam_region = SAMTOOLS_VIEW.out.bam
        .combine(SAMTOOLS_INDEX.out.bai, by: 0)

    emit:
        bam_region = ch_bam_region // channel: [ [id, chr, region], bam, index ]
        versions   = ch_versions   // channel: [ versions.yml ]
}
