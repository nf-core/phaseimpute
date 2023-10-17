include { SAMTOOLS_INDEX as INDEX1     } from '../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_VIEW as VIEW_REGION } from '../../modules/nf-core/samtools/view/main.nf'

workflow BAM_REGION {

    take:
    ch_bam    // channel: [ [id, ref], bam, bai ]
    ch_region // channel: [ [ref, region], val(chr:start-end) ]
    ch_fasta  // channel: [ fasta ]
    main:

    ch_versions = Channel.empty()

    // Add fasta and region to bam channel
    ch_input_region = ch_bam
        .combine(ch_fasta)
        .combine(ch_region)
        .map{ meta, bam, index, fasta, metaR, region ->
            [meta + metaR, bam, index, fasta, region]
        }
        .combine(Channel.of([[]])) // depth parameter

    // Extract region of interest
    VIEW_REGION(ch_input_region, [])
    ch_versions = ch_versions.mix(VIEW_REGION.out.versions.first())

    // Index region of interest
    INDEX1 (VIEW_REGION.out.bam)
    ch_versions = ch_versions.mix(INDEX1.out.versions.first())

    ch_bam_region = VIEW_REGION.bam
        .combine(INDEX1.out.bai, by: 0)

    emit:
        bam_region        = ch_bam_region            // channel: [ metaIR, bam, index ]
        versions          = ch_versions              // channel: [ versions.yml ]
}