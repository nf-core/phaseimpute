include { SAMTOOLS_REHEADER           } from '../../../modules/nf-core/samtools/reheader'
include { SAMTOOLS_INDEX              } from '../../../modules/nf-core/samtools/index'

workflow BAM_CHR_RENAME_SAMTOOLS {
    take:
    ch_bam          // channel: [ [id], bam, index, prefix ]

    main:

    ch_versions = Channel.empty()

    // Rename the chromosome with or without prefix
    SAMTOOLS_REHEADER(
        ch_bam.map{
            meta, bam, index, prefix ->
            def cmd = ""
            if (prefix == "nochr") {
                cmd = 'sed -E "s/^(@SQ.*\\tSN:)chr/\\1/"'
            } else if (prefix == "chr") {
                cmd = 'sed -E "s/^(@SQ.*\\tSN:)([0-9]+|X|Y|MT|M)/\\1chr\\2/"'
            } else {
                error "Invalid chr_prefix: ${prefix}"
            }
            [meta, bam, index, cmd]
        }, // channel: [ [id], bam, index, cmd]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_REHEADER.out.versions.first())

    SAMTOOLS_INDEX(SAMTOOLS_REHEADER.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_bam_renamed = SAMTOOLS_REHEADER.out.bam
        .combine(SAMTOOLS_INDEX.out.bai, by:0)

    emit:
    bam_renamed    = ch_bam_renamed        // [ [id], bam, csi ]
    versions       = ch_versions           // channel: [ versions.yml ]
}
