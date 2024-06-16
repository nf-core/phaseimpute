include { GAWK                        } from '../../../modules/nf-core/gawk'
include { SAMTOOLS_REHEADER           } from '../../../modules/nf-core/samtools/reheader'
include { SAMTOOLS_INDEX              } from '../../../modules/nf-core/samtools/index'

workflow BAM {
    take:
    ch_bam          // channel: [ [id], bam, index ]
    ch_fasta        // channel: [ [id], fasta, fai ]

    main:

    ch_versions = Channel.empty()

    // Generate the chromosome renaming file
    GAWK(ch_fasta.map{ metaG, fasta, fai -> [metaG, fai] }, [])
    ch_versions = ch_versions.mix(GAWK.out.versions)

    // Rename the chromosome without prefix
    SAMTOOLS_REHEADER(
        ch_bam // channel: [ [id], bam, index ]
            .combine(Channel.of([[],[],[]]))
            .combine(GAWK.out.output.map{it[1]})
    )
    ch_versions = ch_versions.mix(SAMTOOLS_REHEADER.out.versions.first())

    SAMTOOLS_INDEX(SAMTOOLS_REHEADER.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_bam_renamed = SAMTOOLS_REHEADER.out.bam
        .combine(SAMTOOLS_INDEX.out.csi, by:0)

    emit:
    bam_renamed    = ch_bam_renamed        // [ [id], bam, csi ]
    versions       = ch_versions           // channel: [ versions.yml ]
}
