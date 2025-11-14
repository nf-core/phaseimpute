include { SAMTOOLS_DEPTH } from '../../../modules/nf-core/samtools/depth'
include { SAMTOOLS_VIEW  } from '../../../modules/nf-core/samtools/view'
include { GAWK           } from '../../../modules/nf-core/gawk'

workflow BAM_DOWNSAMPLE_SAMTOOLS {

    take:
    ch_bam    // channel: [ [id, genome], bam, bai ]
    ch_depth  // channel: [ [depth], depth ]
    ch_fasta  // channel: [ [genome], fasta, fai ]

    main:
    ch_versions      = Channel.empty()

    // Compute mean depth
    SAMTOOLS_DEPTH(ch_bam, [[], []])
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    // Use GAWK to get mean depth
    GAWK(SAMTOOLS_DEPTH.out.tsv, [], false)
    ch_versions = ch_versions.mix(GAWK.out.versions.first())

    // Compute downsampling factor
    ch_depth_factor = GAWK.out.output
        .splitCsv(header: false, sep:'\t')
        .map{ metaICR, row ->
            [ metaICR, row[0] as Float ]
        }
        .combine(ch_depth)
        .map{ metaICR, mean, metaD, depth ->
            [ metaICR, metaICR + metaD, depth as Float / mean ]
        }

    // Add all necessary channel for downsampling
    ch_input_downsample = ch_bam
        .combine(ch_depth_factor, by : 0)
        .map{ _metaICR, bam, index, metaICRD, depth ->
            [ metaICRD, bam, index, [], depth ]
        }

    // Downsample
    SAMTOOLS_VIEW(
        ch_input_downsample,
        ch_fasta,
        [],
        "csi"
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    // Aggregate bam and index
    ch_bam_emul = SAMTOOLS_VIEW.out.bam
        .join(SAMTOOLS_VIEW.out.csi)

    emit:
    bam_emul          = ch_bam_emul                    // channel: [ [id, chr, region, depth], bam, bai ]
    versions          = ch_versions                    // channel: [ versions.yml ]
}
